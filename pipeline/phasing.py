"""Stage: read-backed phasing (WhatsHap) to resolve in-cis vs in-trans.

The monogenic panel flags "possible compound het — review" when a person carries two
different heterozygous alleles in one recessive gene — but a SNV VCF cannot tell whether
they sit on the *same* chromosome copy (**in-cis**: one copy is still intact) or on
*opposite* copies (**in-trans**: both copies hit, i.e. potentially biallelic loss). Only the
reads resolve it. This stage finds disease-panel genes where the person carries ≥2
**predicted-damaging** heterozygous variants (LOF or damaging-missense), phases the reads
around them with WhatsHap, and reports whether those variants are in-cis or in-trans.

Limits (documented in the report):
  - Read-backed phasing only connects variants a read/read-pair spans (~hundreds of bp to
    ~1 kb); damaging variants farther apart land in separate phase blocks and read as
    *unphased* — long reads or trio phasing would extend the range.
  - This is the single-sample read-backed mode. With a sequenced child, the `segregation`
    stage resolves cis/trans by transmission at any distance.

This is the whole-family *moat*: a product without your reads cannot resolve cis/trans.
"""

from __future__ import annotations

import gzip
import subprocess
from pathlib import Path

from . import reportpaths
from . import config
from .call import aligned_cram
from .load import TABLE, connect
from .normalize import normalized_path
from .panel import load_panel, sql_in_list
from .snapshot import latest_snapshot
from .util import log, write_report

_PAD = 1000          # bp window each side of a damaging variant (the read-spanning range)
# A heterozygous variant counts as "potentially damaging" if LOF or damaging-missense.
_DAMAGING_SQL = (
    "(consequence LIKE '%stop_gained%' OR consequence LIKE '%frameshift%' "
    "OR consequence LIKE '%splice_acceptor%' OR consequence LIKE '%splice_donor%' "
    "OR consequence LIKE '%start_lost%' OR (consequence LIKE '%missense%' "
    "AND (sift LIKE 'deleterious%' OR polyphen LIKE 'probably_damaging%')))")


def _short_csq(csq: str) -> str:
    for t in ("stop_gained", "frameshift", "splice_acceptor", "splice_donor", "start_lost",
              "missense"):
        if t in (csq or ""):
            return t
    return (csq or "").split(",")[0]


def _candidate_genes(con, snapshot_id: str, sample: str, panel_genes) -> list[dict]:
    """Panel genes where the sample has ≥2 damaging het variants, with those variants."""
    gl = sql_in_list(panel_genes)
    rows = con.execute(f"""
        SELECT gene, chrom, pos, ref, alt, consequence, clinvar_sig
        FROM {TABLE}
        WHERE snapshot_id = ? AND sample_id = ? AND gene IN ({gl})
          AND zygosity = 'HET' AND {_DAMAGING_SQL}
        ORDER BY gene, pos
    """, [snapshot_id, sample]).fetchall()
    by_gene: dict[tuple[str, str], list] = {}
    for g, c, p, ref, alt, csq, sig in rows:
        by_gene.setdefault((g, c), []).append(
            {"pos": p, "ref": ref, "alt": alt, "csq": _short_csq(csq), "sig": sig})
    return [{"gene": g, "chrom": c, "variants": vs}
            for (g, c), vs in by_gene.items() if len(vs) >= 2]


def _windows_text(candidates: list[dict]) -> str:
    """Merged ±_PAD windows around every candidate variant, as BED text (0-based)."""
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for cg in candidates:
        for v in cg["variants"]:
            by_chrom.setdefault(cg["chrom"], []).append(
                (max(1, v["pos"] - _PAD), v["pos"] + _PAD))
    out = []
    for chrom, spans in by_chrom.items():
        spans.sort()
        s0, e0 = spans[0]
        for s, e in spans[1:]:
            if s <= e0:
                e0 = max(e0, e)
            else:
                out.append((chrom, s0, e0))
                s0, e0 = s, e
        out.append((chrom, s0, e0))
    return "".join(f"{c}\t{s - 1}\t{e}\n" for c, s, e in out)


def _run_whatshap(sample: str, candidates: list[dict], out_dir: Path,
                  force: bool = False) -> Path:
    phased = out_dir / f"{sample}.phased.vcf.gz"
    bed = out_dir / "windows.bed"
    content = _windows_text(candidates)
    # The phased VCF is only valid for the candidate windows it was built from. Reuse it
    # only when those windows are unchanged and the caller isn't forcing — a new snapshot,
    # an extended panel, or a re-normalized VCF must re-phase, not silently read stale
    # phase (new variants would misread as "unphased").
    if phased.exists() and not force and bed.exists() and bed.read_text() == content:
        log.info("Reusing existing phased VCF %s (candidates unchanged).", phased.name)
        return phased
    if phased.exists():
        log.info("Phasing candidates changed (or --force) — re-phasing %s.", sample)
        phased.unlink()
    cram = aligned_cram(sample)
    if not cram.exists():
        raise SystemExit(f"{cram} missing. Run `align {sample}` first.")
    bed.write_text(content)
    sub = out_dir / f"{sample}.regions.vcf.gz"
    # Phase the NORMALIZED VCF — the same representation DuckDB is built from, so the
    # damaging-variant positions line up (modulo the indel anchor handled in _parse_phased).
    subprocess.run(["bcftools", "view", "-R", str(bed), "-Oz", "-o", str(sub),
                    str(normalized_path(sample))], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(sub)], check=True)
    log.info("Phasing %d candidate gene(s) for %s with WhatsHap …", len(candidates), sample)
    with open(config.LOGS_DIR / f"phase_{sample}.log", "wb") as logf:  # no fd leak
        subprocess.run(
            ["uv", "run", "--with", config.WHATSHAP_SPEC, "whatshap", "phase",
             "--reference", str(config.REF_FASTA), "--ignore-read-groups",
             "-o", str(phased), str(sub), str(cram)],
            check=True, stdout=logf, stderr=subprocess.STDOUT)
    return phased


def _parse_phased(phased: Path) -> dict[tuple[str, str], dict]:
    """(chrom,pos) → {ps, alt_hap} for heterozygous variants; alt_hap is the haplotype side
    carrying ALT (0 for '1|0', 1 for '0|1'); both None when the variant wasn't phased.

    Keyed by position alone, because `_verdict` looks these up with DuckDB's VEP-style
    coordinates and the two sides disagree on indel alleles (hence the +1 anchor entry at
    the bottom). That makes the key ambiguous for a split multiallelic — `bcftools norm
    -m -both` turns one 1/2 genotype into two het records at the same POS — and plain
    last-wins assignment then handed both records whichever haplotype came second. A
    genuine compound het with the two alleles on opposite copies was reported as "in-cis,
    other copy intact": the precise inversion of the call this module exists to make.

    Conflicting records at one position are therefore marked ambiguous and read as unphased
    downstream. Refusing to phase a position the key cannot disambiguate is the honest
    answer, and it fails toward "unknown" rather than toward false reassurance.
    """
    out: dict[tuple[str, str], dict] = {}

    def place(key, rec: dict, src: tuple, anchor: bool = False) -> None:
        prior = out.get(key)
        if prior is None:
            out[key] = {**rec, "src": src, "anchor": anchor}
            return
        if anchor:
            return                       # an anchor entry never displaces a real record
        if prior["anchor"]:
            out[key] = {**rec, "src": src, "anchor": False}   # ...but a real record displaces it
            return
        if prior["src"] == src:
            return                       # the same record seen twice
        if prior["ps"] != rec["ps"] or prior["alt_hap"] != rec["alt_hap"]:
            out[key] = {"ps": None, "alt_hap": None, "src": None,
                        "anchor": False, "ambiguous": True}   # sticky: never un-marked

    opener = gzip.open if str(phased).endswith(".gz") else open
    with opener(phased, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, pos, ref, alt = f[0], f[1], f[3], f[4]
            fmt, sample = f[8].split(":"), f[9].split(":")
            try:
                gt = sample[fmt.index("GT")]
            except (ValueError, IndexError):
                continue
            a = gt.replace("|", "/").split("/")
            if len(a) != 2 or a[0] == a[1]:
                continue
            # Guarded like GT above: "PS" in fmt only proves the FORMAT key exists, not that
            # the sample column has that many colon-separated cells — a truncated record
            # raises IndexError here otherwise.
            ps = None
            if "PS" in fmt and "|" in gt:
                pi = fmt.index("PS")
                ps = sample[pi] if pi < len(sample) else None
            alt_hap = (0 if gt.split("|")[0] != "0" else 1) if "|" in gt else None
            rec = {"ps": ps, "alt_hap": alt_hap, "ambiguous": False}
            place((chrom, pos), rec, (ref, alt))
            if len(ref) != len(alt):       # indel: VEP/DuckDB anchor one base to the right
                place((chrom, str(int(pos) + 1)), rec, (ref, alt), anchor=True)
    return out


def _verdict(cg: dict, phased: dict) -> dict:
    """Resolve the haplotype configuration of a gene's damaging het variants."""
    chrom = cg["chrom"]
    placed = []          # (pos, ps, alt_hap, variant)
    n_ambiguous = 0
    for v in cg["variants"]:
        rec = phased.get((chrom, str(v["pos"])))
        if rec is None:
            placed.append((v["pos"], None, None, v))           # not in phased output
        elif rec.get("ambiguous"):
            n_ambiguous += 1                                   # two ALTs share this key
            placed.append((v["pos"], None, None, v))
        else:
            placed.append((v["pos"], rec["ps"], rec["alt_hap"], v))
    blocks: dict[str, list] = {}
    n_phased = 0
    for _pos, ps, hap, _v in placed:
        if ps is not None and hap is not None:
            n_phased += 1
            blocks.setdefault(ps, []).append(hap)
    # cis/trans can only be called when ≥2 damaging variants land in the SAME phase block.
    in_trans = any(len(set(h)) == 2 for h in blocks.values())
    one_sided = any(len(h) >= 2 and len(set(h)) == 1 for h in blocks.values())
    # "the other copy is intact" is only supportable when EVERY damaging variant in the gene
    # was placed, in one block, on one haplotype. The old test asked merely whether *some*
    # block held ≥2 variants on one side, so a third damaging variant left unphased —
    # potentially sitting on the other copy — still produced the reassuring verdict.
    fully_resolved = n_phased == len(placed) and len(blocks) == 1
    resolved_cis = one_sided and fully_resolved
    if in_trans:
        verdict = "⚠️ in-trans (both copies hit — review)"
    elif resolved_cis:
        verdict = "in-cis (≥2 variants on one copy; other intact)"
    elif one_sided:
        verdict = (f"partially resolved (≥2 variants share one copy; "
                   f"{len(placed) - n_phased} unplaced — other copy NOT established)")
    else:
        verdict = "unphased (variants in separate blocks — relative phase unknown)"
    if n_ambiguous:
        verdict += (f" · {n_ambiguous} multiallelic position(s) not separable by coordinate")
    return {"gene": cg["gene"], "chrom": chrom, "n": len(cg["variants"]),
            "n_phased": n_phased, "blocks": len(blocks), "in_trans": in_trans,
            "n_ambiguous": n_ambiguous, "verdict": verdict, "placed": placed}


def _hap_label(ps, hap) -> str:
    if ps is None or hap is None:
        return "unphased"
    return f"hap-{'A' if hap == 0 else 'B'} (block {ps})"


def phase_report(sample: str, force: bool = False) -> Path:
    snapshot_id = latest_snapshot()
    if not snapshot_id:
        raise SystemExit("No snapshot loaded — run `scan` first.")
    con = connect()
    panel_genes = {g for d in load_panel() for g in d.genes}
    out_dir = config.PHASE_DIR / sample
    summary = reportpaths.genome_report("phase", sample)
    if summary.exists() and not force:
        log.info("%s phasing already reported (%s) — skipping (use --force).", sample, summary)
        con.close()
        return summary
    out_dir.mkdir(parents=True, exist_ok=True)

    candidates = _candidate_genes(con, snapshot_id, sample, panel_genes)
    con.close()
    results = []
    if candidates:
        phased = _parse_phased(_run_whatshap(sample, candidates, out_dir, force=force))
        results = [_verdict(cg, phased) for cg in candidates]
    results.sort(key=lambda r: (not r["in_trans"], r["gene"]))
    trans = [r for r in results if r["in_trans"]]

    detail = []
    for r in results:
        detail.append(f"\n### {r['gene']} — {r['verdict']}\n")
        detail.append("| Variant | Consequence | ClinVar | Haplotype |")
        detail.append("|---------|-------------|---------|-----------|")
        for pos, ps, hap, v in sorted(r["placed"], key=lambda x: x[0]):
            detail.append(f"| {r['chrom']}:{pos} {v['ref']}>{v['alt']} | {v['csq']} | "
                          f"{v['sig'] or '—'} | {_hap_label(ps, hap)} |")

    lines = [
        f"# Read-backed phasing (cis/trans) — {sample}",
        "",
        f"- **Snapshot:** `{snapshot_id}`  ·  **Tool:** WhatsHap (read-backed, single-sample)",
        f"- **Panel genes with ≥2 damaging hets:** {len(results)}",
        f"- **Resolved in-trans (both copies hit — review):** {len(trans)}",
        "",
        "## Summary",
        "",
        "| Gene | Damaging hets | Phased | Blocks | Verdict |",
        "|------|--------------:|-------:|-------:|---------|",
        *([f"| {r['gene']} | {r['n']} | {r['n_phased']} | {r['blocks']} | {r['verdict']} |"
           for r in results] or ["| _none_ | | | | |"]),
        "",
        "## Per-gene haplotypes",
        *detail,
        "",
        "## Method & caveats",
        "- **In-trans** means the two damaging variants sit on opposite chromosome copies — "
        "in a recessive gene that potentially disables *both* copies (review with the panel / "
        "ACMG reports). **In-cis** means they share one copy, leaving the other intact; it is "
        "claimed only when *every* damaging variant in the gene was placed in a single phase "
        "block. When some remain unplaced the verdict is **partially resolved**: the shared "
        "copy is established, the other one is not.",
        "- A **multiallelic position** (two different damaging ALT alleles at one coordinate) "
        "cannot be told apart by the position-keyed lookup this report uses, so it is left "
        "unphased rather than guessed at.",
        "- **Configuration ≠ pathogenicity:** \"damaging\" here is LOF or "
        "deleterious-missense (SIFT/PolyPhen); confirm the variants are truly damaging.",
        "- **Read-backed range** is ~hundreds of bp to ~1 kb; variants farther apart read as "
        "*unphased* (separate blocks), which is **not** evidence of in-cis. Long reads or "
        "trio phasing would connect them.",
        "- **With a sequenced child**, the `segregation` report resolves cis/trans by "
        "transmission at any distance; this report is read-backed only.",
        "- Research-grade, not a clinical diagnosis.",
        "",
    ]
    write_report(summary, lines)
    log.info("Phasing report: %s (%d candidate genes, %d in-trans)",
             summary, len(results), len(trans))
    return summary
