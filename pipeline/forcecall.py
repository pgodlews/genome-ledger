"""Stage: force-call (all-sites genotyping at fixed positions).

HaplotypeCaller emits only *variant* sites, so a marker absent from a sample's VCF is
ambiguous: genuine homozygous-reference, or simply never looked at. The single-SNP lookup
engines (traits/HLA) and the small (per-variant) PRS scores all need the genotype at
*fixed* positions **including hom-ref** — today they paper over the gap by assuming any
absent position is hom-ref (the `--absent-to-ref` convention).

This stage makes that rigorous. It force-calls the **union** of the positions those engines
care about off the CRAM:

    union of {traits.tsv, hla_tags.tsv, small PGS scoring files}  (curated, ≤ FC_MAX_LIST)
        │   build sites VCF (REF from the reference FASTA, ALT = the marker's other allele)
        ▼   GATK HaplotypeCaller  -L sites  --alleles sites.vcf  --output-mode EMIT_ALL_ACTIVE_SITES
    forcecalled/<sample>.GRCh38.vcf.gz   — explicit GT + DP at every requested site

An *assumed* hom-ref becomes an *evidenced* `0/0`; a position the reads can't resolve
becomes an explicit no-call (`./.`, low coverage) instead of being silently scored as
reference. The genome-wide PRS demonstrator (millions of variants) stays on the plink2
path — force-calling it would be heavy and is largely redundant.

The sites VCF is sample-independent and cached (hashed by its inputs) under FC_SITES_DIR,
so only the per-sample GATK pass repeats. The stage is idempotent (skips when its output
exists) and reads off the CRAM, like `call`.
"""

from __future__ import annotations

import hashlib
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from . import reportpaths
from . import config
from .call import _java_env, aligned_cram
from .util import (faidx_bases, file_token, log, read_tsv, require_tools, run,
                   write_report)

# Acceptable single-base alleles (SNV markers only — indels in the curated lists are rare
# and force-calling them reliably needs the original ALT representation, not a faidx base).
_ACGT = frozenset("ACGT")


# --- paths -------------------------------------------------------------------
def forcecalled_vcf(sample: str) -> Path:
    return config.FORCECALL_DIR / f"{sample}.{config.BUILD}.vcf.gz"


def available(sample: str) -> bool:
    """True when a force-called VCF exists for the sample (consumers prefer it)."""
    return forcecalled_vcf(sample).exists()


# --- site collection ---------------------------------------------------------
def _tsv_sites(path: Path) -> list[tuple[str, str, str, str]]:
    """(chrom, pos, ref, effect) rows from a traits/HLA rules TSV (shared schema)."""
    return [(r["chrom"], r["pos"], r["ref"], r["effect"]) for r in read_tsv(path)]


def _prs_sites() -> list[tuple[str, str, str, str]]:
    """(chrom, pos, effect, other) from every PGS scoring file small enough to force-call.

    Sourced from `prs._panel_scores()` — the SAME list prs_report scores. It used to read
    config.PRS_SCORES, which is the 3-score fallback list, while the report had long since
    moved to the 44-score COMPLEX_PANEL. The two share no PGS ids at all, so the force-call
    union covered essentially none of the markers the report actually scores: measured
    at 65 of 35,175 small-score markers. Combined with prs_report substituting
    the force-called VCF for the called one, that turned every small score for a force-called
    sample into an assumed-reference genome. The substitution is fixed separately; this
    closes the gap that made it so total.

    The per-score size test is unchanged and is the same FC_MAX_LIST == PRS_PLINK_MIN
    threshold prs uses to route a score to plink2, so the two stay in step by construction.

    Imported lazily — prs pulls in network/plink2 helpers we don't otherwise need here,
    and a missing scoring file should not break force-calling the curated TSV lists."""
    from . import prs as prs_mod

    out: list[tuple[str, str, str, str]] = []
    for pid in [r["pgs_id"] for r in prs_mod._panel_scores()]:
        try:
            variants = prs_mod._load_variants(prs_mod._scoring_file(pid))
        except Exception as e:  # noqa: BLE001 — best-effort; offline runs skip PRS sites
            log.warning("force-call: skipping PGS %s (%s)", pid, e)
            continue
        if len(variants) > config.FC_MAX_LIST:
            log.info("force-call: PGS %s is genome-wide (%d variants) — left to plink2.",
                     pid, len(variants))
            continue
        out += [(c, p, ea, oa) for c, p, ea, oa, _w in variants]
    return out


def _collect_sites() -> dict[tuple[str, str], set[str]]:
    """(chrom, pos) → the set of alleles the marker can take, across all curated lists.

    Allele sets are unioned so one force-called record serves every engine that wants the
    position. REF is resolved later from the reference FASTA (authoritative), so the two
    alleles a list supplies just seed which ALT to force."""
    rows = (_tsv_sites(config.TRAITS_FILE) + _tsv_sites(config.HLA_FILE) + _prs_sites())
    sites: dict[tuple[str, str], set[str]] = {}
    for chrom, pos, a, b in rows:
        if not (chrom and pos):
            continue
        sites.setdefault((chrom, str(pos)), set()).update(
            x.upper() for x in (a, b) if x)
    return sites


# --- sites VCF (cached, sample-independent) ----------------------------------
def _ref_bases(positions: list[tuple[str, str]]) -> dict[tuple[str, str], str]:
    """Reference base at each position via samtools faidx (the authoritative REF).
    Batched (see util.faidx_bases) — FC_MAX_LIST allows up to 50k PGS positions into the
    force-call union, which would otherwise blow the OS exec() argument-list limit in one
    unbatched `samtools faidx` call (a failure already hit once with a large marker panel)."""
    return faidx_bases(positions)


def _fai_lengths() -> dict[str, int]:
    fai = config.REF_FASTA.with_suffix(".fa.fai")
    out: dict[str, int] = {}
    for line in fai.read_text().splitlines():
        name, length = line.split("\t")[:2]
        out[name] = int(length)
    return out


def _sites_records(sites) -> list[tuple[str, int, str, str]]:
    """(chrom, pos, ref, alt) force-call records, coordinate-sorted.

    REF is the reference base; ALT is the marker's other allele(s) — those supplied by a
    list that differ from REF. SNV-only and ACGT-only; a position whose alleles are all
    reference (no derivable ALT) is dropped from --alleles (the reads still get a shot at
    it via -L/EMIT_ALL_ACTIVE_SITES, and the consumer falls back to hom-ref)."""
    refbases = _ref_bases(list(sites))
    recs: list[tuple[str, int, str, str]] = []
    for (chrom, pos), alleles in sites.items():
        ref = refbases.get((chrom, pos))
        if not ref or ref not in _ACGT:
            continue
        alts = sorted(a for a in alleles if a in _ACGT and a != ref)
        if not alts:
            continue
        recs.append((chrom, int(pos), ref, ",".join(alts)))
    # Sort in reference (.fai) order so the sites VCF/intervals are sequence-dict-concordant.
    order = {c: i for i, c in enumerate(_fai_lengths())}
    recs.sort(key=lambda r: (order.get(r[0], len(order)), r[0], r[1]))
    return recs


def _sites_hash(recs) -> str:
    h = hashlib.sha256()
    for c, p, ref, alt in recs:
        h.update(f"{c}\t{p}\t{ref}\t{alt}\n".encode())
    return h.hexdigest()[:12]


def build_sites(force: bool = False) -> tuple[Path, Path]:
    """Build (cached) the force-call sites VCF + interval list. Returns (vcf, intervals).

    Cache key is a hash of the records, so editing a rules file or adding a PGS score
    transparently produces a new cache entry (and forces a re-call downstream)."""
    require_tools("bgzip", "tabix", "samtools")  # _ref_bases shells out to samtools faidx
    sites = _collect_sites()
    recs = _sites_records(sites)
    if not recs:
        raise SystemExit("force-call: no SNV sites collected from the curated lists.")
    key = _sites_hash(recs)
    config.FC_SITES_DIR.mkdir(parents=True, exist_ok=True)
    vcf = config.FC_SITES_DIR / f"sites_{key}.vcf.gz"
    intervals = config.FC_SITES_DIR / f"sites_{key}.list"
    if vcf.exists() and intervals.exists() and not force:
        return vcf, intervals

    lengths = _fai_lengths()
    used_contigs = list(dict.fromkeys(c for c, *_ in recs))  # recs already in .fai order
    header = ["##fileformat=VCFv4.2",
              *[f"##contig=<ID={c},length={lengths[c]}>" for c in used_contigs],
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    body = [f"{c}\t{p}\t.\t{ref}\t{alt}\t.\t.\t." for c, p, ref, alt in recs]
    raw = config.FC_SITES_DIR / f"sites_{key}.vcf"
    raw.write_text("\n".join(header + body) + "\n")
    run(["bgzip", "-f", str(raw)])              # → sites_<key>.vcf.gz
    run(["tabix", "-f", "-p", "vcf", str(vcf)])
    intervals.write_text("".join(f"{c}:{p}-{p}\n" for c, p, *_ in recs))
    log.info("force-call sites: %d SNV positions (%d contigs) → %s",
             len(recs), len(used_contigs), vcf.name)
    return vcf, intervals


# --- force-call --------------------------------------------------------------
def force_call(sample: str, force: bool = False) -> Path:
    """GATK force-call the curated union of positions off the CRAM → all-sites VCF."""
    require_tools("bcftools", "tabix", "bgzip", "samtools")
    if not config.GATK_JAR.exists():
        raise SystemExit(f"GATK not installed at {config.GATK_JAR}. Run `setup`.")
    cram = aligned_cram(sample)
    if not cram.exists():
        raise SystemExit(f"{cram} missing. Run `align {sample}` first.")
    out = forcecalled_vcf(sample)
    out.parent.mkdir(parents=True, exist_ok=True)

    vcf_sites, intervals = build_sites(force=force)
    # The output is only valid for the site set it was called against. The sites VCF is
    # content-hashed (sites_<hash>), so key the per-sample output to that hash — editing a
    # rules file or adding a PGS score must force a re-call, not silently reuse stale GTs.
    # Key on the CRAM as well as the site set. The output is a function of both, but only the
    # sites were keyed, so a re-aligned sample silently kept the genotypes called off its
    # previous alignment — the site set had not changed, so nothing looked stale.
    sites_key = (vcf_sites.name.removeprefix("sites_").removesuffix(".vcf.gz")
                 + "|cram=" + file_token(cram))
    key_file = config.FORCECALL_DIR / f"{sample}.{config.BUILD}.sites"
    tbi = out.with_suffix(".gz.tbi")
    if out.exists() and tbi.exists() and not force \
            and key_file.exists() and key_file.read_text().strip() == sites_key:
        log.info("%s already force-called (%s) — skipping (use --force).", sample, out.name)
        return out
    if out.exists():
        log.info("%s force-call site set changed (or --force) — re-calling.", sample)

    call_at_sites(sample, cram, vcf_sites, intervals, out, sites_key, tag="forcecall",
                  force=force)
    key_file.write_text(sites_key + "\n")
    _qc_report(sample, out)
    log.info("Wrote %s. Consumers (`traits`, `prs`) now prefer it over absent→ref.", out)
    return out


def call_at_sites(sample: str, cram: Path, vcf_sites: Path, intervals: Path, out: Path,
                  sites_key: str, tag: str = "forcecall", force: bool = False) -> Path:
    """GATK-genotype one sample at a fixed site list off its CRAM → all-sites VCF at `out`.

    The site-set-agnostic core of this stage, shared with `denovo` (which force-calls the
    parents at the child's candidate positions — a site list built per run, not the curated
    union). `sites_key` keys the resume shards so a changed site set can never be served
    from a stale part file; `tag` separates the two stages' scratch dirs and logs."""
    env = _java_env()
    parts_dir = out.parent / f".parts_{tag}_{sample}"
    parts_dir.mkdir(parents=True, exist_ok=True)
    # Resume shards are only valid for the same site set — a changed key invalidates them.
    # --force must actually re-call. The shards are keyed by sites_key, so forcing with an
    # unchanged site set found every part present and returned them untouched: the CRAM was
    # never re-read and --force did nothing at all.
    parts_key = parts_dir / ".sites_key"
    if force or not parts_key.exists() or parts_key.read_text().strip() != sites_key:
        for stale in parts_dir.glob("*"):
            stale.unlink()
    parts_key.write_text(sites_key + "\n")

    # Scatter by contig so sites across the genome call in parallel; each shard restricts
    # -L to that contig's slice of the (sorted) interval list. Resume reuses parts.
    by_contig = _intervals_by_contig(intervals)
    log.info("%s: HaplotypeCaller force-call for %s over %d contigs, %d in parallel.",
             tag, sample, len(by_contig), config.HC_WORKERS)

    def _shard(item: tuple[str, list[str]]) -> Path:
        contig, ivals = item
        part = parts_dir / f"{sample}.{contig}.vcf.gz"
        if part.exists() and part.with_suffix(".gz.tbi").exists():
            return part
        ilist = parts_dir / f"{sample}.{contig}.list"
        ilist.write_text("".join(f"{i}\n" for i in ivals))
        run(["java", f"-Xmx{config.HC_JAVA_MEM}", "-jar", config.GATK_JAR,
             "HaplotypeCaller", "-R", config.REF_FASTA, "-I", cram,
             "-L", ilist, "--alleles", vcf_sites,
             "--output-mode", "EMIT_ALL_ACTIVE_SITES",
             "--sample-ploidy", "2", "--native-pair-hmm-threads", "1", "-O", part],
            env=env,
            stdout=open(config.LOGS_DIR / f"{tag}_{sample}_{contig}.log", "wb"),
            stderr=subprocess.STDOUT)
        return part

    with ThreadPoolExecutor(max_workers=config.HC_WORKERS) as ex:
        parts = list(ex.map(_shard, by_contig))

    try:
        run(["bcftools", "concat", "-a", "-Oz", "-o", str(out), *[str(p) for p in parts]])
        run(["tabix", "-f", "-p", "vcf", str(out)])
    except Exception:
        out.unlink(missing_ok=True)  # don't leave a truncated VCF to be skipped-to later
        raise
    for p in parts_dir.glob("*"):
        p.unlink()
    parts_dir.rmdir()
    return out


def _intervals_by_contig(intervals: Path) -> list[tuple[str, list[str]]]:
    """Group the sorted interval list (`chrom:pos-pos` lines) by contig, in file order."""
    groups: dict[str, list[str]] = {}
    for ln in intervals.read_text().splitlines():
        if ln:
            groups.setdefault(ln.split(":")[0], []).append(ln)
    return list(groups.items())


# --- consumer API ------------------------------------------------------------
def _decode(ref: str, alt: str, gt: str, dp: str) -> dict:
    """Turn one force-called record's (REF, ALT, GT, DP) into a genotype dict.

    A `./.` GT is a genuine no-call (the reads couldn't resolve the site) — kept distinct
    from an evidenced `0/0` hom-ref. A PARTIAL no-call (`./1`) is also uncertainty, not a
    confident hom-alt. A haploid record (unexpected for these autosomal/chr6 sites) is
    treated as homozygous."""
    coding = [ref, *alt.split(",")] if alt not in (".", "") else [ref]
    idx = gt.replace("|", "/").split("/")
    if any(g in (".", "") for g in idx):
        alleles = []  # (partial) no-call — uncertainty, not evidence of any genotype
    else:
        try:
            alleles = [coding[int(g)] for g in idx]
        except (ValueError, IndexError):
            alleles = []
    nocall = len(alleles) == 0
    if len(alleles) == 1:            # haploid record → treat as homozygous
        alleles *= 2
    return {"ref": ref, "alt": alt, "gt": gt,
            "alleles": None if nocall else alleles,
            "dp": int(dp) if dp.isdigit() else 0, "nocall": nocall}


def read_genotypes(sample: str) -> dict[tuple[str, str], dict]:
    """(chrom, pos) → {ref, alt, gt, alleles, dp, nocall} from the force-called VCF.

    `alleles` is the pair of called bases (None when no-call). Returns {} when the sample
    hasn't been force-called, so callers degrade to their old absent→hom-ref path."""
    vcf = forcecalled_vcf(sample)
    if not vcf.exists():
        return {}
    # check=True: a failed query (truncated VCF, missing bcftools) must raise, not return
    # an empty dict — callers read {} as "not force-called" and fall back to assuming
    # hom-ref, the exact silent behavior this stage exists to eliminate.
    out = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%DP]\n", str(vcf)],
        capture_output=True, text=True, check=True).stdout
    geno: dict[tuple[str, str], dict] = {}
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 6:
            continue
        chrom, pos, ref, alt, gt, dp = f[:6]
        geno[(chrom, pos)] = _decode(ref, alt, gt, dp)
    return geno


# --- QC ----------------------------------------------------------------------
def _qc_report(sample: str, vcf: Path) -> Path:
    geno = read_genotypes(sample)
    n = len(geno)
    nocall = sum(1 for g in geno.values() if g["nocall"])
    called = [g for g in geno.values() if not g["nocall"]]
    homref = sum(1 for g in called if set(g["alleles"]) == {g["ref"]})
    variant = len(called) - homref
    dps = sorted(g["dp"] for g in called)
    med = dps[len(dps) // 2] if dps else 0
    low = sum(1 for d in dps if d < 8)
    out_dir = config.FORCECALL_DIR / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    pct = (100 * len(called) / n) if n else 0.0
    lines = [
        f"# Force-call (all-sites genotyping) — {sample}",
        "",
        f"- **Sites genotyped:** {n}  ·  source: curated single-SNP + small-PRS union",
        f"- **Called:** {len(called)} ({pct:.1f}%)  —  evidenced hom-ref {homref}, "
        f"variant {variant}",
        f"- **No-calls (`./.`, low/zero coverage):** {nocall}",
        f"- **Median depth at called sites:** {med}×  ·  sites with DP<8: {low}",
        "",
        "## What this sharpens",
        "- An *assumed* hom-ref is now an *evidenced* `0/0`; the `traits`/`hla` and small-PRS "
        "engines read genotypes from this VCF when present.",
        "- No-call sites are reported as indeterminate (or dropped from a score) instead of "
        "being silently counted as homozygous reference.",
        "- Genome-wide PRS is unaffected — it stays on the plink2 path by design.",
        "",
    ]
    summary = reportpaths.genome_report("forcecall", sample)
    write_report(summary, lines)
    log.info("Force-call QC: %d sites, %d no-call, median DP %d× (%s)",
             n, nocall, med, summary.name)
    return summary
