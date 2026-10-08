"""Stage: trio segregation + transmission-phased compound heterozygosity.

The monogenic panel flags "possible compound het — review" whenever someone carries two
different damaging heterozygous variants in one recessive gene, because a SNV VCF cannot say
whether they sit on the *same* chromosome copy (**in-cis** — one copy still intact) or on
opposite copies (**in-trans** — both copies hit). `phase` (WhatsHap) answers that from reads,
but only for variants close enough for a read or read-pair to span: a few hundred bases to
about a kilobase. Most compound-het candidates are further apart than that and come back
honestly "unphased".

A trio answers it at any distance and without the child's reads at all. Each of the child's
chromosome copies came from one parent, so:

    variant A inherited from the father, variant B from the mother   → in-trans
    both A and B inherited from the same parent                      → in-cis

The catch is the same one the `denovo` stage exists to handle, and it bites harder here: the
parents' VCFs are **variant-only**, so "the mother does not carry A" might mean she is hom-ref,
or that nobody looked. Getting that wrong flips the verdict — an in-cis pair (reassuring, one
good copy) reads as in-trans (both copies disabled) or the reverse. A child's
apparently-private variants are often carried by a parent whose VCF missed them, so this is
a common error, not a corner case.

So every "parent does not carry it" that a verdict depends on is confirmed against the
parents' **reads** — GATK force-call plus a caller-independent pileup, the same two-pass check
`denovo` uses — and any site the reads cannot resolve makes the verdict `unresolved` rather
than guessing. Runs off the parents' CRAMs; the child needs only a VCF.
"""

from __future__ import annotations

from pathlib import Path

from . import reportpaths
from . import config, pedigree
from .call import aligned_cram
from .denovo import (CONFIRMED, INHERITED, _pileup_parent, _read_parent, _write_sites,
                     parent_callsets, parent_verdict)
from .forcecall import call_at_sites
from .load import TABLE, connect
from .panel import load_panel, sql_in_list
from .phasing import _DAMAGING_SQL, _short_csq
from .snapshot import latest_snapshot
from .util import log, write_report

# Transmission of one variant the child carries.
PATERNAL = "paternal"
MATERNAL = "maternal"
BOTH = "both parents"
DE_NOVO = "neither parent (de-novo candidate)"
UNRESOLVED = "unresolved"

IN_TRANS = "in-trans — both gene copies hit"
IN_CIS = "in-cis — one copy intact"


def report_path(child: str) -> Path:
    return reportpaths.genome_report("segregation", child)


def _work_dir(child: str) -> Path:
    return config.SEGREGATION_DIR / child


# --- data -------------------------------------------------------------------
def _child_candidates(con, snap: str, child: str, panel_genes) -> list[dict]:
    """Panel genes where the child carries ≥2 damaging het variants, with those variants.

    Same definition `phase` uses, so the two stages answer the same question about the same
    variant pairs — one from reads, one from transmission."""
    gl = sql_in_list(panel_genes)
    rows = con.execute(f"""
        SELECT gene, chrom, pos, ref, alt, consequence, clinvar_sig
        FROM {TABLE}
        WHERE snapshot_id = ? AND sample_id = ? AND gene IN ({gl})
          AND zygosity = 'HET' AND {_DAMAGING_SQL}
        ORDER BY gene, pos
    """, [snap, child]).fetchall()
    by_gene: dict[str, list[dict]] = {}
    for gene, chrom, pos, ref, alt, csq, sig in rows:
        by_gene.setdefault(gene, []).append(
            {"gene": gene, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
             "csq": _short_csq(csq), "sig": sig})
    return [v for vs in by_gene.values() if len(vs) >= 2 for v in vs]


def _parent_calls(con, snap: str, callsets: list[str], variants) -> set[tuple]:
    """{(callset, chrom, pos, ref, alt)} the parents' callsets record as carried.

    Every callset ingested for a parent counts (`parent` *and* `parent_t2t`): one person
    called twice is one genome, and each caller's false negatives are largely its own."""
    if not variants or not callsets:
        return set()
    chroms = sql_in_list({v["chrom"] for v in variants})
    positions = ", ".join(str(pos) for pos in sorted({v["pos"] for v in variants}))
    rows = con.execute(f"""
        SELECT sample_id, chrom, pos, ref, alt FROM {TABLE}
        WHERE snapshot_id = ? AND sample_id IN ({sql_in_list(callsets)})
          AND chrom IN ({chroms}) AND pos IN ({positions})
    """, [snap]).fetchall()
    return {tuple(r) for r in rows}


def _same_event(v: dict, ref: str, alt: str) -> bool:
    """Does a VCF record's (ref, alt) spell the same event as VEP's (v['ref'], v['alt'])?

    VEP drops the shared anchor base and writes '-' for the emptied side, so the VCF form is
    recoverable exactly — and anything short of an exact reconstruction is a *different*
    variant. Matching on "is an indel" alone (`len(ref) != len(alt)`, which this used to do)
    binds the first indel record at the position: `normalize` runs `bcftools norm -m -both`,
    which splits multiallelics into separate records at one POS, and ~14% of indel records in
    a real normalized VCF share a position with another one. That coin-flip then force-calls
    the parents at the wrong allele and the cis/trans verdict rests on it, so a miss here is
    worse than the unresolved-and-dropped path a non-match takes."""
    vr, va = v["ref"].upper(), v["alt"].upper()
    ref, alt = ref.upper(), alt.upper()
    if va == "-":                     # deletion: VCF is <anchor><deleted> -> <anchor>
        return alt == ref[:1] and ref[1:] == vr
    if vr == "-":                     # insertion: VCF is <anchor> -> <anchor><inserted>
        return ref == alt[:1] and alt[1:] == va
    # delins, where VEP trims the anchor but leaves both sides non-empty (VCF ATT>AG at p is
    # VEP TT>G at p+1). These used to take the exact-match branch, which can never match the
    # VCF's anchored form, so they were silently dropped from read-level confirmation.
    return ref[1:] == vr and alt[1:] == va


def _vcf_alleles(child: str, variants: list[dict]) -> int:
    """Attach the true VCF position/REF/ALT to each candidate. Returns how many resolved.

    DuckDB stores VEP-style alleles — a pure insertion or deletion has '-' on one side, and
    the record is anchored one base to the RIGHT of the VCF record (the same offset
    `phasing._parse_phased` aliases in the other direction). A sites VCF for GATK needs real
    alleles, so each candidate is looked up in the child's own normalized VCF at its position
    and one base left. Candidates that cannot be resolved are dropped from the read-level
    confirmation rather than emitted as a malformed record."""
    import subprocess

    from .normalize import normalized_path

    regions = "\n".join(f"{v['chrom']}\t{max(1, v['pos'] - 1)}\t{v['pos']}"
                        for v in variants)
    proc = subprocess.run(
        ["bcftools", "query", "-R", "-", "-f", "%CHROM\t%POS\t%REF\t%ALT\n",
         str(normalized_path(child))],
        input=regions, capture_output=True, text=True)
    # Fails safe (candidates simply go unresolved and are dropped), but silently: a missing
    # or unindexed VCF looks exactly like "no candidate matched". Say which it was.
    if proc.returncode != 0:
        log.warning("bcftools query failed over %s (exit %d) — every candidate will look "
                    "unresolvable and be dropped from read-level confirmation: %s",
                    normalized_path(child).name, proc.returncode,
                    " | ".join(proc.stderr.strip().splitlines()[-2:]))
    found: dict[tuple[str, int], list[tuple[str, str]]] = {}
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) >= 4:
            found.setdefault((f[0], int(f[1])), []).append((f[2], f[3]))

    n = 0
    for v in variants:
        # VEP strips the anchor base a VCF indel record carries, so the two spellings never
        # compare equal directly and the record sits one base to the LEFT of VEP's position.
        anchored = "-" in (v["ref"], v["alt"]) or len(v["ref"]) != len(v["alt"])
        hit = None
        if not anchored:
            for ref, alt in found.get((v["chrom"], v["pos"]), []):
                if ref == v["ref"] and alt == v["alt"]:
                    hit = (v["pos"], ref, alt)
                    break
        else:
            for pos in (v["pos"] - 1, v["pos"]):
                for ref, alt in found.get((v["chrom"], pos), []):
                    if _same_event(v, ref, alt):
                        hit = (pos, ref, alt)
                        break
                if hit:
                    break
        if hit:
            v["vcf_pos"], v["vcf_ref"], v["vcf_alt"] = hit
            n += 1
        else:
            v["vcf_pos"] = None
    return n


# --- read-level confirmation of the absences a verdict rests on --------------
def _confirm_absences(child: str, parents: list[str], variants: list[dict],
                      carried: dict[str, set], force: bool = False) -> dict:
    """For each (parent, variant) the VCFs call *absent*, ask the parent's reads.

    Only the absences matter: a parent recorded as carrying a variant needs no confirmation,
    but "does not carry" from a variant-only VCF is precisely the claim that is unsafe, and
    every cis/trans verdict is built out of those claims."""
    wd = _work_dir(child)
    wd.mkdir(parents=True, exist_ok=True)
    need = [v for v in variants
            if any(v["key"] not in carried[p] for p in parents) and v.get("vcf_pos")]
    unresolved = [v for v in variants
                  if any(v["key"] not in carried[p] for p in parents) and not v.get("vcf_pos")]
    if unresolved:
        log.warning("segregation: %d candidate(s) could not be matched back to a record in "
                    "%s's normalized VCF — their absences stay unconfirmed.",
                    len(unresolved), child)
    if not need:
        return {}
    missing_cram = [p for p in parents if not aligned_cram(p).exists()]
    if missing_cram:
        log.warning("segregation: no CRAM for %s — absences stay unconfirmed and any verdict "
                    "resting on one is reported as unresolved.", ", ".join(missing_cram))
        return {}

    sites = [{"chrom": v["chrom"], "pos": v["vcf_pos"],
              "ref": v["vcf_ref"], "alt": v["vcf_alt"]} for v in need]
    vcf_sites, intervals, key = _write_sites(wd, sites)
    out: dict = {}
    for parent in parents:
        fc = wd / f"{parent}.parental.{key}.vcf.gz"
        if not (fc.exists() and fc.with_suffix(".gz.tbi").exists()) or force:
            call_at_sites(parent, aligned_cram(parent), vcf_sites, intervals, fc, key,
                          tag="segregation")
        calls = _read_parent(fc)
        regions = wd / f"sites_{key}.regions.txt"
        regions.write_text("".join(f"{v['chrom']}\t{v['vcf_pos']}\n" for v in need))
        pile = _pileup_parent(parent, regions, key, wd, force)
        for v in need:
            gt, dp, alt_ad = calls.get((v["chrom"], v["vcf_pos"]), ("", 0, 0))
            status, detail = parent_verdict(gt, dp, alt_ad)
            # denovo's vocabulary is about the *child's* variant being new; here the same
            # evidence answers "does this parent carry it", so translate once, at the source.
            status = {INHERITED: "carries", CONFIRMED: "absent"}.get(status, "unknown")
            nonref = pile.get((v["chrom"], v["vcf_pos"]))
            if nonref is not None and nonref >= config.DN_PILEUP_MIN_ALT:
                status, detail = "carries", f"{nonref} non-ref reads (pileup)"
            out[(parent, v["key"])] = (status, detail)
    return out


def transmission(dad: str, mum: str) -> str:
    """Which parent transmitted the child's alternate allele, from each parent's status.

    `dad`/`mum` are 'carries', 'absent' (evidenced hom-ref) or 'unknown' (reads could not
    resolve, or were never consulted). Anything involving 'unknown' stays UNRESOLVED — the
    whole point of the stage is that a guess here silently flips a cis/trans call."""
    if "unknown" in (dad, mum):
        return UNRESOLVED
    if dad == "carries" and mum == "carries":
        return BOTH
    if dad == "carries":
        return PATERNAL
    if mum == "carries":
        return MATERNAL
    return DE_NOVO


def resolve(variants: list[dict]) -> tuple[str, str]:
    """(verdict, why) for one gene's damaging het variants, from their transmissions."""
    tr = [v["transmission"] for v in variants]
    if any(t == UNRESOLVED for t in tr):
        return UNRESOLVED, "at least one variant's parental origin could not be established"
    if any(t == DE_NOVO for t in tr):
        return UNRESOLVED, ("a variant appears in neither parent — a de-novo event sits on "
                            "one copy, but which one needs the `denovo` stage's confirmation")
    if any(t == BOTH for t in tr):
        return UNRESOLVED, ("both parents carry at least one of these, so the child's two "
                            "copies cannot be told apart by transmission alone")
    if PATERNAL in tr and MATERNAL in tr:
        return IN_TRANS, "one variant from each parent, so they sit on opposite copies"
    parent = "father" if tr[0] == PATERNAL else "mother"
    return IN_CIS, f"every variant came from the {parent}, so they share one copy"


# --- report ------------------------------------------------------------------
def _report(child: str, parents: list[str], genes: dict, snap: str) -> Path:
    trans = {g: r for g, r in genes.items() if r["verdict"] == IN_TRANS}
    lines = [
        f"# Trio segregation & compound heterozygosity — {child}",
        "",
        f"- **Trio:** child `{child}`  ·  father `{parents[0]}`  ·  mother `{parents[1]}`",
        f"- **Snapshot:** `{snap}`",
        f"- **Panel genes with ≥2 damaging hets:** {len(genes)}",
        f"- **Resolved in-trans (both copies hit — review):** {len(trans)}",
        "",
        "## Verdict by gene",
        "",
        "| Gene | Damaging hets | Verdict | Basis |",
        "|------|--------------:|---------|-------|",
        *([f"| {g} | {r['n']} | {r['verdict']} | {r['why']} |"
           for g, r in sorted(genes.items())] or ["| _none_ | | | |"]),
        "",
        "## Transmission detail",
        "",
    ]
    for gene, r in sorted(genes.items()):
        lines += [f"### {gene} — {r['verdict']}", "",
                  "| Variant | Consequence | ClinVar | Father | Mother | Origin |",
                  "|---------|-------------|---------|--------|--------|--------|"]
        for v in r["variants"]:
            lines.append(
                f"| {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']} | {v['csq']} | "
                f"{v['sig'] or '—'} | {v['dad_detail']} | {v['mum_detail']} | "
                f"{v['transmission']} |")
        lines.append("")
    lines += [
        "## Method & caveats",
        "",
        "- **Why this exists alongside `phase`.** Read-backed phasing resolves cis/trans only "
        "for variants a read or read-pair spans (~1 kb); transmission resolves them at any "
        "distance, and needs no reads from the child. Where both stages have an opinion they "
        "should agree — a disagreement is worth investigating, not averaging.",
        "- **Absence is confirmed against reads, not assumed.** Both parents' VCFs are "
        "variant-only, so a missing record is either hom-ref or never called. Every absence a "
        "verdict depends on is re-checked against the parent's CRAM (GATK force-call plus a "
        "caller-independent raw pileup, as in `denovo`). Without a parental CRAM, or where the "
        "reads cannot resolve the site, the verdict is **unresolved** rather than guessed.",
        "- **`in-trans` is the reviewable outcome:** both copies of a recessive-disease gene "
        "carry a damaging variant. **`in-cis`** means both sit on one copy, leaving the other "
        "intact — usually reassuring.",
        "- **Damaging ≠ pathogenic.** \"Damaging\" here is LOF or deleterious-missense "
        "(SIFT/PolyPhen), the same definition the panel and `phase` use. Confirm the variants "
        "themselves before acting on a configuration.",
        "- **Both parents carrying a variant makes it unresolvable by transmission** — the "
        "child's two copies are then indistinguishable by origin. Reads (`phase`) are the "
        "only recourse for those.",
        "- Research-grade, not a clinical diagnosis.",
        "",
    ]
    out = report_path(child)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_report(out, lines)
    return out


# --- entry point --------------------------------------------------------------
def segregation_report(child: str, force: bool = False) -> Path:
    """Resolve the child's panel-gene compound-het candidates by parental transmission."""
    fam = pedigree.load(sequenced=pedigree.sequenced_ids())
    trio = fam.trio(child)
    if not trio:
        raise SystemExit(
            f"No trio for {child}: both parents must be set in {config.FAMILY_FILE} and "
            f"sequenced. Edit the family graph, then re-run.")
    _c, father, mother = trio
    parents = [father, mother]

    out = report_path(child)
    if out.exists() and not force:
        log.info("%s already has a segregation report (%s) — skipping (use --force).",
                 child, out)
        return out

    snap = latest_snapshot()
    if not snap:
        raise SystemExit("No snapshot yet. Run `snapshot` first.")
    panel_genes = {g for d in load_panel() for g in d.genes}
    con = connect()
    variants = _child_candidates(con, snap, child, panel_genes)
    if not variants:
        con.close()
        log.info("segregation: %s carries no panel gene with ≥2 damaging hets — "
                 "nothing to resolve.", child)
        return _report(child, parents, {}, snap)

    for v in variants:
        v["key"] = (v["chrom"], v["pos"], v["ref"], v["alt"])
    callsets = {p: parent_callsets(p) for p in parents}
    recorded = _parent_calls(con, snap, [c for cs in callsets.values() for c in cs], variants)
    con.close()
    carried = {p: {v["key"] for v in variants
                   if any((cs, *v["key"]) in recorded for cs in callsets[p])}
               for p in parents}
    resolved = _vcf_alleles(child, variants)
    log.info("segregation: %d damaging het(s) across %d gene(s) (%d matched back to VCF "
             "records); confirming absences off the parents' reads.",
             len(variants), len({v["gene"] for v in variants}), resolved)
    confirmed = _confirm_absences(child, parents, variants, carried, force)

    for v in variants:
        statuses = {}
        for parent in parents:
            if v["key"] in carried[parent]:
                statuses[parent] = ("carries", "in their VCF")
                continue
            statuses[parent] = confirmed.get((parent, v["key"]),
                                             ("unknown", "not confirmed against reads"))
        v["dad_detail"] = f"{statuses[father][0]} ({statuses[father][1]})"
        v["mum_detail"] = f"{statuses[mother][0]} ({statuses[mother][1]})"
        v["transmission"] = transmission(statuses[father][0], statuses[mother][0])

    genes: dict[str, dict] = {}
    for v in variants:
        genes.setdefault(v["gene"], {"variants": []})["variants"].append(v)
    for gene, r in genes.items():
        r["n"] = len(r["variants"])
        r["verdict"], r["why"] = resolve(r["variants"])

    report = _report(child, parents, genes, snap)
    n_trans = sum(1 for r in genes.values() if r["verdict"] == IN_TRANS)
    log.info("Segregation report: %s (%d gene(s), %d in-trans)", report, len(genes), n_trans)
    return report
