"""Stage: HLA risk via tag SNPs.

Reproduces the HLA-driven report items — **HLA-B27** (ankylosing spondylitis) and the
**celiac** DQ2.5 / DQ8 susceptibility haplotypes — from validated tag SNPs read out of
the called VCF, the same pragmatic approach consumer reports use. This is NOT full
classical HLA typing (no 4-digit alleles, no loci beyond B/DQA1/DQB1) — see
`pipeline/hla_typing.py` (`hla-type` stage, arcasHLA) for that.
"""

from __future__ import annotations

import shutil

from . import reportpaths
from . import config, forcecall
from .traits import _alleles_at, _called_vcf
from .util import log, read_tsv, require_tools, write_report

# The tag rules, under the short names the typed report compares them by.
_TAG_TRAITS = {"B27": "HLA-B27 (ankylosing spondylitis)",
               "DQ2.5": "Celiac DQ2.5 haplotype",
               "DQ8": "Celiac DQ8 haplotype"}


def _read_tags(rules: list[dict], vcf, forced: dict) -> tuple[list, dict, set]:
    """(report rows, trait -> carrier, traits with a no-call) for the tag SNPs in `rules`."""
    rows, present, nocalls = [], {}, set()
    for r in rules:
        g = forced.get((r["chrom"], r["pos"]))
        pair = None if (g and g["nocall"]) else (
            (g["alleles"][0], g["alleles"][1]) if g
            else _alleles_at(vcf, r["chrom"], r["pos"], r["ref"]))
        if pair is None:
            rows.append((r["trait"], "./.", r["rsid"],
                         "Indeterminate — no-call (low coverage at this position)"))
            present[r["trait"]] = False
            nocalls.add(r["trait"])
            log.info("  %-38s ./.  → no-call (low coverage)", r["trait"])
            continue
        a1, a2 = pair
        n = [a1, a2].count(r["effect"])
        rows.append((r["trait"], f"{a1}/{a2}", r["rsid"], r[f"interp{n}"]))
        present[r["trait"]] = n > 0
        log.info("  %-38s %s/%s → %s", r["trait"], a1, a2, r[f"interp{n}"])
    return rows, present, nocalls


def tag_status(sample: str) -> dict[str, bool | None]:
    """{B27, DQ2.5, DQ8} -> carrier, non-carrier, or None for a no-call. Empty when the tag
    SNPs cannot be read here, so a caller can tell "no cross-check" from a negative."""
    rules = read_tsv(config.HLA_FILE)
    vcf = _called_vcf(sample)
    if not rules or not vcf.exists() or shutil.which("bcftools") is None:
        return {}
    _, present, nocalls = _read_tags(rules, vcf, forcecall.read_genotypes(sample))
    return {key: None if trait in nocalls else present[trait]
            for key, trait in _TAG_TRAITS.items() if trait in present}


def hla_report(sample: str, force: bool = False):
    require_tools("bcftools")
    rules = read_tsv(config.HLA_FILE)
    if not rules:
        raise SystemExit(f"HLA tag rules missing: {config.HLA_FILE}")
    vcf = _called_vcf(sample)
    if not vcf.exists():
        raise SystemExit(f"{vcf} missing. Run `call {sample}` first.")

    out_dir = config.HLA_DIR / sample
    summary = reportpaths.genome_report("hla", sample)
    if summary.exists() and not force:
        log.info("%s HLA already reported (%s) — skipping (use --force).", sample, summary)
        return summary
    out_dir.mkdir(parents=True, exist_ok=True)

    # Prefer force-called genotypes exactly as `traits` does. forcecall._collect_sites()
    # already puts HLA_FILE's tag SNPs into the force-call union, so without this the stage
    # computed evidenced genotypes and then threw them away: an unresolved tag SNP fell back
    # to absent→hom-ref, and a DQ2.5/DQ8 no-call was reported as "neither haplotype present"
    # — a false negative in precisely the case the indeterminate wording below exists for.
    forced = forcecall.read_genotypes(sample)
    if forced:
        log.info("hla: using force-called genotypes (%d sites)", len(forced))

    rows, present, nocalls = _read_tags(rules, vcf, forced)

    dq25 = present.get("Celiac DQ2.5 haplotype", False)
    dq8 = present.get("Celiac DQ8 haplotype", False)
    if not dq25 and not dq8:
        if nocalls & {"Celiac DQ2.5 haplotype", "Celiac DQ8 haplotype"}:
            celiac = ("Celiac-susceptibility haplotype status is **indeterminate** — a "
                      "no-call at a DQ2.5/DQ8 tag SNP means absence cannot be concluded. "
                      "Re-run after improving coverage (or use full HLA typing) before "
                      "reading anything into a negative result.")
        else:
            # A tag negative is narrower than "no risk haplotype". The DQ2.5 tag marks the cis
            # haplotype only, so DQ2.5 formed in trans and DQ2.2 both read as absent here.
            celiac = ("Celiac susceptibility is **not ruled out** — neither tagged haplotype "
                      "was detected. The DQ2.5 tag sees only the cis haplotype "
                      "(`DQA1*05:01`-`DQB1*02:01`); it cannot see DQ2.5 formed in trans "
                      "(`DQA1*05:05` with `DQB1*02:02`) or DQ2.2. Use full HLA typing "
                      "(`hla-type`) to settle it.")
    else:
        hap = " + ".join(h for h, p in [("DQ2.5", dq25), ("DQ8", dq8)] if p)
        celiac = (f"Carries a celiac-susceptibility haplotype (**{hap}**). Necessary but "
                  "not sufficient — most carriers never develop celiac; pursue serological/"
                  "clinical testing only if symptomatic.")

    lines = [
        f"# HLA risk (tag SNPs) — {sample}",
        "",
        "| Marker | Genotype | rsID | Result |",
        "|--------|:--------:|------|--------|",
        *[f"| {t} | {g} | {rs} | {res} |" for t, g, rs, res in rows],
        "",
        "## Celiac disease — HLA summary",
        celiac,
        "",
        "## Caveats",
        "- **Tag-SNP based, not full classical HLA typing.** Tag SNPs infer the common "
        "risk haplotypes well but don't give 4-digit alleles. For full HLA-A/B/C/DQ/DR "
        "typing use arcasHLA (from the CRAM) or HIBAG — a heavier follow-up.",
        ("- Genotypes are read from the **force-called** VCF (evidenced hom-ref; a site the "
         "reads could not resolve is reported as indeterminate, not as absence of risk)."
         if forced else
         "- Positions absent from the variant-only VCF are taken as homozygous reference. "
         "Run `force-call` to make this rigorous — until then a no-call at a tag SNP is "
         "indistinguishable from a genuine negative."),
        "- Research-grade; HLA-B27 / celiac status warrants clinical confirmation if relevant.",
        "",
    ]
    write_report(summary, lines)
    log.info("HLA report: %s", summary)
    return summary
