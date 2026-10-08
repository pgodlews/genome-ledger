"""Stage: single-SNP traits / wellness lookups.

Reads a curated rsID → genotype → interpretation rules file (`traits/traits.tsv`,
GRCh38 forward-strand coords from Ensembl) and reports the sample's result for each.
The sample's genotype at each position is read straight from the called VCF; positions
absent from the (variant-only) VCF are taken as homozygous reference — fine for the
high-coverage WGS here, the same assumption PharmCAT's `--absent-to-ref` makes.

Deterministic, well-established single-SNP traits only (eye colour, earwax, lactase,
caffeine, Factor V Leiden, …); polygenic traits belong in the PRS engine (P3).
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

from . import reportpaths
from . import config, forcecall
from .util import PASS_ONLY, log, read_tsv, require_tools, write_report


def _called_vcf(sample: str) -> Path:
    """This sample's VCF — the pipeline's own call, or the normalized ingested one."""
    from .call import sample_vcf

    return sample_vcf(sample)


def _alleles_at(vcf: Path, chrom: str, pos: str, ref: str) -> tuple[str, str] | None:
    """The sample's two alleles at chrom:pos, or None on an explicit no-call. GATK emits
    only variant sites, so an absent position means homozygous reference — but a record
    that IS present with a ./. genotype is uncertainty, not evidence of reference."""
    # PASS-filtered, like every other engine that reads genotypes out of a called VCF. Traits
    # and HLA queried unfiltered, so a call the pipeline's own hard filters rejected (MQ40,
    # QD2, SOR3 — ~275k records, ~5.6% of a typical callset) was read here as a real genotype
    # and could decide a trait or an HLA tag. See util.PASS_ONLY.
    proc = subprocess.run(
        ["bcftools", "query", "-r", f"{chrom}:{pos}", *PASS_ONLY,
         "-f", "%POS\t%REF\t%ALT[\t%GT]\n", str(vcf)],
        capture_output=True, text=True)
    # A failed query must not fall through to the `return ref, ref` below: an unreadable VCF,
    # a missing index or an unknown contig would otherwise be reported as a confident
    # homozygous-reference genotype. Only a query that actually SUCCEEDED and found no record
    # licenses the absent⇒hom-ref convention.
    if proc.returncode != 0:
        log.warning("bcftools query failed at %s:%s in %s (exit %d) — reporting the trait as "
                    "indeterminate rather than hom-ref: %s", chrom, pos, vcf.name,
                    proc.returncode, " | ".join(proc.stderr.strip().splitlines()[-2:]))
        return None
    out = proc.stdout.strip()
    if not out:
        alt_chrom = f"chr{chrom}" if not chrom.startswith("chr") else chrom.removeprefix("chr")
        proc2 = subprocess.run(
            ["bcftools", "query", "-r", f"{alt_chrom}:{pos}", *PASS_ONLY,
             "-f", "%POS\t%REF\t%ALT[\t%GT]\n", str(vcf)],
            capture_output=True, text=True)
        if proc2.returncode == 0 and proc2.stdout.strip():
            out = proc2.stdout.strip()
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 4:   # malformed/truncated bcftools output line — skip, don't crash
            continue
        p, vref, valt, gt = f[:4]
        if p != pos:
            continue
        idx = re.split(r"[/|]", gt)
        if any(g in (".", "") for g in idx):
            return None  # explicit (partial) no-call — not evidence of hom-ref
        coding = [vref, *valt.split(",")]
        try:
            alleles = [coding[int(g)] for g in idx]
        except (ValueError, IndexError):
            return None  # malformed record (e.g. '*' allele w/o matching ALT) — don't guess
        if len(alleles) >= 2:
            return alleles[0], alleles[1]
        if alleles:  # haploid (not expected for autosomal traits)
            return alleles[0], alleles[0]
    return ref, ref  # absent → hom-ref


def traits_report(sample: str, force: bool = False) -> Path:
    require_tools("bcftools")
    rules = read_tsv(config.TRAITS_FILE)
    if not rules:
        raise SystemExit(f"Traits rules file empty/missing: {config.TRAITS_FILE}")
    vcf = _called_vcf(sample)
    if not vcf.exists():
        raise SystemExit(f"{vcf} missing. Run `call {sample}` first.")

    out_dir = config.TRAITS_DIR / sample
    summary = reportpaths.genome_report("traits", sample)
    if summary.exists() and not force:
        log.info("%s traits already reported (%s) — skipping (use --force).", sample, summary)
        return summary
    out_dir.mkdir(parents=True, exist_ok=True)

    # Prefer force-called genotypes when present: an evidenced 0/0 (not an assumed one),
    # and an explicit no-call where the reads couldn't resolve the site. Without a
    # force-call VCF, fall back to the called VCF + absent→hom-ref convention.
    forced = forcecall.read_genotypes(sample)
    if forced:
        log.info("traits: using force-called genotypes (%d sites)", len(forced))

    rows = []
    for r in rules:
        g = forced.get((r["chrom"], r["pos"]))
        if g and g["nocall"]:
            rows.append((r["trait"], "./.", r["rsid"],
                         "Indeterminate — no-call (low coverage at this position)"))
            log.info("  %-38s ./.  → no-call (low coverage)", r["trait"])
            continue
        pair = (g["alleles"][0], g["alleles"][1]) if g else \
            _alleles_at(vcf, r["chrom"], r["pos"], r["ref"])
        if pair is None:
            rows.append((r["trait"], "./.", r["rsid"],
                         "Indeterminate — no-call (low coverage at this position)"))
            log.info("  %-38s ./.  → no-call (low coverage)", r["trait"])
            continue
        a1, a2 = pair
        n_effect = [a1, a2].count(r["effect"])
        interp = r[f"interp{n_effect}"]
        rows.append((r["trait"], f"{a1}/{a2}", r["rsid"], interp))
        log.info("  %-38s %s/%s → %s", r["trait"], a1, a2, interp)

    lines = [
        f"# Traits & wellness — {sample}",
        "",
        f"- **Single-SNP lookups:** {len(rows)}  ·  rules: `{config.TRAITS_FILE.name}`",
        "",
        "| Trait | Genotype | rsID | Result |",
        "|-------|:--------:|------|--------|",
        *[f"| {t} | {g} | {rs} | {res} |" for t, g, rs, res in rows],
        "",
        "## Caveats",
        ("- Genotypes are read from the **force-called** VCF (evidenced hom-ref; no-calls "
         "shown as indeterminate)." if forced else
         "- Positions absent from the variant-only VCF are taken as homozygous reference "
         "(high-coverage WGS assumption). Run `force-call` to make this rigorous."),
        "- Single-SNP, well-established associations only — not polygenic predictions (see PRS).",
        "- Research-grade. Actionable items (e.g. Factor V Leiden) warrant clinical confirmation.",
        "",
    ]
    write_report(summary, lines)
    log.info("Traits report: %s (%d traits)", summary, len(rows))
    return summary
