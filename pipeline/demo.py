"""Demo-only reports for the stages a synthetic family cannot actually run.

`run.py demo` promises zero downloads and under a minute. Five report types cannot be
produced honestly inside that promise, because each needs something the demo family does not
have:

  - repeat expansions   — ExpansionHunter sizes repeats from the reads (a CRAM)
  - pharmacogenomics    — PharmCAT (a Java pipeline) over a genome-wide VCF
  - haplogroups         — Haplogrep / Yleaf over real MT and Y calls
  - polygenic scores    — PGS Catalog scoring files and a 1000 Genomes reference panel
  - (HLA tag SNPs need none of that, and are run for real by the demo itself)

So instead of the tools, the demo ships their *outputs*: tests/fixtures/demo/toy_results.json
holds hand-written sample values, which this module writes out in each tool's own file
format and hands to the pipeline's real parsers, interpretation tables, consent gate and
report writers. What a demo user sees is therefore the genuine report, with genuine
thresholds and wording, over values that are illustrative — and every such page says so
under its title (config.DEMO_NOTICE).
"""

from __future__ import annotations

import json
from contextlib import contextmanager
from pathlib import Path
from statistics import NormalDist

from . import ancestry, config, pgx, prs, repeats, reportpaths
from .util import log

NOTICE = ("**Demo only — synthetic family.** These people do not exist: every genotype was "
          "written by hand to exercise the pipeline. Nothing on this page is a finding about "
          "a real person.")

DELTA_NOTICE = ("**Demo only — synthetic family, synthetic history.** The comparison below is "
                "computed for real by the pipeline's diff and triage code, but both snapshots "
                "are invented: the earlier one is the same family with three variant "
                "classifications rolled back, to imitate ClinVar moving between releases. "
                "Nothing here is a finding about a real person or a real ClinVar change.")

EARLIER_NOTICE = ("**Demo only — synthetic family, synthetic history.** This is the demo's "
                  "invented *earlier* snapshot: the same made-up family, with three variant "
                  "classifications rolled back to imitate an older ClinVar release. Nothing "
                  "here is a finding about a real person.")


def _toy_notice(tool: str) -> str:
    return ("**Demo only — illustrative values.** The demo family is synthetic and has no "
            f"sequencing reads, and the demo downloads nothing, so {tool} was not run. The "
            "values below are bundled sample results passed through the pipeline's real "
            "interpretation and report code, to show what this report looks like. They are "
            "not computed from the demo genomes and describe no real person.")


@contextmanager
def notice(text: str):
    """Set config.DEMO_NOTICE for the duration of a block, restoring what was there."""
    saved = config.DEMO_NOTICE
    config.DEMO_NOTICE = text
    try:
        yield
    finally:
        config.DEMO_NOTICE = saved


def _repeats(sample: str, toy: dict) -> None:
    out_dir = config.REPEATS_DIR / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    # The same shape ExpansionHunter writes: one record per locus, VARID/RU in INFO and the
    # per-allele repeat counts in the REPCN format field.
    body = ["##fileformat=VCFv4.1",
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"]
    for i, (locus, unit) in enumerate(toy["loci"].items(), start=1):
        repcn = "/".join(str(n) for n in toy["sizes"][sample][locus])
        body.append(f"chr1\t{i}\t.\tN\t<STR>\t.\tPASS\tVARID={locus};RU={unit}\t"
                    f"GT:REPCN\t0/1:{repcn}")
    vcf = out_dir / f"{sample}.vcf"
    vcf.write_text("\n".join(body) + "\n")
    repeats.report_from_genotypes(sample, vcf)


def _pgx(sample: str, toy: dict) -> None:
    out_dir = config.PGX_DIR / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    genes = {}
    for gene, label in toy["diplotypes"][sample].items():
        genes[gene] = {"recommendationDiplotypes": [
            {"label": label, "phenotypes": [toy["phenotypes"][gene][label]]}]}
    genes["CYP2D6"] = {"recommendationDiplotypes": []}   # never callable from a plain VCF
    report_json = out_dir / f"{sample}.chr.report.json"
    report_json.write_text(json.dumps({"pharmcatVersion": toy["pharmcatVersion"],
                                       "dataVersion": toy["dataVersion"], "genes": genes}))
    pgx._write_summary(sample, report_json, reportpaths.genome_report("pgx", sample),
                       "there are no reads in the demo")


def _ancestry(sample: str, toy: dict) -> None:
    out_dir = config.ANCESTRY_DIR / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    hg, quality = toy[sample]["mt"]
    mt_file = out_dir / f"{sample}.mt.hg"
    mt_file.write_text('"SampleID"\t"Haplogroup"\t"Rank"\t"Quality"\t"Range"\n'
                       f'"{sample}"\t"{hg}"\t"1"\t"{quality}"\t"1-16569"\n')
    y = None
    if "y" in toy[sample]:
        yhg, qc = toy[sample]["y"]
        y_file = out_dir / "hg_prediction.hg"
        y_file.write_text("Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\tQC-score\n"
                          f"{sample}\t{yhg}\t{yhg.split('-')[-1]}\tNA\t1200\t{qc}\n")
        y = ancestry._parse_yleaf(sample, y_file)
    ancestry.write_ancestry(sample, ancestry._parse_haplogrep(sample, mt_file), y)


def _prs(sample: str, toy: dict) -> None:
    # The real panel and the real consent/tier filter decide which traits appear; only the
    # scores themselves are bundled.
    panel, notes = prs._panel_filters(
        [e for e in prs._panel_scores() if e["pgs_id"] in toy["traits"]], sample)
    rows = []
    for entry in panel:
        pid = entry["pgs_id"]
        pct = float(toy["eur_percentile"][sample][pid])
        dist = toy["traits"][pid]
        total = int(entry["n_variants"])
        rows.append({**entry, "pid": pid,
                     "method": "plink2" if total > config.PRS_PLINK_MIN else "exact",
                     "raw": dist["mean"] + NormalDist().inv_cdf(pct / 100) * dist["sd"],
                     "used": total, "total": total, "pe": pct, "pa": None})
    prs.write_prs(sample, rows, notes, [], forced=False)


def toy_reports(samples: list[str], fixtures: Path) -> None:
    """Write the repeats, PGx, haplogroup and PRS reports for the demo family."""
    toy = json.loads((fixtures / "toy_results.json").read_text())
    for sample in samples:
        with notice(_toy_notice("ExpansionHunter")):
            _repeats(sample, toy["repeats"])
        with notice(_toy_notice("PharmCAT")):
            _pgx(sample, toy["pgx"])
        with notice(_toy_notice("Haplogrep / Yleaf")):
            _ancestry(sample, toy["haplogroups"])
        with notice(_toy_notice("polygenic scoring (PGS Catalog files against a 1000 "
                                "Genomes reference panel)")):
            _prs(sample, toy["prs"])
    log.info("Demo: wrote illustrative repeats / PGx / haplogroup / PRS reports for %d "
             "people.", len(samples))
