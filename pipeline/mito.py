"""Stage: mitochondrial heteroplasmy via GATK Mutect2 mito-mode.

Haploid HaplotypeCaller (the `call` stage) gives a present/absent germline mtDNA call, good
enough for the haplogroup. It cannot say *what fraction* of a person's mtDNA molecules carry
a variant — the **heteroplasmy** level. That fraction is the clinically meaningful number for
mitochondrial disease, which is dosage-dependent: a pathogenic allele only causes disease
above a tissue-specific threshold, and the same mutation can be silent at 20% and severe at 80%.

GATK **Mutect2 --mitochondria-mode** treats the mtDNA like a somatic sample, estimating the
alt-allele fraction (AF) at each site straight off the CRAM; **FilterMutectCalls
--mitochondria-mode** flags low-confidence calls. We split multiallelics, read AF as the
heteroplasmy fraction, classify each PASS variant homoplasmic vs heteroplasmic, and flag any
that match a small curated table of well-established pathogenic mtDNA point mutations.

Runs off the CRAM only (no ClinVar), so the report is filed per genome, outside the dated snapshots.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from . import reportpaths
from . import config
from .call import _java_env, aligned_cram
from .util import log, read_tsv, require_tools, run, write_report


def _curated() -> dict[tuple[str, str, str], dict]:
    """(pos, ref, alt) → row from the curated pathogenic-mtDNA table."""
    return {(r["pos"], r["ref"], r["alt"]): r for r in read_tsv(config.MITO_VARIANTS_FILE)}


def _classify(af: float) -> str:
    return "homoplasmic" if af >= config.MITO_HOMOPLASMIC_AF else "heteroplasmic"


def _parse_calls(split_vcf: Path) -> list[dict]:
    """Read (pos, ref, alt, filter, af, dp) per ALT from the split, filtered Mutect2 VCF.
    The VCF is biallelic-per-line (post `bcftools norm -m-`), so AF/DP are single values."""
    out = subprocess.run(
        ["bcftools", "query", "-f",
         "%POS\t%REF\t%ALT\t%FILTER[\t%AF\t%DP]\n", str(split_vcf)],
        capture_output=True, text=True, check=True).stdout
    calls = []
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 6:
            continue
        pos, ref, alt, filt, af, dp = f[:6]
        try:
            af_f = float(af)
        except ValueError:
            continue
        calls.append({"pos": pos, "ref": ref, "alt": alt, "filter": filt,
                      "af": af_f, "dp": int(dp) if dp.isdigit() else None})
    return calls


def _run_mutect2(sample: str, cram: Path, out_dir: Path, env: dict) -> Path:
    """Mutect2 mito-mode → FilterMutectCalls → split multiallelics. Returns the split VCF."""
    raw = out_dir / f"{sample}.mt.raw.vcf.gz"
    filtered = out_dir / f"{sample}.mt.filtered.vcf.gz"
    split = out_dir / f"{sample}.mt.split.vcf.gz"
    logf = config.LOGS_DIR / f"mito_{sample}.log"

    log.info("Mutect2 --mitochondria-mode on %s (%s) …", sample, config.MT_CONTIG)
    # --mitochondria-mode tunes Mutect2 for the single high-depth haploid contig (low LOD
    # thresholds → sensitive to low heteroplasmy). We keep its default per-start downsampling,
    # matching the Broad mtDNA best-practice pipeline; raising --max-reads-per-alignment-start
    # would buy a little low-AF sensitivity at a large runtime cost (the more so here, where
    # this Mac lacks the AVX-native PairHMM and falls back to the slow Java implementation —
    # like every GATK step, mito runs faster on the Linux box).
    run(["java", f"-Xmx{config.HC_JAVA_MEM}", "-jar", str(config.GATK_JAR), "Mutect2",
         "-R", str(config.REF_FASTA), "-I", str(cram), "-L", config.MT_CONTIG,
         "--mitochondria-mode", "-O", str(raw)],
        env=env, stdout=open(logf, "wb"), stderr=subprocess.STDOUT)

    log.info("FilterMutectCalls --mitochondria-mode …")
    run(["java", f"-Xmx{config.HC_JAVA_MEM}", "-jar", str(config.GATK_JAR),
         "FilterMutectCalls", "-R", str(config.REF_FASTA), "--mitochondria-mode",
         "-V", str(raw), "-O", str(filtered)],
        env=env, stdout=open(logf, "ab"), stderr=subprocess.STDOUT)

    # Split multiallelics so AF (Number=A) becomes one heteroplasmy fraction per record.
    run(["bcftools", "norm", "-m-", "-f", str(config.REF_FASTA), "-Oz", "-o", str(split),
         str(filtered)], stdout=open(logf, "ab"), stderr=subprocess.STDOUT)
    return split


def mito_report(sample: str, force: bool = False) -> Path:
    require_tools("java", "bcftools")
    if not config.GATK_JAR.exists():
        raise SystemExit(f"GATK not installed at {config.GATK_JAR}. Run `setup`.")
    cram = aligned_cram(sample)
    if not cram.exists():
        raise SystemExit(f"{cram} missing. Run `align {sample}` first.")

    out_dir = config.MITO_DIR / sample
    summary = reportpaths.genome_report("mito", sample)
    if summary.exists() and not force:
        log.info("%s mito heteroplasmy already reported (%s) — skipping (use --force).",
                 sample, summary)
        return summary
    out_dir.mkdir(parents=True, exist_ok=True)

    split = _run_mutect2(sample, cram, out_dir, _java_env())
    calls = _parse_calls(split)
    curated = _curated()
    for c in calls:
        c["type"] = _classify(c["af"])
        c["known"] = curated.get((c["pos"], c["ref"], c["alt"]))

    passing = [c for c in calls if c["filter"] in ("PASS", ".")]
    passing.sort(key=lambda c: c["af"], reverse=True)
    het = [c for c in passing
           if c["type"] == "heteroplasmic" and c["af"] >= config.MITO_MIN_HET_AF]
    homo = [c for c in passing if c["type"] == "homoplasmic"]
    low = [c for c in passing if c["af"] < config.MITO_MIN_HET_AF]
    pathogenic = [c for c in passing if c["known"]]
    filtered_out = len(calls) - len(passing)

    def _row(c: dict) -> str:
        m = c["known"]
        name = m["mutation"] if m else f"m.{c['pos']}{c['ref']}>{c['alt']}"
        note = f"⚠️ **{m['locus']} — {m['disease']}**" if m else ""
        dp = f"{c['dp']}x" if c["dp"] is not None else "?"
        return (f"| {name} | {c['af'] * 100:.1f}% | {c['type']} | {dp} | {note} |")

    lines = [
        f"# Mitochondrial heteroplasmy — {sample}",
        "",
        f"- **Tool:** GATK Mutect2 {config.GATK_VERSION} `--mitochondria-mode` "
        f"+ FilterMutectCalls (contig `{config.MT_CONTIG}`, rCRS)",
        f"- **PASS variants:** {len(passing)} — {len(homo)} homoplasmic, "
        f"{len(het)} heteroplasmic (≥{int(config.MITO_MIN_HET_AF * 100)}%)",
        f"- **Known pathogenic mtDNA variants detected:** {len(pathogenic)}",
        "",
    ]

    if pathogenic:
        lines += [
            "## ⚠️ Known pathogenic mtDNA variants",
            "_Curated well-established disease mutations. **Heteroplasmy level matters** — "
            "pathogenicity is threshold-dependent, so confirm clinically and interpret the "
            "fraction in the affected tissue (blood may under-represent other tissues)._",
            "",
            "| Variant | Heteroplasmy | Type | Depth | Disease |",
            "|---------|-------------:|------|------:|---------|",
            *[_row(c) for c in pathogenic],
            "",
        ]

    lines += [
        "## Heteroplasmic variants",
        ("_Variants present in only a fraction of mtDNA molecules — the novel signal over a "
         "haploid germline call._" if het else "_None above the reporting threshold._"),
        "",
    ]
    if het:
        lines += [
            "| Variant | Heteroplasmy | Type | Depth | Note |",
            "|---------|-------------:|------|------:|------|",
            *[_row(c) for c in het],
            "",
        ]

    lines += [
        f"## Homoplasmic variants ({len(homo)})",
        "_Present in ~all mtDNA copies — mostly haplogroup-defining polymorphisms (differences "
        "from the rCRS reference), not disease. Listed for completeness._",
        "",
        "| Variant | Heteroplasmy | Type | Depth | Note |",
        "|---------|-------------:|------|------:|------|",
        *[_row(c) for c in homo],
        "",
        "## Caveats",
        "- **Heteroplasmy is tissue- and time-variable.** This is the fraction in the "
        "sequenced sample (blood); other tissues (muscle, brain) can differ substantially, "
        "and blood levels of some pathogenic variants decline with age.",
        f"- **Low-level calls (<{int(config.MITO_MIN_HET_AF * 100)}%)** are excluded from the "
        f"headline counts ({len(low)} such PASS call(s)); at high mtDNA depth very low "
        "fractions are easily artefacts (NUMT misalignment, sequencing error).",
        f"- **{filtered_out} call(s)** did not PASS FilterMutectCalls and are omitted.",
        "- **NUMTs** (nuclear-mitochondrial segments) can masquerade as low-level heteroplasmy; "
        "Mutect2 mito-mode mitigates but does not eliminate this.",
        "- **Control-region homopolymers** (the poly-C / AC-repeat tracts around m.302–315 and "
        "m.16182–16194) routinely show apparent length heteroplasmy (insertions/deletions of a "
        "few bases) that is an alignment artefact, not biology — discount indel \"heteroplasmy\" "
        "in these tracts.",
        "- Research-grade, not a clinical diagnosis.",
        "",
    ]
    write_report(summary, lines)
    log.info("Mito report: %s (%d PASS: %d homoplasmic, %d heteroplasmic, %d known-pathogenic)",
             summary, len(passing), len(homo), len(het), len(pathogenic))
    return summary
