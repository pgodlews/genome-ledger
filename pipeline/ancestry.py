"""Stage: ancestry haplogroups.

Two lineage markers from a sample's existing calls:
  - **mtDNA haplogroup** (maternal line) via Haplogrep 3 over the MT calls.
  - **Y haplogroup** (paternal line, males only) via Yleaf (hg38-native) over the Y calls.

Both depend only on the genome (not ClinVar), so the report is filed per genome, not per snapshot.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from . import reportpaths
from . import config
from .call import _java_env
from .util import find_row, log, require_tools, run, write_report


def _called_vcf(sample: str) -> Path:
    """This sample's VCF — the pipeline's own call, or the normalized ingested one."""
    from .call import sample_vcf

    return sample_vcf(sample)


def _sex(sample: str) -> str:
    row = find_row(config.SAMPLES_TSV, sample_id=sample) or \
        find_row(config.INCOMING_MANIFEST, sample=sample) or {}
    return row.get("sex", "").lower()


def _subset(vcf: Path, contig: str, dest: Path) -> Path:
    run(["bcftools", "view", str(vcf), contig, "-Oz", "-o", str(dest)])
    run(["tabix", "-f", "-p", "vcf", str(dest)])
    return dest


def _mt_haplogroup(sample: str, vcf: Path, out_dir: Path) -> str:
    mt_vcf = _subset(vcf, config.MT_CONTIG, out_dir / f"{sample}.MT.vcf.gz")
    result = out_dir / f"{sample}.mt.hg"
    run(["java", "-jar", str(config.HAPLOGREP_JAR), "classify",
         "--tree", config.HAPLOGREP_TREE, "--in", str(mt_vcf), "--out", str(result)],
        env=_java_env(),
        stdout=open(config.LOGS_DIR / f"ancestry_{sample}_mt.log", "wb"),
        stderr=subprocess.STDOUT)
    return _parse_haplogrep(sample, result)


def _parse_haplogrep(sample: str, result: Path) -> str:
    """Haplogrep's classification file → the report's mtDNA line, or an honest 'unresolved'."""
    # TSV: SampleID, Haplogroup, Rank, Quality, Range (values are quoted).
    lines = result.read_text().splitlines()
    if len(lines) < 2:
        return "unresolved"
    cells = [c.strip('"') for c in lines[1].split("\t")]
    hg, quality = cells[1], (cells[3] if len(cells) > 3 else "?")
    log.info("%s mtDNA haplogroup: %s (Q=%s)", sample, hg, quality)
    try:
        q = float(quality)
    except ValueError:
        q = 0.0
    if q < config.MT_MIN_QUALITY:
        # Haplogrep always returns *a* haplogroup; with too few informative variants that
        # answer is the tree root, and it looks exactly as confident as a real call. Report
        # the failure instead — a wrong maternal lineage is worse than none.
        log.warning("%s mtDNA haplogroup %s has quality %s (< %s) — reporting unresolved. "
                    "Too few informative MT variants: check the MT calls, and if the sample "
                    "was lifted from another build, check chrM specifically (build-to-build "
                    "chrM liftover is unreliable).", sample, hg, quality, config.MT_MIN_QUALITY)
        return (f"unresolved — haplogrep quality {quality} below {config.MT_MIN_QUALITY} "
                f"(too few informative MT variants; nearest call was {hg})")
    return f"{hg} (quality {quality})"


def _y_haplogroup(sample: str, vcf: Path, out_dir: Path) -> str:
    y_vcf = _subset(vcf, "Y", out_dir / f"{sample}.Y.vcf.gz")
    yl_dir = out_dir / "yleaf"
    run(["uv", "run", "--with", config.YLEAF_SPEC, "--with", "pysam", "--with", "pandas",
         "Yleaf", "-vcf", str(y_vcf), "-rg", "hg38", "-o", str(yl_dir), "-force"],
        stdout=open(config.LOGS_DIR / f"ancestry_{sample}_y.log", "wb"),
        stderr=subprocess.STDOUT)
    return _parse_yleaf(sample, yl_dir / "hg_prediction.hg")


def _parse_yleaf(sample: str, pred: Path) -> str:
    """Yleaf's prediction file → the report's Y-DNA line."""
    lines = pred.read_text().splitlines() if pred.exists() else []
    if len(lines) < 2:
        return "unresolved"
    cells = lines[1].split("\t")  # Sample, Hg, Hg_marker, ...
    hg, qc = cells[1], (cells[5] if len(cells) > 5 else "?")
    log.info("%s Y haplogroup: %s (QC=%s)", sample, hg, qc)
    return f"{hg} (QC {qc})"


def ancestry_report(sample: str, force: bool = False) -> Path:
    require_tools("bcftools", "tabix", "uv")
    if not config.HAPLOGREP_JAR.exists():
        raise SystemExit(f"Haplogrep not installed at {config.HAPLOGREP_JAR}. Run `setup`.")
    vcf = _called_vcf(sample)
    if not vcf.exists():
        raise SystemExit(f"{vcf} missing. Run `call {sample}` first.")

    out_dir = config.ANCESTRY_DIR / sample
    summary = reportpaths.genome_report("ancestry", sample)
    if summary.exists() and not force:
        log.info("%s ancestry already reported (%s) — skipping (use --force).",
                 sample, summary)
        return summary
    out_dir.mkdir(parents=True, exist_ok=True)

    mt = _mt_haplogroup(sample, vcf, out_dir)
    sex = _sex(sample)
    male = sex in ("male", "m", "1")
    y = _y_haplogroup(sample, vcf, out_dir) if male else None
    return write_ancestry(sample, mt, y)


def write_ancestry(sample: str, mt: str, y: str | None) -> Path:
    """Write the haplogroup report. `y` is None for a person without a Y chromosome."""
    out_dir = config.ANCESTRY_DIR / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = reportpaths.genome_report("ancestry", sample)
    male = y is not None
    if not male:
        y = "— (Y haplogroup is paternal-line; not applicable to females)"

    lines = [
        f"# Ancestry — {sample}",
        "",
        f"- **mtDNA haplogroup** (maternal line): **{mt}**  ·  Haplogrep 3 / "
        f"{config.HAPLOGREP_TREE}",
        f"- **Y-DNA haplogroup** (paternal line): **{y}**  ·  Yleaf (hg38)"
        if male else f"- **Y-DNA haplogroup**: {y}",
        "",
        "## Notes",
        "- Haplogroups are deep-lineage markers, not a continental 'ancestry %' "
        "(which this pipeline does not estimate).",
        "- Research-grade.",
        "",
    ]
    write_report(summary, lines)
    log.info("Ancestry report: %s", summary)
    return summary
