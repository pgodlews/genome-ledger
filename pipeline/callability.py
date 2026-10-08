"""Stage: callability / coverage report via mosdepth.

A consumer report implies completeness; ours should be honest about what the sequencing
could and could not see. SNV calling emits only differences from the reference, so a gene
absent from the report is ambiguous — truly reference, or never adequately covered?
mosdepth measures read depth genome-wide and per key gene, so "variant absent" can be told
apart from "region not callable". Complements the SV stage (C1): SV + coverage together =
an honest account of what we could and couldn't assess.

mosdepth is a single static binary (no Python 2 / Docker), so this runs anywhere the CRAM +
reference live. The per-gene table reuses the SV gene panel (the clinically critical genes
where coverage matters most).
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path

from . import reportpaths
from . import config
from .call import aligned_cram
from .util import log, read_tsv, write_report


def _mosdepth_bin() -> str | None:
    hits = list(config.MOSDEPTH_DIR.glob("mosdepth"))
    return str(hits[0]) if hits else shutil.which("mosdepth")


def _write_panel_bed(panel: list[dict], bed: Path) -> None:
    """Panel TSV coords are 1-based inclusive; BED is 0-based half-open."""
    bed.write_text("".join(
        f"{g['chrom']}\t{int(g['start']) - 1}\t{int(g['end'])}\t{g['gene']}\n"
        for g in panel))


def _mean_depth(summary: Path) -> float | None:
    """Genome-wide mean depth from mosdepth's summary (the 'total' row, col 'mean')."""
    for line in summary.read_text().splitlines():
        f = line.split("\t")
        if f[0] == "total":
            return float(f[3])
    return None


def _global_callable_frac(dist: Path, min_depth: int) -> float | None:
    """Fraction of the genome covered at >= min_depth (mosdepth global dist: the cumulative
    proportion of bases at each depth; the 'total' row aggregates all contigs)."""
    for line in dist.read_text().splitlines():
        chrom, depth, prop = line.split("\t")
        if chrom == "total" and int(depth) == min_depth:
            return float(prop)
    return None


def _gene_callability(thresholds_gz: Path, regions_gz: Path, min_depth: int) -> list[dict]:
    """Per-gene callable fraction (bases >= min_depth / gene length) + mean depth."""
    means: dict[str, float] = {}
    with gzip.open(regions_gz, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:                      # chrom start end region mean
                means[f[3]] = float(f[4])
    rows = []
    with gzip.open(thresholds_gz, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            ci = header.index(f"{min_depth}X")    # e.g. "10X" column
        except ValueError:
            ci = None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom, start, end, region = f[0], int(f[1]), int(f[2]), f[3]
            length = end - start
            ge = int(f[ci]) if ci is not None and ci < len(f) else 0
            rows.append({
                "gene": region, "chrom": chrom,
                "mean": means.get(region),
                "frac": (ge / length) if length else 0.0,
            })
    return rows


def callability_report(sample: str, force: bool = False) -> Path:
    mosdepth = _mosdepth_bin()
    if not mosdepth:
        raise SystemExit(f"mosdepth not found (under {config.MOSDEPTH_DIR} or PATH). "
                         "Run `setup` (or `brew install mosdepth` on macOS).")
    cram = aligned_cram(sample)
    if not cram.exists():
        raise SystemExit(f"{cram} missing. Run `align {sample}` first.")

    out_dir = config.CALLABLE_DIR / sample
    summary_md = reportpaths.genome_report("callable", sample)
    if summary_md.exists() and not force:
        log.info("%s callability already reported (%s) — skipping (use --force).",
                 sample, summary_md)
        return summary_md
    out_dir.mkdir(parents=True, exist_ok=True)

    panel = read_tsv(config.SV_GENES_FILE)
    bed = out_dir / "panel.bed"
    _write_panel_bed(panel, bed)
    prefix = out_dir / sample
    md = config.CALLABLE_MIN_DEPTH

    summary = Path(f"{prefix}.mosdepth.summary.txt")
    if force or not summary.exists():
        log.info("Running mosdepth on %s …", sample)
        with open(config.LOGS_DIR / f"callable_{sample}.log", "wb") as run_log:
            subprocess.run(
                [mosdepth, "-t", "4", "--no-per-base", "--by", str(bed),
                 "--thresholds", f"1,{md},20", "--fasta", str(config.REF_FASTA),
                 str(prefix), str(cram)],
                check=True, stdout=run_log, stderr=subprocess.STDOUT,
                cwd=str(out_dir))

    mean = _mean_depth(summary)
    gcall = _global_callable_frac(Path(f"{prefix}.mosdepth.global.dist.txt"), md)
    genes = _gene_callability(Path(f"{prefix}.thresholds.bed.gz"),
                              Path(f"{prefix}.regions.bed.gz"), md)
    disease = {g["gene"]: g.get("disease", "") for g in panel}
    genes.sort(key=lambda r: r["frac"])          # worst-covered first
    low = [g for g in genes if g["frac"] < config.CALLABLE_MIN_FRAC]

    def _pct(x: float | None) -> str:
        return f"{x * 100:.1f}%" if x is not None else "?"

    def _row(g: dict) -> str:
        status = "✅ callable" if g["frac"] >= config.CALLABLE_MIN_FRAC else "⚠️ low coverage"
        mean_s = f"{g['mean']:.0f}x" if g["mean"] is not None else "?"
        return (f"| {g['gene']} | {disease.get(g['gene'], '')} | {mean_s} | "
                f"{_pct(g['frac'])} | {status} |")

    lines = [
        f"# Callability & coverage — {sample}",
        "",
        f"- **Tool:** mosdepth {config.MOSDEPTH_VERSION}",
        f"- **Genome mean depth:** {mean:.1f}x" if mean is not None else "- **Genome mean depth:** ?",
        f"- **Genome callable (≥{md}x):** {_pct(gcall)}",
        f"- **Key genes below {int(config.CALLABLE_MIN_FRAC * 100)}% callable:** {len(low)}",
        "",
        f"## Key-gene callability (fraction of gene ≥{md}x)",
        "",
        "| Gene | Disease | Mean depth | Callable | Status |",
        "|------|---------|-----------|----------|--------|",
        *[_row(g) for g in genes],
        "",
        "## Caveats",
        f"- **\"Callable\" = ≥{md}x.** Below that, a \"no variant\" result for that base is "
        "not reliable — absence of a finding there is uninformative, not reassuring.",
        "- Segmental-duplication genes (SMN1, PMS2, CYP21A2) can show inflated or ambiguous "
        "depth from multi-mapping reads — interpret their coverage cautiously.",
        "- Genome-wide callable % includes hard-to-map regions (telomeres, centromeres) that "
        "depress the number without affecting most genes.",
        "- Research-grade, not a clinical diagnosis.",
        "",
    ]
    write_report(summary_md, lines)
    log.info("Callability report: %s (mean %.1fx, %d low-coverage genes)",
             summary_md, mean or 0.0, len(low))
    return summary_md
