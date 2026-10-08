"""Stage: structural-variant / CNV calling via Manta.

SNV/indel calling is blind to large deletions and duplications — a whole exon of BRCA1
missing, the SMN1 exon-7 deletion behind spinal muscular atrophy, DMD exon-level CNVs.
The monogenic panel flags exactly these genes as "not assessable from SNV calls"; this
stage starts to close that gap.

Manta calls germline SVs (DEL/DUP/INS/INV/BND) straight from the CRAM, but its workflow
scripts (`configManta.py`/`runWorkflow.py`) need Python 2 — so we run it inside a Docker
container (biocontainers), which also means it runs natively on the Linux box where Docker
and the CRAMs live. This stage runs Manta, then interprets the PASS calls that overlap the
curated SV-relevant gene panel (`panels/sv_genes.tsv`).
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
from pathlib import Path

from . import reportpaths
from . import config
from .call import aligned_cram
from .util import log, read_tsv, require_tools, run, write_report

# SV types Manta emits that have a genomic span we can intersect with genes. BND
# (translocation breakends) have no END and are reported separately as a count.
_SPANNED = ("DEL", "DUP", "INV", "INS")

# A focal panel-gene hit must be smaller than this. Manta also emits chromosome-arm-scale
# events (tens of Mb) that span a gene incidentally — those aren't focal findings and are
# usually artifacts or whole-arm CNVs, so they're set aside, not reported as a gene hit.
# 3 Mb still preserves a whole-gene deletion of the largest panel gene (DMD, ~2.2 Mb).
_FOCAL_MAX_BP = 3_000_000


def _info(field: str, info: str) -> str:
    for kv in info.split(";"):
        if kv == field:           # flag-style INFO key (e.g. IMPRECISE)
            return "true"
        if kv.startswith(field + "="):
            return kv.split("=", 1)[1]
    return ""


def _genotype(fmt: str, sample: str) -> str:
    """The GT cell, phase-normalized to '/' separators.

    Manta emits unphased genotypes, so this is defensive — but the carried-check below
    compares against literal "0/0" and "./.", which a phased "0|0" would slip past and be
    counted as a carried SV."""
    keys, vals = fmt.split(":"), sample.split(":")
    try:
        return vals[keys.index("GT")].replace("|", "/")
    except (ValueError, IndexError):
        return "."


def _parse_sv_vcf(vcf_gz: Path) -> tuple[list[dict], int]:
    """Return (spanned SV records, breakend count). Only PASS calls are kept."""
    svs: list[dict] = []
    bnd = 0
    with gzip.open(vcf_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, pos, _id, _ref, _alt, _qual, filt, info, fmt, sample = f[:10]
            if filt not in ("PASS", "."):
                continue
            svtype = _info("SVTYPE", info)
            gt = _genotype(fmt, sample)
            if gt in ("0/0", "./.", "."):   # not carried / no genotype
                continue
            if svtype == "BND":
                bnd += 1
                continue
            if svtype not in _SPANNED:
                continue
            start = int(pos)
            end_s = _info("END", info)
            end = int(end_s) if end_s else start
            svlen = _info("SVLEN", info)
            svs.append({
                "chrom": chrom, "start": start, "end": max(start, end),
                # VCF POS/END are 1-based inclusive, so the span is end-start+1.
                "svtype": svtype, "svlen": abs(int(svlen)) if svlen.lstrip("-").isdigit()
                else (max(start, end) - start + 1), "gt": gt,
            })
    return svs, bnd


def _overlapping_genes(sv: dict, panel: list[dict]) -> list[dict]:
    """Panel genes whose span intersects this SV (same contig, ranges overlap)."""
    hits = []
    for g in panel:
        if g["chrom"] != sv["chrom"]:
            continue
        if sv["start"] <= int(g["end"]) and int(g["start"]) <= sv["end"]:
            hits.append(g)
    return hits


def _run_manta(cram: Path, run_dir: Path) -> Path:
    """Run Manta in Docker; return the diploidSV.vcf.gz path. Mounts GENOMES_ROOT at the
    same path inside the container so all of CRAM/reference/output resolve unchanged."""
    require_tools("docker")
    root = config.GENOMES_ROOT
    uid, gid = os.getuid(), os.getgid()
    threads = str(os.cpu_count() or 4)
    inner = (
        f"configManta.py --bam {cram} --referenceFasta {config.REF_FASTA} "
        f"--runDir {run_dir} && {run_dir}/runWorkflow.py -m local -j {threads}"
    )
    run(["docker", "run", "--rm",
         "-v", f"{root}:{root}",
         "-w", str(run_dir),
         "-e", f"HOME={run_dir}",
         "--user", f"{uid}:{gid}",
         config.MANTA_IMAGE, "bash", "-c", inner],
        stdout=open(config.LOGS_DIR / f"sv_{cram.stem}.log", "wb"),
        stderr=subprocess.STDOUT)
    return run_dir / "results" / "variants" / "diploidSV.vcf.gz"


def sv_report(sample: str, force: bool = False) -> Path:
    cram = aligned_cram(sample)
    if not cram.exists():
        raise SystemExit(f"{cram} missing. Run `align {sample}` first.")

    out_dir = config.SV_DIR / sample
    summary = reportpaths.genome_report("sv", sample)
    if summary.exists() and not force:
        log.info("%s SV already reported (%s) — skipping (use --force).", sample, summary)
        return summary
    out_dir.mkdir(parents=True, exist_ok=True)

    run_dir = out_dir / "manta"
    diploid = run_dir / "results" / "variants" / "diploidSV.vcf.gz"
    if force or not diploid.exists():
        if run_dir.exists():
            shutil.rmtree(run_dir)        # Manta refuses to reuse a populated runDir
        run_dir.mkdir(parents=True, exist_ok=True)
        log.info("Running Manta (Docker: %s) on %s …", config.MANTA_IMAGE, sample)
        diploid = _run_manta(cram, run_dir)

    panel = read_tsv(config.SV_GENES_FILE)
    svs, bnd = _parse_sv_vcf(diploid)

    flagged = []
    large_spanning = 0      # arm-scale events overlapping a panel gene — set aside, not focal
    for sv in svs:
        genes = _overlapping_genes(sv, panel)
        if not genes:
            continue
        if sv["svlen"] > _FOCAL_MAX_BP:
            large_spanning += 1
            continue
        flagged.append((sv, genes))
    flagged.sort(key=lambda fg: (fg[1][0]["gene"], fg[0]["start"]))

    by_type: dict[str, int] = {}
    for sv in svs:
        by_type[sv["svtype"]] = by_type.get(sv["svtype"], 0) + 1
    type_summary = ", ".join(f"{k}: {v}" for k, v in sorted(by_type.items())) or "none"

    def _kb(n: int) -> str:
        return f"{n/1000:.1f} kb" if n >= 1000 else f"{n} bp"

    def _flag_row(fg) -> str:
        sv, genes = fg
        gene_s = ", ".join(g["gene"] for g in genes)
        disease = "; ".join(sorted({g["disease"] for g in genes}))
        return (f"| {gene_s} | {disease} | {sv['svtype']} | "
                f"{sv['chrom']}:{sv['start']}-{sv['end']} | {_kb(sv['svlen'])} | {sv['gt']} |")

    lines = [
        f"# Structural variants — {sample}",
        "",
        f"- **Caller:** Manta {config.MANTA_VERSION} (Docker)",
        f"- **PASS SVs (carried):** {len(svs)} ({type_summary}); breakends: {bnd}",
        f"- **Focal panel-gene hits:** {len(flagged)}"
        + (f" ({large_spanning} arm-scale event(s) spanning a panel gene set aside as "
           f"non-focal — likely whole-arm CNV/artifact)" if large_spanning else ""),
        "",
        "## Panel-gene hits (review)",
        "",
    ]
    if flagged:
        lines += [
            "| Gene | Disease | Type | Locus (GRCh38) | Size | GT |",
            "|------|---------|------|----------------|------|----|",
            *[_flag_row(fg) for fg in flagged],
        ]
    else:
        lines.append("_No PASS structural variant overlaps a curated SV-panel gene._")
    lines += [
        "",
        "## Caveats",
        "- **Screening, not diagnostic.** Short-read SV calls (esp. in segmental-duplication "
        "regions like SMN1, PMS2, CYP21A2) are error-prone — confirm any flag with MLPA / "
        "array-CGH / long-read sequencing.",
        "- Manta calls DEL/DUP/INS/INV/BND; copy-number-neutral events and exact CNV dosage "
        "need a dedicated CNV caller (a possible follow-up).",
        "- Panel-gene spans are whole-gene-body approximations; an overlap means \"review\", "
        "not \"pathogenic\".",
        "- Research-grade, not a clinical diagnosis.",
        "",
    ]
    write_report(summary, lines)
    log.info("SV report: %s (%d PASS SVs, %d panel hits)", summary, len(svs), len(flagged))
    return summary
