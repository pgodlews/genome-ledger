"""Stage: pharmacogenomics report via PharmCAT.

PharmCAT turns a sample's variants into CPIC/DPWG star-allele diplotypes, metabolizer
phenotypes, and guideline drug guidance. It is a Java jar + a Python preprocessor, so we
rename our Ensembl-named VCF to the UCSC (chr1) naming PharmCAT expects, run its pipeline
(its Python deps supplied ephemerally via `uv run --with`, kept out of this project's
env), and summarize the JSON into markdown alongside PharmCAT's own HTML.

PGx depends on the genome + PharmCAT version (not ClinVar), so the report is filed per
genome, outside the dated snapshots; PharmCAT's own files stay in PGX_DIR.

Caveats surfaced in the report:
  - WGS VCFs carry only variant sites, so PGx positions absent from the VCF are assumed
    reference (`--absent-to-ref`) — fine for high-coverage WGS; the rigorous fix is a
    force-call step at the PGx positions.
  - CYP2D6 is not callable from a plain VCF (CNV/hybrid/pseudogene) → "No Result"; a
    Cyrius/Aldy outside-call off the CRAM is the follow-up.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from . import reportpaths
from . import config
from .call import _java_env, aligned_cram  # openjdk@17 on PATH; CRAM path for Cyrius
from .util import log, require_tools, run, write_report

# Ensembl (1, MT) → UCSC (chr1, chrM) for PharmCAT's chr-prefixed reference.
_ENS2CHR = {**{str(i): f"chr{i}" for i in range(1, 23)},
            "X": "chrX", "Y": "chrY", "MT": "chrM"}
# Phenotype substrings that are NOT actionable (used to flag the rest).
_NORMAL_HINTS = ("normal", "indeterminate", "no result", "n/a", "—", "favorable",
                 "reference", "non-responsive")


def _with(specs) -> list[str]:
    """Expand pinned dependency specs into repeated `--with` arguments for `uv run`."""
    out = []
    for spec in specs:
        out += ["--with", spec]
    return out


def _called_vcf(sample: str) -> Path:
    """This sample's VCF — the pipeline's own call, or the normalized ingested one."""
    from .call import sample_vcf

    return sample_vcf(sample)


def _cyrius_cyp2d6(sample: str, out_dir: Path) -> tuple[str | None, str]:
    """Genotype CYP2D6 from the CRAM with Cyrius. Returns (diplotype, note). diplotype is
    None when Cyrius/the CRAM is unavailable or the call is empty — the report then keeps
    the standard "CYP2D6 unresolved" caveat. Best-effort: never breaks the PGx report."""
    star = config.CYRIUS_DIR / "star_caller.py"
    cram = aligned_cram(sample)
    if not star.exists():
        return None, "Cyrius not installed (run `setup`)"
    if not cram.exists():
        return None, f"CRAM {cram.name} not present on this host"
    manifest = out_dir / "cyrius_manifest.txt"
    manifest.write_text(f"{cram}\n")
    log.info("Running Cyrius (CYP2D6) on %s …", sample)
    subprocess.run(
        ["uv", "run", *_with(config.PGX_STAR_DEPS), "python", str(star),
         "--manifest", str(manifest), "--genome", "38",
         "--reference", str(config.REF_FASTA), "--prefix", sample,
         "--outDir", str(out_dir), "--threads", "4"],
        check=True, cwd=str(config.CYRIUS_DIR),
        stdout=open(config.LOGS_DIR / f"cyrius_{sample}.log", "wb"),
        stderr=subprocess.STDOUT)
    tsv = out_dir / f"{sample}.tsv"
    rows = [ln.split("\t") for ln in tsv.read_text().splitlines() if ln and ln[0] != "#"]
    data = [r for r in rows if r and r[0] != "Sample"]   # skip header
    if not data or len(data[0]) < 2:
        return None, "Cyrius produced no genotype"
    diplotype, filt = data[0][1], (data[0][2] if len(data[0]) > 2 else "?")
    if diplotype in ("", "None", ".", "no_call"):
        return None, f"Cyrius could not resolve CYP2D6 (filter {filt})"
    return diplotype, f"Cyrius {diplotype} (filter {filt})"


def pgx_report(sample: str, force: bool = False) -> Path:
    require_tools("bcftools", "tabix", "uv")
    if not (config.PHARMCAT_JAR.exists() and config.PHARMCAT_PIPELINE.exists()):
        raise SystemExit(f"PharmCAT not installed under {config.PHARMCAT_DIR}. Run `setup`.")
    vcf = _called_vcf(sample)
    if not vcf.exists():
        raise SystemExit(f"{vcf} missing. Run `call {sample}` first.")

    out_dir = config.PGX_DIR / sample
    summary = reportpaths.genome_report("pgx", sample)
    if summary.exists() and not force:
        log.info("%s PGx already reported (%s) — skipping (use --force).", sample, summary)
        return summary
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. rename contigs Ensembl → UCSC, as PharmCAT expects.
    rename_map = out_dir / "ens2chr.txt"
    rename_map.write_text("".join(f"{e}\t{u}\n" for e, u in _ENS2CHR.items()))
    chr_vcf = out_dir / f"{sample}.chr.vcf.gz"
    run(["bcftools", "annotate", "--rename-chrs", str(rename_map),
         "-Oz", "-o", str(chr_vcf), str(vcf)])
    run(["tabix", "-f", "-p", "vcf", str(chr_vcf)])

    # 2. PharmCAT pipeline (preprocess → match → phenotype → report). Its Python deps
    # come from `uv run --with`; its jar is found relative to the preprocessor cwd.
    log.info("Running PharmCAT %s on %s …", config.PHARMCAT_VERSION, sample)
    run(["uv", "run", *_with(config.PHARMCAT_DEPS), "python", str(config.PHARMCAT_PIPELINE),
         str(chr_vcf), "--absent-to-ref", "-o", str(out_dir),
         "-reporterJson", "-reporterHtml"],
        env=_java_env(), cwd=str(config.PHARMCAT_DIR / "preprocessor"),
        stdout=open(config.LOGS_DIR / f"pgx_{sample}.log", "wb"),
        stderr=subprocess.STDOUT)

    # 3. CYP2D6 via Cyrius (off the CRAM) → PharmCAT outside call. Best-effort: if Cyrius or
    # the CRAM isn't available the report keeps the standard "CYP2D6 unresolved" caveat.
    cyp2d6_note = "not run"
    try:
        diplotype, cyp2d6_note = _cyrius_cyp2d6(sample, out_dir)
    except Exception as e:  # noqa: BLE001 — never let CYP2D6 break the rest of the report
        diplotype, cyp2d6_note = None, f"Cyrius step failed: {e}"
        log.warning("Cyrius CYP2D6 step failed (%s) — PGx report keeps the unresolved caveat.", e)
    if diplotype:
        # Re-run phenotyper + reporter with the outside call merged in (reuses the matcher
        # JSON the pipeline already produced; regenerates report.json/html WITH CYP2D6).
        oc = out_dir / "cyp2d6_outside_call.tsv"
        oc.write_text(f"CYP2D6\t{diplotype}\n")
        log.info("Merging CYP2D6 %s into the PharmCAT report (outside call) …", diplotype)
        run(["java", "-jar", str(config.PHARMCAT_JAR),
             "-phenotyper", "-pi", str(out_dir / f"{sample}.chr.match.json"), "-po", str(oc),
             "-reporter", "-reporterJson", "-reporterHtml",
             "-bf", f"{sample}.chr", "-o", str(out_dir)],
            env=_java_env(), cwd=str(out_dir),  # keep PharmCAT's pharmcat.log out of the repo
            stdout=open(config.LOGS_DIR / f"pgx_cyp2d6_{sample}.log", "wb"),
            stderr=subprocess.STDOUT)

    _write_summary(sample, out_dir / f"{sample}.chr.report.json", summary, cyp2d6_note)
    log.info("PGx report: %s (+ %s.chr.report.html); CYP2D6: %s", summary, sample, cyp2d6_note)
    return summary


def _write_summary(sample: str, report_json: Path, out: Path,
                   cyp2d6_note: str = "not run") -> None:
    rep = json.loads(report_json.read_text())
    genes = rep.get("genes", {})
    rows = []
    for gene, g in sorted(genes.items()):
        dips = g.get("recommendationDiplotypes") or g.get("sourceDiplotypes") or []
        if dips:
            call = dips[0].get("label") or "?"
            pheno = "; ".join(dips[0].get("phenotypes") or []) or "—"
        else:
            call, pheno = "(no call)", "—"
        rows.append((gene, call, pheno))

    flagged = [r for r in rows if not any(h in r[2].lower() for h in _NORMAL_HINTS)]
    flagged_md = [f"| {g} | {c} | {p} |" for g, c, p in flagged] or ["| _none_ | | |"]
    all_md = [f"| {g} | {c} | {p} |" for g, c, p in rows]

    lines = [
        f"# Pharmacogenomics — {sample}",
        "",
        f"- **Tool:** PharmCAT {rep.get('pharmcatVersion', config.PHARMCAT_VERSION)} "
        f"(data {rep.get('dataVersion', '?')})  ·  **Genes called:** {len(rows)}",
        f"- **Non-normal metabolizer / function calls:** {len(flagged)}",
        "",
        "## Actionable / non-normal calls",
        "",
        "| Gene | Diplotype | Phenotype / function |",
        "|------|-----------|----------------------|",
        *flagged_md,
        "",
        "## All gene calls",
        "",
        "| Gene | Diplotype | Phenotype / function |",
        "|------|-----------|----------------------|",
        *all_md,
        "",
        "## Caveats",
        "- PGx positions absent from the WGS VCF were assumed reference "
        "(`--absent-to-ref`); rigorous calling needs a force-call step at the PGx sites.",
        # CYP2D6 isn't callable from a plain VCF (CNV/hybrid/pseudogene); we resolve it from
        # the CRAM with Cyrius and merge it as a PharmCAT outside call when available.
        (f"- **CYP2D6** resolved from the CRAM via {cyp2d6_note} and merged as a PharmCAT "
         "outside call."
         if cyp2d6_note.startswith("Cyrius ") else
         f"- **CYP2D6** not resolvable from a plain VCF (CNV/hybrid/pseudogene); Cyrius "
         f"outside-call unavailable here ({cyp2d6_note})."),
        f"- Research-grade. Full CPIC/DPWG drug guidance is in `{sample}.chr.report.html`; "
        "confirm clinically before any prescribing change.",
        "",
    ]
    write_report(out, lines)
