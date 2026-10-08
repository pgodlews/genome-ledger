"""Stage 3: normalize a sample once (split multiallelics, left-align, index).

Normalization is what makes annotations comparable across snapshots: every run
sees the same canonical variant representation, so a later diff reflects real
database changes rather than representational drift.
"""

from __future__ import annotations

from . import config
from .util import find_row, log, require_tools, run

PASS_ONLY = True  # keep only FILTER=PASS records


def normalized_path(sample: str):
    return config.NORMALIZED_DIR / f"{sample}.norm.vcf.gz"


def normalize(sample: str, force: bool = False) -> None:
    require_tools("bcftools", "tabix")
    row = find_row(config.SAMPLES_TSV, sample_id=sample)
    if not row:
        raise SystemExit(f"Unknown sample {sample}. Ingest it first.")

    raw = row["raw_path"]
    out = normalized_path(sample)
    out.parent.mkdir(parents=True, exist_ok=True)
    tbi = out.with_suffix(".gz.tbi")
    # Require the index too: a killed run can leave a truncated .vcf.gz (or one whose
    # tabix never ran), and an exists()-only check would then treat it as complete.
    if out.exists() and tbi.exists() and not force:
        log.info("%s already normalized — skipping (use --force to redo).", sample)
        return

    if not config.REF_FASTA.exists():
        raise SystemExit(f"Reference FASTA missing: {config.REF_FASTA}. Run `setup`.")

    rename_file = _write_rename_map()

    log.info("Normalizing %s …", sample)
    # Pipe: PASS filter → rename contigs (chr1→1) → split multiallelics +
    # left-align against the reference → bgzip.
    view = ["bcftools", "view"]
    if PASS_ONLY:
        # "PASS,." keeps unfiltered records too. Our own calls are hard-filtered so every
        # record carries PASS or a filter name, but a legacy/vendor VCF can have FILTER="."
        # throughout — against a bare "-f PASS" that silently normalizes to an EMPTY file,
        # and the emptiness only surfaces much later as a sample with no variants.
        view += ["-f", "PASS,."]
    view += [raw, "-Ou"]
    # -x ID alongside the contig rename: the DuckDB loader parses VEP's Uploaded_variation
    # as "chrom_pos_ref/alt", but VEP emits the *supplied* identifier there whenever the
    # input VCF's ID column is populated ("rs334"), and load.py's `WHERE pos IS NOT NULL`
    # then dropped those rows without a word. Clearing ID here makes VEP synthesize the
    # coordinate form for every record, so the loader's assumption holds by construction.
    # Nothing downstream reads ID from the normalized VCF (ClinVar identity comes from the
    # snapshot's own ALLELEID), so this loses no information.
    rename = ["bcftools", "annotate", "--rename-chrs", str(rename_file), "-x", "ID", "-Ou"]
    norm = ["bcftools", "norm", "-f", str(config.REF_FASTA), "-m", "-both",
            "-Oz", "-o", str(out)]

    import subprocess
    try:
        p1 = subprocess.Popen([str(c) for c in view], stdout=subprocess.PIPE)
        p2 = subprocess.Popen([str(c) for c in rename], stdin=p1.stdout,
                              stdout=subprocess.PIPE)
        p1.stdout.close()
        subprocess.run([str(c) for c in norm], stdin=p2.stdout, check=True)
        p2.stdout.close()
        for p, name in ((p1, "view"), (p2, "annotate")):
            p.wait()
            if p.returncode:
                raise SystemExit(f"bcftools {name} failed for {sample}")

        # -f, like every other tabix call in the pipeline: without it, re-indexing an
        # existing .tbi is a hard error ("the index file exists"), and the cleanup below
        # then deletes the VCF that was just built. Because the skip-guard above only lets
        # us reach here with a stale index when --force was passed, `normalize --force`
        # destroyed the sample's normalized VCF every time — and left the .tbi behind, so
        # each subsequent plain `normalize` rebuilt the VCF and deleted it again.
        run(["tabix", "-f", "-p", "vcf", str(out)])
    except BaseException:
        # BaseException, not Exception: the pipe-stage check above raises SystemExit, which
        # derives from BaseException — an `except Exception` here let the partial output
        # survive the very failure this cleanup exists for. Also covers KeyboardInterrupt,
        # where deleting the half-written VCF is equally what we want.
        out.unlink(missing_ok=True)  # don't leave a truncated VCF to be skipped-to later
        tbi.unlink(missing_ok=True)  # ...nor an index pointing at a VCF that no longer exists
        raise
    log.info("Wrote %s (+ .tbi). Next: annotate against a snapshot.", out)


def _write_rename_map():
    """Write the chr→Ensembl contig map (idempotent) for bcftools annotate."""
    path = config.NORMALIZED_DIR / "contig_rename.txt"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{old}\t{new}\n"
                            for old, new in config.CONTIG_RENAME.items()))
    return path
