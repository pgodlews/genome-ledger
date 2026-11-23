"""Stage 3: normalize a sample once (split multiallelics, left-align, index).

Normalization is what makes annotations comparable across snapshots: every run
sees the same canonical variant representation, so a later diff reflects real
database changes rather than representational drift.
"""

from __future__ import annotations

import re
import subprocess
from collections import Counter

from . import config
from .util import _fai_contigs, find_row, log, require_tools, run

PASS_ONLY = True  # keep only FILTER=PASS records


def normalized_path(sample: str):
    return config.NORMALIZED_DIR / f"{sample}.norm.vcf.gz"


def _reference_targets(sample: str, raw: str) -> list[str]:
    """`bcftools view` arguments that drop records on contigs the reference does not have.

    `bcftools norm -f` aborts on the first such record, so one call on an alt or decoy contig
    (`chr1_KI270706v1_random`, which CONTIG_RENAME rightly leaves alone) used to stop the whole
    sample. Empty unless a record really is on a foreign contig, so an input that normalized
    before gets exactly the same command, and the same output, as it always did. Every drop
    is counted in the log. A file where most records would go is refused instead: that is a
    naming mismatch (RefSeq accessions, say), not a few alt contigs.
    """
    reference = _fai_contigs()
    if not reference:
        return []
    header = subprocess.run(["bcftools", "view", "-h", raw], capture_output=True, text=True,
                            check=True).stdout
    declared = re.findall(r"^##contig=<ID=([^,>]+)", header, flags=re.M)
    if all(config.CONTIG_RENAME.get(c, c) in reference for c in declared):
        return []

    query = ["bcftools", "query", "-f", "%CHROM\n", raw]
    if PASS_ONLY:
        query[2:2] = ["-i", 'FILTER="PASS" || FILTER="."']
    counts: Counter[str] = Counter()
    with subprocess.Popen(query, stdout=subprocess.PIPE, text=True) as proc:
        for line in proc.stdout:
            counts[line.rstrip("\n")] += 1
    if proc.returncode:
        raise SystemExit(f"bcftools query failed for {sample}")

    keep = {c: n for c, n in counts.items() if config.CONTIG_RENAME.get(c, c) in reference}
    drop = {c: n for c, n in counts.items() if c not in keep}
    if not drop:
        return []
    n_keep, n_drop = sum(keep.values()), sum(drop.values())
    if n_drop >= n_keep:
        raise SystemExit(
            f"{sample}: {n_drop:,} of {n_keep + n_drop:,} records are on contigs absent from "
            f"{config.REF_FASTA.name} (e.g. {', '.join(sorted(drop)[:3])}). This looks like a "
            "contig naming or build mismatch, so nothing was normalized.")
    log.warning("%s: dropping %d record(s) on %d contig(s) absent from %s (kept %d): %s",
                sample, n_drop, len(drop), config.REF_FASTA.name, n_keep,
                ", ".join(f"{c}={n}" for c, n in sorted(drop.items())[:10]))
    return ["-t", ",".join(keep)]


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
    view += [*_reference_targets(sample, raw), raw, "-Ou"]
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
