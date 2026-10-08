"""Stage 0: FASTQ → VCF front-end (align + call), GRCh38.

The rest of the pipeline starts at a per-sample VCF (`ingest`). This stage produces
that VCF from the sequencer's raw reads, so the genome is owned end-to-end rather
than trusting a vendor's calls:

    incoming/<vendor>_{1,2}.fq.gz
        │   align  — fastp (trim) │ bwa-mem │ fixmate │ sort │ markdup → CRAM
        ▼
    aligned/<sample>.cram
        │   call   — GATK HaplotypeCaller, scattered by contig, then concat
        ▼
    called/<sample>.GRCh38.vcf.gz   →  ingest → normalize → annotate → …

Both subcommands are idempotent (skip when their output exists) and resumable: a
half-finished `call` reuses the per-contig parts it already produced. Everything is
CPU-bound and streamed — no intermediate SAM/FASTQ is ever written to disk.
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from . import config
from .util import find_row, log, require_tools, run


# --- paths -------------------------------------------------------------------
def aligned_cram(sample: str) -> Path:
    return config.ALIGNED_DIR / f"{sample}.cram"


def called_vcf(sample: str) -> Path:
    return config.CALLED_DIR / f"{sample}.{config.BUILD}.vcf.gz"


def sample_vcf(sample: str) -> Path:
    """The per-sample VCF the genotype-lookup engines (traits/HLA/PRS/ancestry) should read.

    `called_vcf` is this pipeline's own FASTQ→VCF output, which only exists for samples it
    called itself. A sample **ingested** as a finished VCF (a vendor callset, another lab's
    DeepVariant run) never has one, and those engines would refuse to run at all.

    The fallback is the *normalized* VCF, not the ingested raw file, and that matters: the
    raw file carries whatever contig naming its source used (the lifted DeepVariant callsets
    are chr-prefixed), while every one of these engines looks positions up by Ensembl-style
    coordinate. Reading the raw file would silently match nothing and — because they all
    treat an absent position as homozygous reference — report confident wrong genotypes
    genome-wide. `normalize` has already renamed contigs, split multiallelics and
    left-aligned, so it is the representation these lookups actually assume.

    One difference to know about: the normalized VCF is PASS-filtered, so a variant the
    caller filtered out reads as hom-ref rather than as a call. That is the same
    absent→hom-ref convention these engines already apply, applied to a slightly larger set.
    """
    from .normalize import normalized_path

    own = called_vcf(sample)
    if own.exists():
        return own
    norm = normalized_path(sample)
    if norm.exists():
        log.info("%s has no called VCF (ingested sample) — reading %s.", sample, norm.name)
        return norm
    return own          # neither exists: callers report the missing called-VCF path


def _qc_dir() -> Path:
    d = config.ALIGNED_DIR / "qc"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _tmp_dir() -> Path:
    d = config.ALIGNED_DIR / ".tmp"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _incoming_row(sample: str) -> dict:
    row = find_row(config.INCOMING_MANIFEST, sample=sample)
    if not row:
        raise SystemExit(
            f"No entry for '{sample}' in {config.INCOMING_MANIFEST}. Add a row: "
            f"sample\\tvendor_id\\trelation\\tsex"
        )
    return row


# What follows the manifest's vendor_id in a FASTQ file name: an optional run of `_`/`.`
# fields (sample number, lane), then the read number, then the extension. Covers
#   ACME123_1.fq.gz                     (MGI/DNBSEQ and many providers)
#   ACME123_R1.fastq.gz
#   ACME123_S1_L001_R1_001.fastq.gz     (Illumina bcl2fastq/BCL Convert, one file per lane)
#   ACME123_L01_1.fq.gz                 (MGI, one file per lane)
# Index reads (`_I1_001`) do not match, and the id must be followed by `_` or `.`, so
# `ACME12` never picks up `ACME123`'s files.
_FASTQ_NAME = re.compile(r"^(?P<lane>(?:[_.][^/]*?)?)[_.]R?(?P<read>[12])(?P<tail>(?:_\d{3})?)"
                         r"\.f(?:ast)?q(?:\.gz)?$")


def _fastq_files(vendor_id: str) -> tuple[list[Path], list[Path]]:
    """Every R1 file and its R2 mate for this vendor_id, in lane order.

    It used to accept exactly `<vendor_id>_1.fq.gz` + `_2.fq.gz` — one provider's naming —
    so reads named the Illumina way, or split across lanes as Illumina and MGI both deliver
    them, stopped `align` with "FASTQ not found". Refuses an R1 without its R2 (or the
    reverse) rather than aligning a lane as single-end."""
    by_read: dict[str, dict[str, Path]] = {"1": {}, "2": {}}
    for path in sorted(config.INCOMING_DIR.glob(f"{vendor_id}*")):
        m = _FASTQ_NAME.match(path.name[len(vendor_id):])
        if m and path.is_file():
            by_read[m["read"]][m["lane"] + m["tail"]] = path
    r1, r2 = by_read["1"], by_read["2"]
    if not r1 and not r2:
        raise SystemExit(
            f"FASTQ not found for '{vendor_id}' in {config.INCOMING_DIR}. Expected paired files "
            f"such as {vendor_id}_1.fq.gz / _2.fq.gz, {vendor_id}_R1.fastq.gz / _R2, or "
            f"{vendor_id}_S1_L001_R1_001.fastq.gz / _R2_001 (one pair per lane).")
    unpaired = sorted({k for k in r1.keys() ^ r2.keys()})
    if unpaired:
        names = [p.name for k in unpaired for p in (r1.get(k), r2.get(k)) if p]
        raise SystemExit(f"FASTQ without a mate for '{vendor_id}': {', '.join(names)}.")
    lanes = sorted(r1)
    return [r1[k] for k in lanes], [r2[k] for k in lanes]


def _lane_stream(files: list[Path], fifo: Path):
    """Decompress lane files, in order, into a named pipe that fastp reads as one input.

    fastp takes a single file per mate; concatenating 60-90 GB of lanes on disk first would
    double the footprint for nothing. `gzip -dcf` passes plain-text FASTQ through as-is."""
    os.mkfifo(fifo)
    return _popen(["sh", "-c", 'out=$1; shift; exec gzip -dcf -- "$@" > "$out"', "sh",
                   fifo, *files], stderr=subprocess.DEVNULL)


def _java_env() -> dict:
    """GATK is pure Java but needs a JRE on PATH; openjdk@17 is keg-only, so prepend
    its bin rather than relying on a system Java (there is none on this Mac).

    Homebrew-specific, so it is skipped where brew isn't the JRE source — a Linux
    host often has a system Java on PATH, where `brew --prefix` would raise instead, breaking
    every GATK stage that imports this (force-call, mito, ancestry)."""
    env = dict(os.environ)
    if shutil.which("brew") is None:
        return env
    try:
        prefix = subprocess.check_output(
            ["brew", "--prefix", config.JDK_FORMULA], text=True).strip()
    except subprocess.CalledProcessError:   # formula not installed — fall back to PATH java
        return env
    env["PATH"] = f"{prefix}/bin:{env['PATH']}"
    return env


# --- align -------------------------------------------------------------------
# Read-name formats, matched against the first record of R1. MGI/DNBSEQ: flowcell, then
# L<lane>C<column>R<row> and a running number — "@V300012345L1C001R0010000001/1".
# Illumina (CASAVA 1.8+): seven colon-separated fields — "@A00123:8:HABCDXX:1:1101:1000:2000";
# before 1.8: "@HWUSI-EAS100R:6:73:941:1973#0/1".
_PLATFORM_READ_NAMES = (
    ("DNBSEQ", re.compile(r"^@[A-Z0-9]+L\d{1,2}C\d{3}R\d{3}\d+(/[12])?(\s|$)")),
    ("ILLUMINA", re.compile(r"^@[^\s:]+:\d+:[^\s:]+:\d+:\d+:\d+:\d+(\s|$)")),
    ("ILLUMINA", re.compile(r"^@[^\s:]+:\d+:\d+:\d+:\d+#")),
)


def read_platform(fastq: Path) -> str:
    """The @RG PL value for these reads: READ_PLATFORM if set, else detected from the first
    read name, else ILLUMINA with a warning.

    The tag is not cosmetic. Base-quality recalibration and duplicate marking treat platforms
    differently, so a default that suits one sequencer silently mislabels another's reads.
    Reads re-exported from SRA or renamed by another tool carry neither format; that is the
    case the warning is for."""
    if config.READ_PLATFORM:
        return config.READ_PLATFORM
    opener = gzip.open if str(fastq).endswith(".gz") else open
    try:
        with opener(fastq, "rt") as fh:
            first = fh.readline()
    except (OSError, EOFError, UnicodeDecodeError):
        first = ""
    for platform, pattern in _PLATFORM_READ_NAMES:
        if pattern.match(first):
            log.info("Sequencing platform from the read names in %s: %s.",
                     Path(fastq).name, platform)
            return platform
    log.warning("Could not tell the sequencing platform from the read names in %s (%r) — "
                "tagging the read group PL:ILLUMINA. Set READ_PLATFORM if that is wrong.",
                Path(fastq).name, first.strip()[:60])
    return "ILLUMINA"


def align(sample: str, force: bool = False) -> Path:
    """fastp-trim → bwa-mem → fixmate → coordinate-sort → markdup → CRAM.

    Streamed as one pipe: bwa emits mates adjacently (interleaved `-p` input), so
    `fixmate -m` works directly with no intermediate name-sort, and nothing but the
    final CRAM touches disk."""
    require_tools("fastp", "bwa", "samtools")
    row = _incoming_row(sample)
    r1s, r2s = _fastq_files(row["vendor_id"])
    out = aligned_cram(sample)
    out.parent.mkdir(parents=True, exist_ok=True)
    # Require the index too, as `call`/`normalize`/`force-call` do: samtools index runs
    # after markdup, so a kill in between leaves a complete-looking CRAM with no .crai and
    # an exists()-only check would skip straight past it.
    if out.exists() and Path(f"{out}.crai").exists() and not force:
        log.info("%s already aligned (%s) — skipping (use --force).", sample, out.name)
        return out

    ref = config.REF_FASTA
    if not (ref.with_suffix(".fa.bwt").exists() or Path(str(ref) + ".bwt").exists()):
        raise SystemExit(f"bwa index missing for {ref}. Run `setup`.")

    qc, tmp = _qc_dir(), _tmp_dir()
    logp = config.LOGS_DIR / f"align_{sample}.log"
    platform = read_platform(r1s[0])
    rg = (f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:{platform}"
          f"\\tLB:{sample}")

    writers = []
    if len(r1s) == 1:
        r1, r2 = r1s[0], r2s[0]
    else:
        lane_dir = tmp / f"lanes_{sample}"
        shutil.rmtree(lane_dir, ignore_errors=True)
        lane_dir.mkdir(parents=True)
        r1, r2 = lane_dir / "R1.fastq", lane_dir / "R2.fastq"
        log.info("%s: %d lanes, streamed in order: %s", sample, len(r1s),
                 ", ".join(p.name for p in r1s))
    # Poly-G tails are an artefact of Illumina's two-colour chemistry, and fastp turns the trim
    # on by itself when the read names say NovaSeq/NextSeq. Forcing it for every Illumina run
    # would also clip real G stretches on four-colour instruments. DNBSEQ keeps the forced trim
    # this pipeline has always applied to it, so those alignments stay reproducible.
    poly_g = ["--trim_poly_g"] if platform == "DNBSEQ" else []
    fastp_cmd = ["fastp", "--in1", r1, "--in2", r2, "--stdout",
                 "--thread", "4", "--detect_adapter_for_pe", *poly_g,
                 "--json", qc / f"{sample}.fastp.json",
                 "--html", qc / f"{sample}.fastp.html"]
    # -K 100M: fixed batch size → output independent of thread count (reproducible).
    # -Y: soft-clip supplementary alignments (GATK-friendly). -p: interleaved stdin.
    bwa_cmd = ["bwa", "mem", "-p", "-t", config.ALIGN_THREADS, "-K", "100000000",
               "-Y", "-R", rg, ref, "-"]
    fixmate_cmd = ["samtools", "fixmate", "-m", "-u", "-", "-"]
    sort_cmd = ["samtools", "sort", "-u", "-@", config.SORT_THREADS,
                "-m", config.SORT_MEM, "-T", tmp / f"sort_{sample}", "-"]
    # Pin CRAM 3.0: samtools defaults to 3.1, which GATK's bundled htsjdk can't read.
    markdup_cmd = ["samtools", "markdup", "-@", config.SORT_THREADS,
                   "--reference", ref, "--output-fmt", "cram,version=3.0",
                   "-f", qc / f"{sample}.markdup.txt",
                   "-T", tmp / f"markdup_{sample}", "-", out]

    log.info("Aligning %s (%s, %s, %s) — bwa-mem on %s threads → %s",
             sample, row["vendor_id"], platform, config.BUILD, config.ALIGN_THREADS, out.name)
    with open(logp, "wb") as lf:
        if len(r1s) > 1:
            writers = [_lane_stream(r1s, r1), _lane_stream(r2s, r2)]
        fastp = _popen(fastp_cmd, stdout=subprocess.PIPE, stderr=lf)
        bwa = _popen(bwa_cmd, stdin=fastp.stdout, stdout=subprocess.PIPE, stderr=lf)
        fastp.stdout.close()
        fix = _popen(fixmate_cmd, stdin=bwa.stdout, stdout=subprocess.PIPE, stderr=lf)
        bwa.stdout.close()
        srt = _popen(sort_cmd, stdin=fix.stdout, stdout=subprocess.PIPE, stderr=lf)
        fix.stdout.close()
        mkd = _popen(markdup_cmd, stdin=srt.stdout, stderr=lf)
        srt.stdout.close()
        rcs = {"markdup": mkd.wait(), "sort": srt.wait(), "fixmate": fix.wait(),
               "bwa": bwa.wait(), "fastp": fastp.wait()}
        for i, w in enumerate(writers):
            if fastp.returncode and w.poll() is None:
                w.kill()     # fastp died: nothing will ever read the rest of the pipe
            rcs[f"lane reader R{i + 1}"] = w.wait()
    if writers:
        shutil.rmtree(r1.parent, ignore_errors=True)
    failed = [name for name, rc in rcs.items() if rc]
    if failed:
        out.unlink(missing_ok=True)
        raise SystemExit(f"align {sample} failed in: {', '.join(failed)} "
                         f"(see {logp})")

    run(["samtools", "index", str(out)])
    flag = qc / f"{sample}.flagstat.txt"
    with open(flag, "w") as fh:
        run(["samtools", "flagstat", str(out)], stdout=fh)
    log.info("Wrote %s (+ .crai). QC in %s. Next: `call %s`.", out, qc, sample)
    return out


# --- scatter intervals -------------------------------------------------------
def _fai_lengths() -> dict[str, int]:
    """Contig → length from the reference .fai (col 1 = name, col 2 = length)."""
    fai = config.REF_FASTA.with_suffix(".fa.fai")
    out: dict[str, int] = {}
    for line in fai.read_text().splitlines():
        name, length = line.split("\t")[:2]
        out[name] = int(length)
    return out


def _scatter_tasks() -> list[tuple[str, int, str]]:
    """Sub-chromosomal windows for parallel calling: (region, ploidy, label).

    Diploid contigs are split into ~CALL_WINDOW_BP windows so a slow chromosome
    (chr16's segmental duplications) doesn't gate the whole stage on one core; MT
    stays one haploid shard. Windows are emitted in genomic order so the concatenated
    VCF is coordinate-sorted, and the label encodes zero-padded coordinates so the
    per-shard part files also sort naturally."""
    lengths = _fai_lengths()
    win = config.CALL_WINDOW_BP
    tasks: list[tuple[str, int, str]] = []
    for contig in config.CALL_CONTIGS:
        n = lengths[contig]
        start = 1
        while start <= n:
            end = min(start + win - 1, n)
            tasks.append((f"{contig}:{start}-{end}", 2,
                          f"{contig}_{start:09d}-{end:09d}"))
            start = end + 1
    tasks.append((config.MT_CONTIG, config.MT_PLOIDY, config.MT_CONTIG))  # whole, haploid
    return tasks


def _shard_costs(labels: list[str]) -> dict[str, float]:
    """label -> longest runtime (s) that window has taken on any previously called sample.

    Read from the per-shard logs `_call_shard` already writes: the file is created when the
    shard starts and last written when it ends, so birth->mtime is its wall time. Samples
    sequenced at similar depth with the same library prep take similar time per window, so one
    sample's timings predict the next one's well enough to schedule by.
    """
    seen: list[tuple[str, float]] = []
    for path in config.LOGS_DIR.glob("call_*.log"):
        try:
            st = path.stat()
        except OSError:            # log rotated or deleted mid-scan
            continue
        secs = st.st_mtime - getattr(st, "st_birthtime", st.st_ctime)
        if secs > 0:
            seen.append((path.stem, secs))
    # The filename is call_<sample>_<label>.log and *both* parts can contain "_"
    # (`SAMPLE_t2t`, `16_030000001-060000000`), so match on the label suffix rather
    # than trying to split the stem.
    costs: dict[str, float] = {}
    for label in labels:
        suffix = f"_{label}"
        best = max((s for stem, s in seen if stem.endswith(suffix)), default=0.0)
        if best:
            costs[label] = best
    return costs


# --- call --------------------------------------------------------------------
def call(sample: str, force: bool = False) -> Path:
    """GATK HaplotypeCaller scattered across contigs → concatenated per-sample VCF."""
    require_tools("bcftools", "tabix")
    if not config.GATK_JAR.exists():
        raise SystemExit(f"GATK not installed at {config.GATK_JAR}. Run `setup`.")
    cram = aligned_cram(sample)
    if not cram.exists():
        raise SystemExit(f"{cram} missing. Run `align {sample}` first.")
    out = called_vcf(sample)
    out.parent.mkdir(parents=True, exist_ok=True)
    # Require the index too: a killed run can leave a truncated .vcf.gz (or one whose
    # tabix never ran), and an exists()-only check would then treat it as complete.
    if out.exists() and out.with_suffix(".gz.tbi").exists() and not force:
        log.info("%s already called (%s) — skipping (use --force).", sample, out.name)
        return out

    env = _java_env()
    parts_dir = config.CALLED_DIR / f".parts_{sample}"
    parts_dir.mkdir(parents=True, exist_ok=True)
    if force:
        # The shards below resume on existence alone, which is what makes an interrupted run
        # cheap to restart — but it also meant --force re-concatenated the previous run's
        # shards without re-reading the CRAM. Re-aligning and then calling with --force
        # produced the old calls.
        for stale in parts_dir.glob("*"):
            stale.unlink()

    def _call_shard(task: tuple[str, int, str]) -> Path:
        region, ploidy, label = task
        part = parts_dir / f"{sample}.{label}.vcf.gz"
        if part.exists() and part.with_suffix(".gz.tbi").exists():
            return part  # resume: shard already done
        run(["java", f"-Xmx{config.HC_JAVA_MEM}", "-jar", config.GATK_JAR,
             "HaplotypeCaller", "-R", config.REF_FASTA, "-I", cram,
             "-L", region, "--sample-ploidy", ploidy, "-O", part,
             "--native-pair-hmm-threads", "1"],
            env=env,
            stdout=open(config.LOGS_DIR / f"call_{sample}_{label}.log", "wb"),
            stderr=subprocess.STDOUT)
        return part

    tasks = _scatter_tasks()
    # Schedule longest-first. With HC_WORKERS lanes over len(tasks) windows, a window that
    # runs several times the average has to START early or it becomes a single-core tail
    # with every other lane idle. Genomic order picked chr16:30-60Mb (centromere + 16p11.2
    # segdups) 86th on baseline runs, so it ran for hours alone on one core after every other window
    # was done — pure wall-clock loss. Costs come from previous samples' shard logs; with no
    # history every window scores the same and sorted() is stable, so the order stays
    # genomic and this is a no-op.
    costs = _shard_costs([t[2] for t in tasks])
    typical = sorted(costs.values())[len(costs) // 2] if costs else 0.0
    order = sorted(range(len(tasks)), key=lambda i: -costs.get(tasks[i][2], typical))
    log.info("Calling %s — HaplotypeCaller over %d windows (~%dMb each; MT haploid), "
             "%d in parallel%s.", sample, len(tasks),
             config.CALL_WINDOW_BP // 1_000_000, config.HC_WORKERS,
             f"; longest-first from {len(costs)} known window timings" if costs
             else "; genomic order (no timing history yet)")
    # Submission order is by cost, but `parts` must stay genomic for `bcftools concat` to
    # produce a coordinate-sorted VCF — so results are placed back at their task index.
    parts: list[Path] = [None] * len(tasks)   # type: ignore[list-item]
    with ThreadPoolExecutor(max_workers=config.HC_WORKERS) as ex:
        futures = {ex.submit(_call_shard, tasks[i]): i for i in order}
        for fut in as_completed(futures):
            parts[futures[fut]] = fut.result()   # first failure surfaces here

    raw = config.CALLED_DIR / f"{sample}.{config.BUILD}.raw.vcf.gz"
    log.info("Concatenating %d shard VCFs → %s", len(parts), raw.name)
    try:
        run(["bcftools", "concat", "-Oz", "-o", raw, *[str(p) for p in parts]])
        run(["tabix", "-f", "-p", "vcf", str(raw)])

        # Raw HaplotypeCaller calls are unfiltered (FILTER='.'); apply GATK germline
        # hard-filters so FILTER carries PASS/fail — `normalize` keeps only PASS.
        _hard_filter(sample, raw, out, env)
    except Exception:
        # Don't leave truncated outputs behind — the skip check would treat them as done.
        raw.unlink(missing_ok=True)
        out.unlink(missing_ok=True)
        raise
    for f in (raw, raw.with_suffix(".gz.tbi")):
        f.unlink(missing_ok=True)
    _qc_call(sample, out)
    shutil.rmtree(parts_dir, ignore_errors=True)

    # The manifest lookup raises SystemExit when the sample has no row — which used to abort
    # a fully successful multi-hour call at the very last statement, purely to format a hint.
    # A missing row here costs the hint, not the run.
    row = find_row(config.INCOMING_MANIFEST, sample=sample) or {}
    hint = (f" --relation {row['relation']} --sex {row['sex']}"
            if row.get("relation") and row.get("sex") else
            f" --relation <relation> --sex <sex>   # no row for '{sample}' in "
            f"{config.INCOMING_MANIFEST}")
    log.info("Wrote %s. Next: `ingest %s --sample %s%s`.", out, out, sample, hint)
    return out


# GATK germline hard-filter thresholds (single-sample WGS, too little data for VQSR).
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890471
_SNP_FILTERS = [("QD2", "QD < 2.0"), ("FS60", "FS > 60.0"), ("MQ40", "MQ < 40.0"),
                ("MQRankSum-12.5", "MQRankSum < -12.5"),
                ("ReadPosRankSum-8", "ReadPosRankSum < -8.0"), ("SOR3", "SOR > 3.0")]
_INDEL_FILTERS = [("QD2", "QD < 2.0"), ("FS200", "FS > 200.0"),
                  ("ReadPosRankSum-20", "ReadPosRankSum < -20.0"), ("SOR10", "SOR > 10.0")]


def _gatk(args: list, env: dict, logname: str) -> None:
    run(["java", f"-Xmx{config.HC_JAVA_MEM}", "-jar", config.GATK_JAR, *args], env=env,
        stdout=open(config.LOGS_DIR / logname, "wb"), stderr=subprocess.STDOUT)


def _hard_filter(sample: str, raw: Path, out: Path, env: dict) -> None:
    """Split SNPs/indels, apply each type's hard-filters (FILTER←PASS/fail), remerge."""
    d = config.CALLED_DIR / f".filter_{sample}"
    d.mkdir(parents=True, exist_ok=True)
    ref = config.REF_FASTA
    snp, indel = d / "snp.vcf.gz", d / "indel.vcf.gz"
    snp_f, indel_f = d / "snp.filt.vcf.gz", d / "indel.filt.vcf.gz"

    _gatk(["SelectVariants", "-R", ref, "-V", raw, "--select-type-to-include", "SNP",
           "-O", snp], env, f"filter_{sample}_select_snp.log")
    _gatk(["SelectVariants", "-R", ref, "-V", raw,
           "--select-type-to-include", "INDEL", "--select-type-to-include", "MIXED",
           "--select-type-to-include", "MNP", "-O", indel], env,
          f"filter_{sample}_select_indel.log")

    snp_args = ["VariantFiltration", "-R", ref, "-V", snp, "-O", snp_f]
    for name, expr in _SNP_FILTERS:
        snp_args += ["--filter-name", name, "--filter-expression", expr]
    _gatk(snp_args, env, f"filter_{sample}_filter_snp.log")
    indel_args = ["VariantFiltration", "-R", ref, "-V", indel, "-O", indel_f]
    for name, expr in _INDEL_FILTERS:
        indel_args += ["--filter-name", name, "--filter-expression", expr]
    _gatk(indel_args, env, f"filter_{sample}_filter_indel.log")

    # concat -a (inputs are indexed by GATK) merges the two streams coordinate-sorted.
    run(["bcftools", "concat", "-a", "-Oz", "-o", out, snp_f, indel_f])
    run(["tabix", "-f", "-p", "vcf", str(out)])
    log.info("%s hard-filtered → %s (PASS counts in QC stats)", sample, out.name)
    shutil.rmtree(d, ignore_errors=True)


def _qc_call(sample: str, vcf: Path) -> None:
    """Dump bcftools stats and surface Ti/Tv (~2.0–2.1 for a sane WGS callset)."""
    stats = _qc_dir() / f"{sample}.bcftools_stats.txt"
    with open(stats, "w") as fh:
        run(["bcftools", "stats", str(vcf)], stdout=fh)
    # The "TSTV" data line holds ts, tv, ts/tv, ... — column 5 is the ratio.
    for ln in stats.read_text().splitlines():
        if ln.startswith("TSTV"):
            log.info("%s Ti/Tv = %s (full stats: %s)", sample, ln.split("\t")[4], stats.name)
            break


def _popen(cmd, **kw) -> subprocess.Popen:
    return subprocess.Popen([str(c) for c in cmd], **kw)
