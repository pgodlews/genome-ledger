"""Stage 1: one-time environment setup (tools + reference data + VEP).

Downloads several GB (reference FASTA ~900MB, VEP indexed cache ~25GB) and builds
VEP from source (incl. Bio::DB::HTS). Run it once, interactively, on a good
connection. Every step is idempotent.
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import socket
import subprocess
import urllib.request
import zipfile
from pathlib import Path

from . import config
from .call import _java_env
from .util import ensure_dirs, log, require_tools, run

# VEP is NOT in Homebrew; it is built from source in _install_vep().
# xz provides liblzma, which htslib's CRAM support needs when built from source.
# bwa/fastp: FASTQ→CRAM alignment; openjdk@17: JRE for GATK (keg-only formula).
BREW_PKGS = ["bcftools", "htslib", "samtools", "duckdb", "xz",
             "bwa", "fastp", config.JDK_FORMULA]


def _build_env() -> dict:
    """Env for the from-source htslib build. On macOS, point clang at Homebrew
    headers/libs (lzma.h etc.) which it doesn't search under /opt/homebrew by default;
    on Linux the apt-installed -dev libs are in standard paths, so no flags needed."""
    import platform
    env = dict(os.environ)
    if platform.system() == "Darwin":
        prefix = subprocess.check_output(["brew", "--prefix"], text=True).strip()
        inc, lib = f"{prefix}/include", f"{prefix}/lib"
        env["CPPFLAGS"] = f"-I{inc} {env.get('CPPFLAGS', '')}".strip()
        env["CFLAGS"] = f"-I{inc} {env.get('CFLAGS', '')}".strip()
        env["LDFLAGS"] = f"-L{lib} {env.get('LDFLAGS', '')}".strip()
        env["C_INCLUDE_PATH"] = inc
        env["LIBRARY_PATH"] = lib
    # INSTALL.pl shells out to `unzip` for the Ensembl API; force overwrite so it
    # never prompts (no TTY in background → a prompt would abort the install).
    env["UNZIP"] = "-o"
    env["UNZIPOPT"] = "-o"
    return env


def _download(url: str, dest: Path, label: str) -> None:
    if dest.exists():
        log.info("%s already present (%s) — skipping.", label, dest.name)
        return
    log.info("Downloading %s → %s", label, dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    # Multi-GB downloads over a connection shared with other concurrent jobs occasionally
    # hit a mid-transfer reset (observed: ConnectionResetError on a 2.4GB file) — retry a
    # few times rather than losing the whole download and forcing a manual re-run. Also seen:
    # a silent STALL (socket never errors, just stops delivering bytes) — pass a per-request
    # timeout (NOT socket.setdefaulttimeout, which would mutate global interpreter state).
    last_exc = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as fh:
                shutil.copyfileobj(r, fh)
            break
        except (ConnectionResetError, TimeoutError, OSError) as e:
            last_exc = e
            tmp.unlink(missing_ok=True)
            log.warning("  %s download attempt %d/4 failed (%s)%s", label, attempt + 1, e,
                       " — retrying." if attempt < 3 else "")
    else:
        raise last_exc
    tmp.rename(dest)


# What setup needs already on the PATH on Linux, and the Debian/Ubuntu packages that
# provide it. setup does not install system packages on Linux: that needs root, and which
# package manager to drive is not this script's decision.
LINUX_TOOLS = ("bcftools", "samtools", "tabix", "bgzip", "bwa", "fastp", "java", "perl",
               "git", "cc", "make", "unzip", "curl")
LINUX_APT = ("bcftools samtools tabix bwa fastp default-jre-headless perl git "
             "build-essential unzip curl zlib1g-dev libbz2-dev liblzma-dev "
             "libcurl4-openssl-dev libssl-dev libdbi-perl libarchive-zip-perl "
             "libwww-perl libjson-perl liblist-moreutils-perl cpanminus")


def _system_tools() -> None:
    """The command-line tools everything else shells out to: installed with Homebrew on
    macOS, checked for on Linux."""
    import platform
    if platform.system() == "Darwin":
        log.info("Installing tools via Homebrew: %s", " ".join(BREW_PKGS))
        run(["brew", "install", *BREW_PKGS])
        return
    missing = [t for t in LINUX_TOOLS if shutil.which(t) is None]
    if missing:
        raise SystemExit(
            f"Missing required tool(s): {', '.join(missing)}. On Debian/Ubuntu:\n"
            f"  sudo apt install {LINUX_APT}\nthen re-run `setup`.")
    log.info("System tools present (%s).", ", ".join(LINUX_TOOLS))


def setup() -> None:
    ensure_dirs()
    _system_tools()

    # Reference FASTA — used by both `bcftools norm` and `vep --fasta`.
    fa_gz = config.REFS_DIR / f"{config.BUILD}.primary_assembly.fa.gz"
    _download(config.REF_FASTA_URL, fa_gz, f"{config.BUILD} reference FASTA")
    if not config.REF_FASTA.exists():
        log.info("Decompressing reference FASTA …")
        with gzip.open(fa_gz, "rb") as fi, open(config.REF_FASTA, "wb") as fo:
            shutil.copyfileobj(fi, fo, length=1 << 24)
    if not config.REF_FASTA.with_suffix(".fa.fai").exists():
        run(["samtools", "faidx", str(config.REF_FASTA)])

    _install_alignment_refs()
    _install_gatk()
    _install_pharmcat()
    _install_expansionhunter()
    _install_haplogrep()
    _install_alphamissense()
    _install_spliceai()
    _install_gnomad_constraint()
    _install_manta()
    _install_arcashla()
    _install_mosdepth()
    _install_cyrius()
    _install_plink2()
    _install_liftover_refs()
    _install_vep()
    log.info("Setup complete. Next: `align <sample>` → `call <sample>` → `ingest …`")


def _install_alignment_refs() -> None:
    """Build the bwa index + GATK sequence dictionary for the reference FASTA.
    The bwa index is ~4.5GB and takes ~1h; both steps are idempotent."""
    ref = config.REF_FASTA
    if not Path(str(ref) + ".bwt").exists():
        log.info("Building bwa index for %s (~1h, ~4.5GB) …", ref.name)
        run(["bwa", "index", str(ref)])
    else:
        log.info("bwa index already present — skipping.")

    seqdict = ref.with_suffix(".dict")  # GATK expects <basename>.dict beside the FASTA
    if not seqdict.exists():
        log.info("Building sequence dictionary %s …", seqdict.name)
        run(["samtools", "dict", "-o", str(seqdict), str(ref)])
    else:
        log.info("Sequence dictionary already present — skipping.")


def _install_gatk() -> None:
    """Download + unzip the GATK4 release (pure Java; runs native via openjdk@17)."""
    if config.GATK_BIN.exists():
        log.info("GATK %s already installed — skipping.", config.GATK_VERSION)
        return
    config.GATK_DIR.parent.mkdir(parents=True, exist_ok=True)
    zip_path = config.GATK_DIR.parent / f"gatk-{config.GATK_VERSION}.zip"
    _download(config.GATK_URL, zip_path, f"GATK {config.GATK_VERSION}")
    log.info("Unzipping GATK → %s", config.GATK_DIR)
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(config.GATK_DIR.parent)  # extracts gatk-<ver>/
    config.GATK_BIN.chmod(0o755)
    zip_path.unlink(missing_ok=True)
    log.info("GATK installed at %s", config.GATK_BIN)


def _install_pharmcat() -> None:
    """Download the PharmCAT jar + Python preprocessor (for the `pgx` stage)."""
    if config.PHARMCAT_JAR.exists() and config.PHARMCAT_PIPELINE.exists():
        log.info("PharmCAT %s already installed — skipping.", config.PHARMCAT_VERSION)
        return
    config.PHARMCAT_DIR.mkdir(parents=True, exist_ok=True)
    _download(config.PHARMCAT_JAR_URL, config.PHARMCAT_JAR,
              f"PharmCAT {config.PHARMCAT_VERSION} jar")
    pre_tgz = config.PHARMCAT_DIR / f"pharmcat-preprocessor-{config.PHARMCAT_VERSION}.tar.gz"
    _download(config.PHARMCAT_PREPROCESSOR_URL, pre_tgz, "PharmCAT preprocessor")
    import tarfile
    with tarfile.open(pre_tgz) as t:
        t.extractall(config.PHARMCAT_DIR, filter="data")  # → preprocessor/
    config.PHARMCAT_PIPELINE.chmod(0o755)
    pre_tgz.unlink(missing_ok=True)
    log.info("PharmCAT installed at %s", config.PHARMCAT_DIR)


def _install_expansionhunter() -> None:
    """Download Illumina ExpansionHunter (repeat-expansion genotyping) for this platform."""
    import platform
    if list(config.EXPANSIONHUNTER_DIR.glob("*/bin/ExpansionHunter")):
        log.info("ExpansionHunter already installed — skipping.")
        return
    config.EXPANSIONHUNTER_DIR.mkdir(parents=True, exist_ok=True)
    asset = "macOS" if platform.system() == "Darwin" else "linux_x86_64"
    v = config.EXPANSIONHUNTER_VERSION
    url = (f"https://github.com/Illumina/ExpansionHunter/releases/download/v{v}/"
           f"ExpansionHunter-v{v}-{asset}.tar.gz")
    tgz = config.EXPANSIONHUNTER_DIR / f"eh-{v}.tar.gz"
    _download(url, tgz, f"ExpansionHunter {v} ({asset})")
    import tarfile
    with tarfile.open(tgz) as t:
        t.extractall(config.EXPANSIONHUNTER_DIR, filter="data")
    tgz.unlink(missing_ok=True)
    log.info("ExpansionHunter installed under %s", config.EXPANSIONHUNTER_DIR)


def _install_haplogrep() -> None:
    """Download Haplogrep 3 (for mtDNA haplogroups in the `ancestry` stage). The
    'linux' tarball is just a launcher + the portable haplogrep3.jar (cross-platform).
    Yleaf (Y haplogroup) is fetched on demand by uv, so nothing to install here."""
    if config.HAPLOGREP_JAR.exists():
        log.info("Haplogrep %s already installed — skipping.", config.HAPLOGREP_VERSION)
        return
    config.HAPLOGREP_DIR.mkdir(parents=True, exist_ok=True)
    tgz = config.HAPLOGREP_DIR / f"haplogrep3-{config.HAPLOGREP_VERSION}.tar.gz"
    _download(config.HAPLOGREP_URL, tgz, f"Haplogrep {config.HAPLOGREP_VERSION}")
    import tarfile
    with tarfile.open(tgz) as t:
        t.extractall(config.HAPLOGREP_DIR, filter="data")  # → haplogrep3.jar + launcher + data/
    tgz.unlink(missing_ok=True)
    log.info("Haplogrep installed at %s", config.HAPLOGREP_JAR)


def _install_alphamissense() -> None:
    """Download AlphaMissense (deep missense pathogenicity) and re-compress as bgzip +
    tabix so reports can point-query it. ~1GB download, ~5GB intermediate; optional —
    reports degrade gracefully without it."""
    if config.ALPHAMISSENSE_FILE.exists():
        log.info("AlphaMissense already installed — skipping.")
        return
    require_tools("bgzip", "tabix")
    raw = config.REFS_DIR / "AlphaMissense_hg38.tsv.gz"
    _download(config.ALPHAMISSENSE_URL, raw, "AlphaMissense GRCh38 (~1GB)")
    log.info("Re-compressing AlphaMissense as bgzip + tabix (slow) …")
    with open(config.ALPHAMISSENSE_FILE, "wb") as out:
        gz = subprocess.Popen(["gzip", "-dc", str(raw)], stdout=subprocess.PIPE)
        bg = subprocess.Popen(["bgzip", "-c"], stdin=gz.stdout, stdout=out)
        gz.stdout.close()
        if bg.wait() or gz.wait():
            raise SystemExit("AlphaMissense bgzip failed")
    run(["tabix", "-s", "1", "-b", "2", "-e", "2", str(config.ALPHAMISSENSE_FILE)])
    raw.unlink(missing_ok=True)
    log.info("AlphaMissense ready: %s", config.ALPHAMISSENSE_FILE)


def _install_spliceai() -> None:
    """Tabix-index the precomputed SpliceAI VCFs (deep splice-disruption prediction) so
    reports can point-query them. The VCFs themselves are NOT auto-downloaded — they live on
    Illumina BaseSpace behind a login, with no anonymous URL. Drop them in REFS_DIR and this
    indexes them; otherwise it logs how to obtain them and skips (reports degrade gracefully).
    """
    files = [config.SPLICEAI_SNV_FILE, config.SPLICEAI_INDEL_FILE]
    present = [f for f in files if f.exists()]
    if not present:
        log.info("SpliceAI VCFs not found — skipping (optional). %s",
                 config.SPLICEAI_OBTAIN)
        return
    require_tools("tabix")
    for vcf in present:
        if vcf.with_suffix(vcf.suffix + ".tbi").exists():
            log.info("SpliceAI %s already indexed — skipping.", vcf.name)
            continue
        log.info("Tabix-indexing SpliceAI %s …", vcf.name)
        run(["tabix", "-p", "vcf", str(vcf)])
    log.info("SpliceAI ready: %s", ", ".join(f.name for f in present))


def _hg38_chr_fasta() -> Path:
    """A chr-prefixed copy of REF_FASTA, cached — LiftoverVcf's REFERENCE_SEQUENCE needs to
    match the UCSC chain's target-side contig names ('chr1'), but our own reference is
    Ensembl-named ('1'). Rewriting headers is a text substitution over an already-downloaded
    file, far cheaper than a second ~3GB hg38 download. Restricted to the canonical contigs
    (1-22, X, Y, MT)."""
    out = config.HG38_CHR_FASTA
    if out.exists() and out.with_suffix(".fa.fai").exists():
        return out
    log.info("Building chr-prefixed reference copy for liftover (%s) …", out.name)
    canonical = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
    with open(config.REF_FASTA) as fi, open(out, "w") as fo:
        keep = False
        for line in fi:
            if line.startswith(">"):
                name = line[1:].split()[0]
                keep = name in canonical
                if keep:
                    fo.write(f">{'chrM' if name == 'MT' else 'chr' + name}\n")
            elif keep:
                fo.write(line)
    run(["samtools", "faidx", str(out)])
    seqdict = out.with_suffix(".dict")
    seqdict.unlink(missing_ok=True)
    run(["samtools", "dict", "-o", str(seqdict), str(out)])
    return out


def _install_liftover_refs() -> None:
    """The hg19->GRCh38 chain and a chr-named reference copy, for lifting an older-build
    callset with GATK LiftoverVcf before `ingest` (README: "The VCF must be GRCh38").
    Best-effort: a failure logs and skips, and re-running `setup` retries."""
    try:
        if not config.HG19TOHG38_CHAIN.exists():
            config.LIFTOVER_DIR.mkdir(parents=True, exist_ok=True)
            _download(config.HG19TOHG38_CHAIN_URL, config.HG19TOHG38_CHAIN,
                      "hg19->hg38 liftover chain")
        _hg38_chr_fasta()
    except Exception as e:  # noqa: BLE001 — optional reference data
        log.warning("Liftover references not prepared (%s) — only needed to lift an "
                    "hg19/GRCh37 VCF; re-run `setup` later to retry.", e)


def mirror_1000g(chroms: list[str] | None = None) -> None:
    """One-time local mirror of the 1000G 3202-sample phased autosome VCFs (~34.6GB total),
    downloaded to `config.KG_VCF_LOCAL_DIR`. Once present, `config.kg_vcf_source()` prefers
    the local copy over the remote URL everywhere it's used (prs.py's genome-wide PRS
    percentile + calibration).

    Why: those callers do `bcftools query -R <hundreds-to-thousands of scattered
    positions>` against the remote VCF over plain HTTP. Measured directly: a handful of
    positions returns in under a second, but ~1000 positions scattered across a whole
    chromosome didn't finish in 5 minutes — despite the same server sustaining ~17MB/s on a
    plain sequential GET. It's not the server being slow; it's the per-non-contiguous-region
    HTTP range-request overhead compounding at scale. A local file turns every one of those
    into a fast local tabix seek instead. Idempotent — skips chromosomes already present."""
    config.KG_VCF_LOCAL_DIR.mkdir(parents=True, exist_ok=True)
    for chrom in chroms or [str(c) for c in range(1, 23)]:
        dest = config.KG_VCF_LOCAL_DIR / f"chr{chrom}.vcf.gz"
        if dest.exists() and (dest.with_suffix(".gz.tbi")).exists():
            log.info("1000G chr%s already mirrored — skipping.", chrom)
            continue
        url = config.KG_VCF_URL_TMPL.format(chrom=chrom)
        _download(url, dest, f"1000G chr{chrom} (phased, 3202 samples)")
        _download(url + ".tbi", Path(f"{dest}.tbi"), f"1000G chr{chrom} index")
    log.info("1000G local mirror ready: %s", config.KG_VCF_LOCAL_DIR)


def _reduce_constraint(lines) -> list[tuple[str, str, str]]:
    """Reduce gnomAD's per-transcript constraint metrics to one (gene, pLI, LOEUF) per gene.

    The metrics file has many transcripts per gene; we keep the MANE-select transcript when
    flagged, else the canonical one, else the first seen — a stable, documented choice.
    Pure (operates on an iterable of TSV lines) so it can be unit-tested without the download."""
    it = iter(lines)
    try:
        header = {c: i for i, c in enumerate(next(it).rstrip("\n").split("\t"))}
    except StopIteration:
        return []
    g, mane, canon = header["gene"], header["mane_select"], header["canonical"]
    pli, loeuf = header["lof.pLI"], header["lof.oe_ci.upper"]
    n = max(g, mane, canon, pli, loeuf)
    best: dict[str, tuple[int, str, str]] = {}     # gene → (priority, pLI, LOEUF)
    for line in it:
        f = line.rstrip("\n").split("\t")
        if len(f) <= n or not f[g]:
            continue
        prio = 2 if f[mane] == "true" else (1 if f[canon] == "true" else 0)
        if f[g] not in best or prio > best[f[g]][0]:
            best[f[g]] = (prio, f[pli], f[loeuf])
    return [(gene, p, l) for gene, (_pr, p, l) in sorted(best.items())]


def _install_gnomad_constraint() -> None:
    """Download gnomAD v4.1 per-transcript constraint metrics (~95MB) and reduce them to a
    small gene-keyed table (gene → pLI, LOEUF) the ACMG engine reads to grade PVS1. Optional —
    the engine degrades gracefully (PVS1 unchanged) when the table is absent."""
    if config.GNOMAD_CONSTRAINT_FILE.exists():
        log.info("gnomAD constraint already installed — skipping.")
        return
    raw = config.REFS_DIR / "gnomad.v4.1.constraint_metrics.tsv"
    _download(config.GNOMAD_CONSTRAINT_URL, raw,
              f"gnomAD {config.GNOMAD_CONSTRAINT_VERSION} constraint metrics (~95MB)")
    log.info("Reducing constraint metrics to one row per gene …")
    with open(raw) as fh:
        rows = _reduce_constraint(fh)
    with open(config.GNOMAD_CONSTRAINT_FILE, "w") as out:
        out.write("gene\tpLI\tLOEUF\n")
        out.writelines(f"{gene}\t{p}\t{l}\n" for gene, p, l in rows)
    raw.unlink(missing_ok=True)
    log.info("gnomAD constraint ready: %s (%d genes)", config.GNOMAD_CONSTRAINT_FILE, len(rows))


def _install_mosdepth() -> None:
    """Install mosdepth — a single static binary. Linux
    release binary into MOSDEPTH_DIR; on macOS fall back to Homebrew (the `callable` stage
    usually runs on the Linux box anyway). Located by glob under MOSDEPTH_DIR or PATH."""
    import platform
    import shutil
    import stat
    if shutil.which("mosdepth") or list(config.MOSDEPTH_DIR.glob("mosdepth")):
        log.info("mosdepth already available — skipping.")
        return
    if platform.system() == "Darwin":
        log.info("On macOS install mosdepth via Homebrew: `brew install mosdepth`.")
        return
    config.MOSDEPTH_DIR.mkdir(parents=True, exist_ok=True)
    v = config.MOSDEPTH_VERSION
    dest = config.MOSDEPTH_DIR / "mosdepth"
    _download(f"https://github.com/brentp/mosdepth/releases/download/v{v}/mosdepth",
              dest, f"mosdepth {v} (linux binary)")
    dest.chmod(dest.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    log.info("mosdepth installed: %s", dest)


def _install_plink2() -> None:
    """Download the plink2 binary (genome-wide PRS scoring). Single static binary in a zip;
    mac arm64 or linux. Located by glob under PLINK2_DIR or PATH. Idempotent."""
    import platform
    import shutil
    import zipfile
    if shutil.which("plink2") or list(config.PLINK2_DIR.glob("plink2")):
        log.info("plink2 already available — skipping.")
        return
    config.PLINK2_DIR.mkdir(parents=True, exist_ok=True)
    url = config.PLINK2_URL_MAC if platform.system() == "Darwin" else config.PLINK2_URL_LINUX
    z = config.PLINK2_DIR / "plink2.zip"
    _download(url, z, "plink2")
    with zipfile.ZipFile(z) as zf:
        zf.extractall(config.PLINK2_DIR)
    binp = config.PLINK2_DIR / "plink2"
    binp.chmod(binp.stat().st_mode | 0o111)
    z.unlink(missing_ok=True)
    log.info("plink2 installed: %s", binp)


def _install_cyrius() -> None:
    """Clone Cyrius (CYP2D6 star-allele caller, feeds the PGx report) and reconcile its
    chr-prefixed GRCh38 contigs to our Ensembl naming — in the data files (contig column)
    and the one hardcoded ("chrNN", ...) region tuple in the caller source. Run-time deps
    come from `uv run --with` (numpy/scipy/pysam/statsmodels). Idempotent."""
    import re
    if not (config.CYRIUS_DIR / "star_caller.py").exists():
        require_tools("git")
        config.CYRIUS_DIR.parent.mkdir(parents=True, exist_ok=True)
        run(["git", "clone", "--depth", "1", config.CYRIUS_REPO, str(config.CYRIUS_DIR)])
    data = config.CYRIUS_DIR / "data"
    for name in ("CYP2D6_region_38.bed", "CYP2D6_SNP_38.txt", "CYP2D6_target_variant_38.txt",
                 "CYP2D6_target_variant_homology_region_38.txt", "CYP2D6_haplotype_38.txt"):
        f = data / name
        if f.exists():                       # strip leading "chr" from the contig column
            f.write_text(re.sub(r"(?m)^chr([0-9XYM])", r"\1", f.read_text()))
    cv = config.CYRIUS_DIR / "caller" / "call_variants.py"
    if cv.exists():                          # ("chr22", ...) → ("22", ...)
        cv.write_text(re.sub(r'\("chr([0-9XYM]+)"', r'("\1"', cv.read_text()))
    log.info("Cyrius ready: %s (GRCh38 contigs Ensembl-ified)", config.CYRIUS_DIR)


def _install_manta() -> None:
    """Pre-pull the Manta SV-caller Docker image. Manta's workflow
    scripts need Python 2, so we run it containerized rather than installing a binary —
    best-effort, since the `sv` stage typically runs on the Linux box (Docker + CRAMs),
    not the Mac. Skips quietly if Docker isn't present or the image is already local."""
    import shutil
    if shutil.which("docker") is None:
        log.info("Docker not found — skipping Manta image pull (the `sv` stage needs "
                 "Docker, usually on the Linux box).")
        return
    have = subprocess.run(["docker", "image", "inspect", config.MANTA_IMAGE],
                          capture_output=True)
    if have.returncode == 0:
        log.info("Manta image already present (%s) — skipping.", config.MANTA_IMAGE)
        return
    log.info("Pulling Manta image %s …", config.MANTA_IMAGE)
    pull = subprocess.run(["docker", "pull", config.MANTA_IMAGE])
    if pull.returncode:
        log.warning("Could not pull %s — run `sv` on a box with Docker access, or set "
                    "MANTA_IMAGE.", config.MANTA_IMAGE)
    else:
        log.info("Manta image ready: %s", config.MANTA_IMAGE)


def _install_arcashla() -> None:
    """Pre-pull the arcasHLA Docker image (full classical HLA typing).
    Same rationale as Manta — arcasHLA ships via bioconda, not a portable binary, so we
    run it containerized; best-effort since `hla-type` typically runs on the Linux box
    (Docker + CRAMs), not the Mac. Skips quietly if Docker isn't present."""
    import shutil
    if shutil.which("docker") is None:
        log.info("Docker not found — skipping arcasHLA image pull (the `hla-type` stage "
                 "needs Docker, usually on the Linux box).")
        return
    have = subprocess.run(["docker", "image", "inspect", config.ARCASHLA_IMAGE],
                          capture_output=True)
    if have.returncode == 0:
        log.info("arcasHLA image already present (%s) — skipping.", config.ARCASHLA_IMAGE)
        return
    log.info("Pulling arcasHLA image %s …", config.ARCASHLA_IMAGE)
    pull = subprocess.run(["docker", "pull", config.ARCASHLA_IMAGE])
    if pull.returncode:
        log.warning("Could not pull %s — run `hla-type` on a box with Docker access, "
                    "or set ARCASHLA_IMAGE.", config.ARCASHLA_IMAGE)
    else:
        log.info("arcasHLA image ready: %s", config.ARCASHLA_IMAGE)


def _install_vep() -> None:
    """Clone ensembl-vep and run its installer for the API, htslib, and cache."""
    require_tools("git", "cc", "make")
    cache_marker = config.VEP_CACHE_DIR / "homo_sapiens" / \
        f"{config.VEP_CACHE_VERSION}_{config.BUILD}"
    if config.VEP_BIN.exists() and cache_marker.exists():
        log.info("VEP already installed with %s cache r%s — skipping.",
                 config.BUILD, config.VEP_CACHE_VERSION)
        return

    # If VEP itself is already built (e.g. from a previous build), only the assembly
    # cache is missing — fetch just that ('c') and skip the API/htslib rebuild. The
    # authoritative test is whether Bio::DB::HTS actually loads (the `vep` script and
    # even a partial htslib dir can exist without a usable build).
    probe = subprocess.run(
        ["perl", "-e", "use Bio::DB::HTS::Tabix;"],
        env=dict(os.environ, PERL5LIB=str(config.VEP_DIR)),
        capture_output=True) if config.VEP_BIN.exists() else None
    vep_built = bool(probe and probe.returncode == 0)
    # AUTO letters: a=API, l=Bio::DB::HTS (builds htslib), c=indexed cache, f=species
    # FASTA. We annotate against our own refs/ FASTA, so never fetch VEP's (no 'f').
    auto = "c" if vep_built else "alc"

    if not config.VEP_DIR.exists():
        config.VEP_DIR.parent.mkdir(parents=True, exist_ok=True)
        log.info("Cloning ensembl-vep → %s", config.VEP_DIR)
        run(["git", "clone", "--depth", "1", "--branch",
             f"release/{config.VEP_CACHE_VERSION}",
             config.VEP_REPO_URL, str(config.VEP_DIR)])

    # Clear any partial API download from an interrupted run so unzip starts clean.
    stale_tmp = config.VEP_DIR / "Bio" / "tmp"
    if stale_tmp.exists():
        log.info("Removing stale partial API download: %s", stale_tmp)
        shutil.rmtree(stale_tmp, ignore_errors=True)

    log.info("Running VEP INSTALL.pl (--AUTO %s) for %s cache r%s%s …",
             auto, config.BUILD, config.VEP_CACHE_VERSION,
             " — cache only" if vep_built else " (API + htslib + ~25GB cache; slow)")
    run([
        "perl", "INSTALL.pl",
        "--AUTO", auto,
        "--SPECIES", "homo_sapiens",
        "--ASSEMBLY", config.BUILD,
        "--CACHEDIR", str(config.VEP_CACHE_DIR),
        "--CACHE_VERSION", config.VEP_CACHE_VERSION,
        "--NO_UPDATE", "--NO_TEST",
    ], cwd=str(config.VEP_DIR), env=_build_env())

    if not vep_built:
        _fix_vep_dylibs()
        # Sanity check: Bio::DB::HTS (needed for --custom ClinVar + --fasta) must load.
        log.info("Verifying VEP / Bio::DB::HTS …")
        env = dict(os.environ, PERL5LIB=str(config.VEP_DIR))
        run(["perl", "-e", "use Bio::DB::HTS::Tabix; use Bio::DB::HTS::Faidx;"], env=env)
    # Bio::DB::HTS loading is not the same as VEP running. On a fresh Linux box the build
    # succeeded and this stage reported "ready" while `vep` itself could not start, for want
    # of one Perl module (List::MoreUtils) — found only when the first annotation failed.
    # Starting VEP is the only honest check, and it costs a second.
    probe = subprocess.run([str(config.VEP_BIN), "--help"], capture_output=True, text=True,
                           env=dict(os.environ, PERL5LIB=str(config.VEP_DIR)))
    if probe.returncode != 0:
        missing = re.search(r"Can't locate (\S+)\.pm", probe.stderr or "")
        hint = (f" Perl module {missing.group(1).replace('/', '::')} is missing — install it "
                f"(Debian/Ubuntu: see the package list in the README) and re-run `setup`."
                if missing else "")
        raise SystemExit(f"VEP is installed but does not start.{hint}\n"
                         + (probe.stderr or "").strip()[-600:])
    log.info("VEP %s cache r%s ready at %s", config.BUILD,
             config.VEP_CACHE_VERSION, config.VEP_BIN)


def _fix_vep_dylibs() -> None:
    """The from-source htslib builds libhts with install_name /usr/local/lib/...,
    which doesn't exist on this Mac, so the Bio::DB::HTS XS bundles fail to load.
    Repoint them at @rpath (the bundles already carry an rpath to VEP_DIR/htslib).
    macOS-only — on Linux the ELF .so links resolve normally (no otool/install_name_tool)."""
    import platform
    if platform.system() != "Darwin":
        return
    htslib_dir = config.VEP_DIR / "htslib"
    for bundle in config.VEP_DIR.rglob("*.bundle"):
        out = subprocess.run(["otool", "-L", str(bundle)],
                             capture_output=True, text=True).stdout
        if "/usr/local/lib/libhts" not in out:
            continue
        log.info("Repointing libhts ref in %s", bundle.name)
        subprocess.run(["install_name_tool", "-change",
                        "/usr/local/lib/libhts.2.dylib", "@rpath/libhts.2.dylib",
                        str(bundle)], check=True)
        rpaths = subprocess.run(["otool", "-l", str(bundle)],
                                capture_output=True, text=True).stdout
        if str(htslib_dir) not in rpaths:
            subprocess.run(["install_name_tool", "-add_rpath",
                            str(htslib_dir), str(bundle)], check=False)
