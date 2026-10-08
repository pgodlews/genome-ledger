"""Central configuration: filesystem layout, external endpoints, thresholds.

Everything is rooted at GENOMES_ROOT (default ~/genomes, override with the
GENOMES_ROOT env var). Code lives in the git repo; data lives under this root and
is never committed.
"""

from __future__ import annotations

import os
from pathlib import Path

# --- Working root (data, never in git) ---------------------------------------
GENOMES_ROOT = Path(os.environ.get("GENOMES_ROOT", Path.home() / "genomes")).expanduser()

REFS_DIR = GENOMES_ROOT / "refs"
INCOMING_DIR = GENOMES_ROOT / "incoming"   # raw FASTQ reads from the sequencer
ALIGNED_DIR = GENOMES_ROOT / "aligned"     # per-sample marked-dup CRAMs
CALLED_DIR = GENOMES_ROOT / "called"       # per-sample VCFs from HaplotypeCaller
RAW_DIR = GENOMES_ROOT / "raw"
NORMALIZED_DIR = GENOMES_ROOT / "normalized"
SNAPSHOTS_DIR = GENOMES_ROOT / "db" / "snapshots"
VEP_CACHE_DIR = GENOMES_ROOT / "vep_cache"
ANNOTATED_DIR = GENOMES_ROOT / "annotated"
DUCKDB_DIR = GENOMES_ROOT / "duckdb"
REPORTS_DIR = GENOMES_ROOT / "reports"
MANIFEST_DIR = GENOMES_ROOT / "manifest"
LOGS_DIR = GENOMES_ROOT / "logs"

# The human-editable family graph (source of truth for relationships): one row per person —
# id, display_name, sex, father, mother, partner (parent/partner cells reference other ids).
# `pedigree.ped` is *generated* from this (for Mendelian/trio tooling — pipeline/pedigree.py).
# Sample ids are perspective-neutral given names; "father"/"mother" are no longer ids but
# relationships computed relative to a chosen person (the ego). See pipeline/pedigree.py.
FAMILY_FILE = GENOMES_ROOT / "family" / "family.tsv"
PEDIGREE_FILE = GENOMES_ROOT / "pedigree.ped"
# S1: per-(sample_id, category) opt-in for incidental findings (adult-onset/untreatable
# predictions with no clinical intervention) — one row per grant, default is NOT consented.
# See pipeline/consent.py.
CONSENT_FILE = GENOMES_ROOT / "family" / "consent.tsv"
SAMPLES_TSV = MANIFEST_DIR / "samples.tsv"
RUNS_TSV = MANIFEST_DIR / "runs.tsv"
DUCKDB_FILE = DUCKDB_DIR / "genomes.duckdb"
# Maps a sample id to its sequencer FASTQ files: sample, vendor_id, relation, sex. The
# files under INCOMING_DIR start with the vendor_id — `<vendor_id>_1.fq.gz` / `_2`,
# `_R1.fastq.gz` / `_R2`, or one `_S1_L001_R1_001.fastq.gz` pair per lane (call._fastq_files).
INCOMING_MANIFEST = INCOMING_DIR / "manifest.tsv"

ALL_DIRS = [
    REFS_DIR, INCOMING_DIR, ALIGNED_DIR, CALLED_DIR, RAW_DIR, NORMALIZED_DIR,
    SNAPSHOTS_DIR, VEP_CACHE_DIR, ANNOTATED_DIR, DUCKDB_DIR, REPORTS_DIR,
    MANIFEST_DIR, LOGS_DIR, GENOMES_ROOT / "pgx", GENOMES_ROOT / "ancestry",
    GENOMES_ROOT / "traits", GENOMES_ROOT / "hla", GENOMES_ROOT / "prs",
    GENOMES_ROOT / "repeats", GENOMES_ROOT / "sv", GENOMES_ROOT / "callable",
    GENOMES_ROOT / "phase", GENOMES_ROOT / "forcecalled", GENOMES_ROOT / "mito",
    GENOMES_ROOT / "interpret", GENOMES_ROOT / "hla_typing", GENOMES_ROOT / "denovo",
    GENOMES_ROOT / "segregation",
]


def set_genomes_root(root: Path | str) -> None:
    """Dynamically repoint all GENOMES_ROOT-derived paths to a custom root (e.g. for demo or testing)."""
    global GENOMES_ROOT, REFS_DIR, INCOMING_DIR, ALIGNED_DIR, CALLED_DIR, RAW_DIR
    global NORMALIZED_DIR, SNAPSHOTS_DIR, VEP_CACHE_DIR, ANNOTATED_DIR, DUCKDB_DIR
    global REPORTS_DIR, MANIFEST_DIR, LOGS_DIR, FAMILY_FILE, PEDIGREE_FILE
    global CONSENT_FILE, SAMPLES_TSV, RUNS_TSV, DUCKDB_FILE, INCOMING_MANIFEST
    global ALL_DIRS, REPORT_DIRS
    global MITO_DIR, ANCESTRY_DIR, REPEATS_DIR, SV_DIR, CALLABLE_DIR, PGX_DIR
    global TRAITS_DIR, HLA_DIR, HLA_TYPE_DIR, PRS_DIR, KG_VCF_LOCAL_DIR
    global DENOVO_DIR, SEGREGATION_DIR
    global PHASE_DIR, SHARE_DIR, FORCECALL_DIR, FC_SITES_DIR, INTERPRET_DIR, GATK_DIR, GATK_BIN, GATK_JAR
    global VEP_DIR, HAPLOGREP_DIR, EXPANSIONHUNTER_DIR, MOSDEPTH_DIR, PHARMCAT_DIR
    global CYRIUS_DIR, PGS_DIR, PLINK2_DIR, LIFTOVER_DIR
    global REF_FASTA, HG38_CHR_FASTA, HG19TOHG38_CHAIN
    global ALPHAMISSENSE_FILE, GNOMAD_CONSTRAINT_FILE, SPLICEAI_SNV_FILE, SPLICEAI_INDEL_FILE
    global ENFORMER_SCORES_FILE, ENFORMER_CALIBRATION_FILE

    GENOMES_ROOT = Path(root).expanduser().resolve()
    REFS_DIR = GENOMES_ROOT / "refs"
    INCOMING_DIR = GENOMES_ROOT / "incoming"
    ALIGNED_DIR = GENOMES_ROOT / "aligned"
    CALLED_DIR = GENOMES_ROOT / "called"
    RAW_DIR = GENOMES_ROOT / "raw"
    NORMALIZED_DIR = GENOMES_ROOT / "normalized"
    SNAPSHOTS_DIR = GENOMES_ROOT / "db" / "snapshots"
    VEP_CACHE_DIR = GENOMES_ROOT / "vep_cache"
    ANNOTATED_DIR = GENOMES_ROOT / "annotated"
    DUCKDB_DIR = GENOMES_ROOT / "duckdb"
    REPORTS_DIR = GENOMES_ROOT / "reports"
    MANIFEST_DIR = GENOMES_ROOT / "manifest"
    LOGS_DIR = GENOMES_ROOT / "logs"

    FAMILY_FILE = GENOMES_ROOT / "family" / "family.tsv"
    PEDIGREE_FILE = GENOMES_ROOT / "pedigree.ped"
    CONSENT_FILE = GENOMES_ROOT / "family" / "consent.tsv"
    SAMPLES_TSV = MANIFEST_DIR / "samples.tsv"
    RUNS_TSV = MANIFEST_DIR / "runs.tsv"
    DUCKDB_FILE = DUCKDB_DIR / "genomes.duckdb"
    INCOMING_MANIFEST = INCOMING_DIR / "manifest.tsv"

    MITO_DIR = GENOMES_ROOT / "mito"
    ANCESTRY_DIR = GENOMES_ROOT / "ancestry"
    REPEATS_DIR = GENOMES_ROOT / "repeats"
    SV_DIR = GENOMES_ROOT / "sv"
    CALLABLE_DIR = GENOMES_ROOT / "callable"
    PGX_DIR = GENOMES_ROOT / "pgx"
    TRAITS_DIR = GENOMES_ROOT / "traits"
    HLA_DIR = GENOMES_ROOT / "hla"
    HLA_TYPE_DIR = GENOMES_ROOT / "hla_typing"
    PRS_DIR = GENOMES_ROOT / "prs"
    KG_VCF_LOCAL_DIR = GENOMES_ROOT / "refs" / "1000G_phased"
    DENOVO_DIR = GENOMES_ROOT / "denovo"
    SEGREGATION_DIR = GENOMES_ROOT / "segregation"
    PHASE_DIR = GENOMES_ROOT / "phase"
    SHARE_DIR = GENOMES_ROOT / "share"
    FORCECALL_DIR = GENOMES_ROOT / "forcecalled"
    FC_SITES_DIR = FORCECALL_DIR / "sites"
    INTERPRET_DIR = GENOMES_ROOT / "interpret"
    ENFORMER_SCORES_FILE = INTERPRET_DIR / f"{INTERPRET_MODEL}_scores.tsv"
    ENFORMER_CALIBRATION_FILE = INTERPRET_DIR / f"{INTERPRET_MODEL}_calibration.json"

    GATK_DIR = GENOMES_ROOT / "tools" / f"gatk-{GATK_VERSION}"
    GATK_BIN = GATK_DIR / "gatk"
    GATK_JAR = GATK_DIR / f"gatk-package-{GATK_VERSION}-local.jar"
    VEP_DIR = GENOMES_ROOT / "tools" / "ensembl-vep"
    HAPLOGREP_DIR = GENOMES_ROOT / "tools" / "haplogrep"
    EXPANSIONHUNTER_DIR = GENOMES_ROOT / "tools" / "expansionhunter"
    MOSDEPTH_DIR = GENOMES_ROOT / "tools" / "mosdepth"
    PHARMCAT_DIR = GENOMES_ROOT / "tools" / "pharmcat"
    CYRIUS_DIR = GENOMES_ROOT / "tools" / "Cyrius"
    PGS_DIR = GENOMES_ROOT / "tools" / "pgs"
    PLINK2_DIR = GENOMES_ROOT / "tools" / "plink2"
    LIFTOVER_DIR = GENOMES_ROOT / "refs" / "liftover"

    REF_FASTA = REFS_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    HG38_CHR_FASTA = REFS_DIR / "GRCh38.chr_named.fa"
    HG19TOHG38_CHAIN = LIFTOVER_DIR / "hg19ToHg38.over.chain.gz"
    ALPHAMISSENSE_FILE = REFS_DIR / "AlphaMissense_hg38.tsv.bgz"
    GNOMAD_CONSTRAINT_FILE = REFS_DIR / "gnomad_constraint_by_gene.tsv"
    SPLICEAI_SNV_FILE = REFS_DIR / "spliceai_scores.masked.snv.hg38.vcf.gz"
    SPLICEAI_INDEL_FILE = REFS_DIR / "spliceai_scores.masked.indel.hg38.vcf.gz"

    ALL_DIRS = [
        REFS_DIR, INCOMING_DIR, ALIGNED_DIR, CALLED_DIR, RAW_DIR, NORMALIZED_DIR,
        SNAPSHOTS_DIR, VEP_CACHE_DIR, ANNOTATED_DIR, DUCKDB_DIR, REPORTS_DIR,
        MANIFEST_DIR, LOGS_DIR, PGX_DIR, ANCESTRY_DIR,
        TRAITS_DIR, HLA_DIR, PRS_DIR,
        REPEATS_DIR, SV_DIR, CALLABLE_DIR,
        PHASE_DIR, FORCECALL_DIR, MITO_DIR,
        INTERPRET_DIR, HLA_TYPE_DIR, DENOVO_DIR,
        SEGREGATION_DIR,
    ]

    REPORT_DIRS = (REPORTS_DIR,)

# --- Reference build ----------------------------------------------------------
# The pipeline now calls variants itself from FASTQ on GRCh38 (see pipeline/call.py),
# so GRCh38 is the source of truth. Reverting to GRCh37 is just this block.
BUILD = "GRCh38"
# Ensembl GRCh38 primary assembly FASTA. Used as the alignment + calling reference
# (bwa/GATK) and by `bcftools norm -f` / `vep --fasta`. Contigs are named
# 1..22,X,Y,MT,KI270*/GL000* — Ensembl style, matching ClinVar/VEP GRCh38, so our
# self-called VCFs need no contig renaming downstream.
REF_FASTA = REFS_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
REF_FASTA_URL = (
    "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/"
    "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
)

# Contig rename map for `bcftools annotate --rename-chrs` in normalize. Our own
# GRCh38 calls are already Ensembl-named (1.. style), so this is a no-op for them;
# it is retained so a future chr-prefixed VCF (e.g. a vendor file) still normalizes.
CONTIG_RENAME = {**{f"chr{i}": str(i) for i in range(1, 23)},
                 "chrX": "X", "chrY": "Y", "chrM": "MT"}

# --- Alignment & variant calling (FASTQ → VCF front-end, pipeline/call.py) -----
# Sequencer platform tag (PL) for the BAM @RG line. Unset, `align` reads it off the first
# read name in the FASTQ (call.read_platform): MGI/DNBSEQ and Illumina name reads in formats
# that cannot be confused. Set READ_PLATFORM to override detection — e.g. ELEMENT for AVITI
# reads, which use Illumina's read-name format.
READ_PLATFORM = os.environ.get("READ_PLATFORM", "").strip() or None
# bwa-mem alignment threads (the CPU bottleneck); leave a few cores for the
# concurrent fastp/sort/markdup processes in the same pipe. Override with ALIGN_THREADS.
ALIGN_THREADS = int(os.environ.get("ALIGN_THREADS", max(4, (os.cpu_count() or 8) - 4)))
SORT_THREADS = int(os.environ.get("SORT_THREADS", 4))
SORT_MEM = os.environ.get("SORT_MEM", "3G")  # per sort thread

# GATK4 (HaplotypeCaller) — pure-Java, runs native on Apple Silicon (no Docker).
# Downloaded as a release zip by `setup`; needs a JRE (openjdk@17 via Homebrew).
GATK_VERSION = "4.6.1.0"
GATK_DIR = GENOMES_ROOT / "tools" / f"gatk-{GATK_VERSION}"
GATK_BIN = GATK_DIR / "gatk"  # python launcher; we invoke the jar directly instead
# The `gatk` launcher is a Python script and this Mac has no `python` (only python3),
# so call.py runs `java -jar <local.jar> <Tool>` — fully supported for HaplotypeCaller.
GATK_JAR = GATK_DIR / f"gatk-package-{GATK_VERSION}-local.jar"
GATK_URL = (f"https://github.com/broadinstitute/gatk/releases/download/"
            f"{GATK_VERSION}/gatk-{GATK_VERSION}.zip")
JDK_FORMULA = "openjdk@17"  # keg-only; call.py puts its bin on PATH for gatk

# HaplotypeCaller has no internal multithreading in GATK4, so we scatter by contig
# and run many single-threaded callers in parallel, then bcftools-concat the parts.
# Restricted to the assembled chromosomes; unplaced scaffolds carry no ClinVar sites.
# Diploid contigs (X/Y left diploid — PAR/sex-ploidy modelling is a separate concern).
CALL_CONTIGS = [str(i) for i in range(1, 23)] + ["X", "Y"]
# Mitochondrion is called separately as haploid: it is a single non-recombining,
# high-depth contig, so diploid calling is wrong. ploidy-1 HaplotypeCaller is the
# simple, adequate approach for germline + haplogroup assignment. Heteroplasmy (a mix of
# mtDNA molecules carrying different alleles) is the `mito` stage's job — see below.
# Ensembl GRCh38 names it "MT" (it is the rCRS, so standard m. positions map 1:1).
MT_CONTIG = os.environ.get("MT_CONTIG", "MT")
MT_PLOIDY = 1

# --- Mitochondrial heteroplasmy (Mutect2 mito-mode, pipeline/mito.py) -- C4 ---
# Haploid HaplotypeCaller gives a present/absent germline mtDNA call but no *heteroplasmy*
# fraction — what proportion of mtDNA molecules carry a variant. Mito disease is dosage-
# dependent (a pathogenic allele matters above a tissue-specific threshold), so the fraction
# is the clinically meaningful number. GATK Mutect2 in --mitochondria-mode estimates the
# alt-allele fraction (AF) per site off the CRAM; FilterMutectCalls --mitochondria-mode flags
# low-confidence calls. Runs off the CRAM only (no ClinVar), so output lives outside snapshots.
MITO_DIR = GENOMES_ROOT / "mito"
MITO_VARIANTS_FILE = Path(__file__).resolve().parent.parent / "panels" / "mito_variants.tsv"
# AF at/above this is effectively homoplasmic (whole-population); below it is heteroplasmic.
MITO_HOMOPLASMIC_AF = 0.95
# AF below this is treated as low-level/noise — reported but not headlined (Mutect2's mito
# sensitivity reaches ~1%, but very low fractions are easily artefacts at high depth).
MITO_MIN_HET_AF = 0.03
HC_WORKERS = int(os.environ.get("HC_WORKERS", min(len(CALL_CONTIGS) + 1,
                                                  max(1, (os.cpu_count() or 8) - 2))))
HC_JAVA_MEM = os.environ.get("HC_JAVA_MEM", "4g")  # -Xmx per HaplotypeCaller JVM

# --- 1000G reference-panel builds (pipeline/prs.py) --------------------------
# Each chromosome's panel is independent — its own input VCF, its own output pgen, no shared
# state — so they parallelise the same way HC_WORKERS scatters HaplotypeCaller. Serially,
# building the panels for the 22 genome-wide panel scores measured ~38 h of one core.
# workers x bcftools threads should stay near the core count; past that they contend.
# Scoring-file downloads: bounded so a stalled connection fails instead of hanging a scan.
PGS_DOWNLOAD_TIMEOUT = int(os.environ.get("PGS_DOWNLOAD_TIMEOUT", "300"))
PRS_REF_WORKERS = int(os.environ.get("PRS_REF_WORKERS",
                                     max(1, min(22, (os.cpu_count() or 8) // 2))))
PRS_BCFTOOLS_THREADS = int(os.environ.get("PRS_BCFTOOLS_THREADS", 2))
# Scatter by sub-chromosomal windows rather than whole contigs, so the slow,
# repeat-rich chromosomes (chr16!) parallelize instead of leaving one core to grind
# for hours at the tail. ~30Mb windows → ~110 shards across the genome.
CALL_WINDOW_BP = int(os.environ.get("CALL_WINDOW_BP", 30_000_000))

# VEP is installed from source (no Homebrew/conda formula on this Mac). INSTALL.pl
# downloads the indexed offline cache (gnomAD AF + consequence prediction) into
# VEP_CACHE_DIR and builds Bio::DB::HTS inside VEP_DIR.
VEP_CACHE_VERSION = "112"
VEP_DIR = GENOMES_ROOT / "tools" / "ensembl-vep"
VEP_BIN = VEP_DIR / "vep"
VEP_REPO_URL = "https://github.com/Ensembl/ensembl-vep.git"
# Parallel VEP workers. WGS annotation is CPU-bound; leave a couple of cores free.
# Override with the VEP_FORKS env var.
VEP_FORKS = int(os.environ.get("VEP_FORKS", max(1, min((os.cpu_count() or 4) - 2, 12))))

# --- ClinVar (the version-pinned, fast-moving signal) ------------------------
# The "current" weekly file; we copy it into a dated snapshot on download.
CLINVAR_VCF_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
# Directory of dated weekly releases, used by `check-releases` to detect new versions.
CLINVAR_WEEKLY_INDEX = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/"

# gnomAD AF arrives via the VEP cache, so its version tracks VEP_CACHE_VERSION.
# The GRCh38 r112 cache carries gnomAD genomes v3.1.2 + exomes v2.1.1 (liftover) —
# a much larger cohort than the GRCh37 v2-only cache. Confirm on cache download.
GNOMAD_VERSION = "genomes v3.1.2 + exomes v2.1.1 (GRCh38, via VEP cache r112)"

# --- Ancestry: haplogroups (pipeline/ancestry.py) ----------------------------
# mtDNA via Haplogrep 3 (portable jar; trees auto-download on first classify), Y via
# Yleaf (hg38-native, fetched by uv from a pinned git commit). Depend on the genome
# only (not ClinVar), so output lives in ANCESTRY_DIR outside the dated snapshots.
HAPLOGREP_VERSION = "3.3.2"
HAPLOGREP_DIR = GENOMES_ROOT / "tools" / "haplogrep"
HAPLOGREP_JAR = HAPLOGREP_DIR / "haplogrep3.jar"
HAPLOGREP_URL = (f"https://github.com/genepi/haplogrep3/releases/download/"
                 f"v{HAPLOGREP_VERSION}/haplogrep3-{HAPLOGREP_VERSION}-linux.tar.gz")
HAPLOGREP_TREE = "phylotree-rcrs@17.3"
# Yleaf pinned to a commit for reproducibility; supplied to PGx-style `uv run --with`.
# Ephemeral `uv run --with` dependencies, pinned. They sit outside the project lockfile, so
# unpinned names resolved to whatever was newest on the day — a stage could change its output
# without a single line of this repo changing, and a cached result carried no record of what
# produced it. Versions below are the ones that resolved on 2026-09-07 and are known to work;
# bump deliberately.
WHATSHAP_SPEC = os.environ.get("WHATSHAP_SPEC", "whatshap==2.8")
PGX_STAR_DEPS = tuple(os.environ.get(
    "PGX_STAR_DEPS",
    "numpy==2.5.3,scipy==1.18.1,pysam==0.24.1,statsmodels==0.15.0").split(","))
PHARMCAT_DEPS = tuple(os.environ.get(
    "PHARMCAT_DEPS", "pandas==3.0.5,colorama==0.4.6,packaging~=24.1").split(","))

YLEAF_SPEC = ("yleaf @ git+https://github.com/genid/Yleaf"
              "@a2368eaaf08c0af53063514ed435d9e844b0ac1e")
ANCESTRY_DIR = GENOMES_ROOT / "ancestry"

# --- Repeat expansions (ExpansionHunter, pipeline/repeats.py) ----------------
# Genotypes the disease-associated short-tandem-repeat loci (Huntington HTT, Fragile X
# FMR1, ALS C9orf72, the ataxias …) that SNV calling is blind to. Illumina ships a
# macOS build + a GRCh38 catalog (Ensembl-named, matching our CRAMs). The binary +
# catalog are found by glob under EXPANSIONHUNTER_DIR (platform-agnostic).
EXPANSIONHUNTER_VERSION = "5.0.0"
EXPANSIONHUNTER_DIR = GENOMES_ROOT / "tools" / "expansionhunter"
REPEATS_DIR = GENOMES_ROOT / "repeats"

# --- Structural variants / CNV (Manta, pipeline/sv.py) -----------------------
# SNV/indel calling is blind to large deletions/duplications — a whole exon of BRCA1/2
# missing, the SMN1 exon-7 deletion behind spinal muscular atrophy, DMD exon-level CNVs.
# The monogenic panel flags exactly these genes as "not assessable from SNV calls". Manta
# calls SVs (DEL/DUP/INS/INV/BND) straight from the CRAM, but its workflow scripts need
# Python 2, so it runs inside a Docker container (biocontainers) — typically on the Linux
# box where Docker + the CRAMs live. SV calls depend only on the genome, not ClinVar, so
# output lives outside the dated snapshots in SV_DIR.
MANTA_VERSION = "1.6.0"
MANTA_IMAGE = os.environ.get(
    "MANTA_IMAGE", "quay.io/biocontainers/manta:1.6.0--h9ee0642_1")
SV_DIR = GENOMES_ROOT / "sv"
# Curated SV-relevant gene panel with GRCh38 (Ensembl-named, no "chr") coordinates —
# genes where large deletions/duplications are a real disease mechanism. Repo reference
# data; per-sample output in SV_DIR. A PASS SV overlapping one of these is flagged.
SV_GENES_FILE = Path(__file__).resolve().parent.parent / "panels" / "sv_genes.tsv"

# --- Callability / coverage (mosdepth, pipeline/callability.py) --------------
# Honest coverage reporting: a consumer report implies completeness, ours should not.
# mosdepth measures read depth genome-wide and per key gene so a "variant absent" call can
# be told apart from "region not adequately covered". A single static binary (no Python 2 /
# Docker), so it runs anywhere the CRAM + reference are. Reuses the SV gene panel
# (SV_GENES_FILE) for the per-gene callability table. Depends only on the genome, not ClinVar.
MOSDEPTH_VERSION = "0.3.10"
MOSDEPTH_DIR = GENOMES_ROOT / "tools" / "mosdepth"
CALLABLE_DIR = GENOMES_ROOT / "callable"
CALLABLE_MIN_DEPTH = 10      # a base is "callable" at >= this depth
CALLABLE_MIN_FRAC = 0.95     # flag a gene whose callable fraction falls below this

# --- Pharmacogenomics (PharmCAT, pipeline/pgx.py) ----------------------------
# PharmCAT = a Java jar + a Python preprocessor (its own deps, run via `uv run --with`
# so they stay out of this project's env). PGx results depend on the genome + PharmCAT
# version, NOT on ClinVar, so they live outside the dated snapshots in PGX_DIR.
# PharmCAT uses UCSC (chr1) contig naming, so pgx.py renames our Ensembl-named VCF first.
PHARMCAT_VERSION = "3.2.0"
PHARMCAT_DIR = GENOMES_ROOT / "tools" / "pharmcat"
PHARMCAT_JAR = PHARMCAT_DIR / f"pharmcat-{PHARMCAT_VERSION}-all.jar"
PHARMCAT_PIPELINE = PHARMCAT_DIR / "preprocessor" / "pharmcat_pipeline"
PHARMCAT_JAR_URL = (f"https://github.com/PharmGKB/PharmCAT/releases/download/"
                    f"v{PHARMCAT_VERSION}/pharmcat-{PHARMCAT_VERSION}-all.jar")
PHARMCAT_PREPROCESSOR_URL = (f"https://github.com/PharmGKB/PharmCAT/releases/download/"
                             f"v{PHARMCAT_VERSION}/pharmcat-preprocessor-{PHARMCAT_VERSION}.tar.gz")
PGX_DIR = GENOMES_ROOT / "pgx"

# --- CYP2D6 star-allele caller (Cyrius, feeds the PGx report) -----------------
# CYP2D6 can't be resolved from a plain VCF (CNV/hybrid/pseudogene), so PharmCAT returns
# "No Result". Cyrius genotypes CYP2D6 straight from the CRAM; pgx.py feeds its diplotype to
# PharmCAT as an outside call (-po) so CYP2D6 lands in the report. Cyrius is a cloned Python
# tool (deps via `uv run --with`); its GRCh38 data files + one hardcoded region literal are
# chr-prefixed, so setup strips them to our Ensembl contig naming. Runs where the CRAM is.
CYRIUS_DIR = GENOMES_ROOT / "tools" / "Cyrius"
CYRIUS_REPO = "https://github.com/Illumina/Cyrius.git"

# --- Disease panel (monogenic report, pipeline/panel.py) ---------------------
# Curated disease → gene(s) → inheritance panel. Lives in the repo (reference, not
# genomic data). Seeded from the ACMG SF list + classic carrier-screening genes;
# extend via Genomics England PanelApp. The AD/AR/XL mode drives carrier-vs-affected
# interpretation of a sample's actionable variants in the disease's genes.
PANEL_FILE = Path(__file__).resolve().parent.parent / "panels" / "monogenic.tsv"

# --- Penetrance / expressivity notes (S2, pipeline/penetrance.py) ------------
# Small curated gene -> plain-language penetrance/expressivity note, separate from the
# disease panel so it can be extended without touching panel membership. Optional and
# best-effort: genes not listed here fall back to a generic caveat at render time.
PENETRANCE_FILE = Path(__file__).resolve().parent.parent / "panels" / "penetrance.tsv"

# --- Single-SNP traits / wellness (pipeline/traits.py) -----------------------
# Curated rsID → genotype → interpretation rules with GRCh38 forward-strand coords
# (resolved via Ensembl). Reference data in the repo; per-sample output in TRAITS_DIR.
TRAITS_FILE = Path(__file__).resolve().parent.parent / "traits" / "traits.tsv"
TRAITS_DIR = GENOMES_ROOT / "traits"

# --- HLA risk tag SNPs (pipeline/hla.py) -------------------------------------
# Pragmatic HLA risk (B27, celiac DQ2.5/DQ8) via validated tag SNPs, same rules schema
# as traits. Full classical typing (arcasHLA/HIBAG) is a heavier follow-up.
HLA_FILE = Path(__file__).resolve().parent.parent / "panels" / "hla_tags.tsv"
HLA_DIR = GENOMES_ROOT / "hla"

# --- Full classical HLA typing (arcasHLA, pipeline/hla_typing.py) -----------
# 4-digit HLA-A/B/C/DPA1/DPB1/DQA1/DQB1/DRB1 alleles straight from the CRAM, superseding
# the tag-SNP proxies above. arcasHLA extracts chr6(+HLA alt/decoy) reads with samtools,
# then kallisto-pseudoaligns them against the IMGT/HLA reference to call genotypes — no
# Python-2 issue like Manta, but it ships via bioconda, so (like Manta) we run it
# containerized (biocontainers) rather than adding conda to this project's toolchain.
# Runs wherever the CRAM + Docker live (e.g. a Linux compute server, same as `sv`).
ARCASHLA_VERSION = "0.6.0"
ARCASHLA_IMAGE = os.environ.get(
    "ARCASHLA_IMAGE", "quay.io/biocontainers/arcas-hla:0.6.0--hdfd78af_2")
HLA_TYPE_DIR = GENOMES_ROOT / "hla_typing"
# Genes arcasHLA resolves well from short reads (its own recommended set).
HLA_TYPE_GENES = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

# --- Polygenic risk scores (PGS Catalog, pipeline/prs.py) --------------------
# Weighted dosage scores from PGS Catalog harmonized (GRCh38) scoring files. Small scores use
# the exact per-variant lookup + 1000G percentile calibration; large (genome-wide) scores use
# `plink2 --score`, which scales to millions of variants. PGS000018 is the Inouye CAD metaGRS
# (~1.7M variants) — the genome-wide demonstrator.
PRS_SCORES = [("PGS000805", "Type 2 diabetes"), ("PGS000001", "Breast cancer"),
              ("PGS000018", "Coronary artery disease (genome-wide)")]
# The trait panel proper (docs/polygenic-panel-design.md): one pinned PGS per trait, generated
# from PGS Catalog by objective criteria via panels/build_complex_traits.py. Used in preference
# to PRS_SCORES above, which stays as the offline/demo fallback when the panel is absent.
# Both the EFO term and the score ID are human-reviewed: fuzzy trait matching may *propose* a
# row but must never populate one (it silently mapped "carpal tunnel syndrome" to *body height*).
COMPLEX_PANEL = Path(__file__).resolve().parent.parent / "panels" / "complex_traits.tsv"
# Half the panel is genome-wide, and each such score needs its own cached 1000G reference-panel
# scoring run — intrinsic, not an artefact: for 18 of those 22 traits no small score exists at
# all, and for the other 4 dropping to one costs 0.06-0.095 AUROC. So the panel is meant to be
# run in tiers: this cap skips scores above it, letting the cheap ~22 land in minutes while the
# genome-wide half backfills separately. None = no cap (run everything).
PRS_MAX_VARIANTS: int | None = None
# S1 incidental-findings consent filter (docs/hardware-and-scope.md §5). Rows flagged
# `s1_gate` in the panel are severe *and* non-actionable, so they are withheld by default and
# shown only to a person who has opted in (`run.py consent`, consent.tsv). PRS_S1_GATE=0
# turns the gate off for everyone — full disclosure, a choice an operator makes for
# themselves, not one made for them.
# `check-releases` rescoring: a PRS report records the trait panel it was scored against, and
# the scheduled check re-scores the callsets whose report predates a panel change (a newly
# registered score, a corrected row). Only callsets that already have a report — it never
# starts polygenic scoring where nobody ran `prs`. PRS_AUTO_RESCORE=0 leaves the reports
# stale until `prs <sample>` is run by hand.
PRS_AUTO_RESCORE = os.environ.get("PRS_AUTO_RESCORE", "1").strip().lower() not in (
    "0", "false", "no", "off")
PRS_S1_GATE = os.environ.get("PRS_S1_GATE", "1").strip().lower() not in ("0", "false", "no", "off")
PGS_DIR = GENOMES_ROOT / "tools" / "pgs"
PGS_URL_TMPL = ("https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pid}/"
                "ScoringFiles/Harmonized/{pid}_hmPOS_GRCh38.txt.gz")
# Scores with more variants than this are scored with plink2 (genome-wide); smaller ones use
# the exact per-variant path. plink2 is a single binary, downloaded by setup.
PRS_PLINK_MIN = 50_000
PLINK2_DIR = GENOMES_ROOT / "tools" / "plink2"
PLINK2_URL_MAC = ("https://s3.amazonaws.com/plink2-assets/alpha7/"
                  "plink2_mac_arm64_20260504.zip")
PLINK2_URL_LINUX = "https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_latest.zip"
PRS_DIR = GENOMES_ROOT / "prs"
# 1000 Genomes (NYGC high-coverage, GRCh38) for PRS percentile calibration. Per-population
# allele frequencies (AF_EUR, AF_AFR, …) live in the INFO, so we remote-query just the PGS
# positions and compute the population score distribution analytically (HWE + CLT) rather
# than genotyping all 3202 individuals. The effect-allele frequencies are cached per PGS in
# PGS_DIR, so calibration is a one-time network cost. Contigs are UCSC-named (chr1).
KG_VCF_URL_TMPL = ("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
                   "1000G_2504_high_coverage/working/20201028_3202_phased/"
                   "CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}."
                   "filtered.shapeit2-duohmm-phased.vcf.gz")
# Local mirror of the above (see setup._mirror_1000g / run.py mirror-1kg): bcftools -R
# region-restricted queries over plain HTTP do one small range request per non-contiguous
# position, which is fast for a handful of positions but falls off a cliff at hundreds+
# scattered across a whole chromosome (measured: 5 positions <1s, 1000 positions >300s,
# despite the same server sustaining 17MB/s on a plain bulk GET) — a local file turns that
# into a fast local tabix seek. ~34.6GB for all 22 autosomes; every KG_VCF_URL_TMPL caller
# should go through `kg_vcf_source()` below rather than the template directly.
KG_VCF_LOCAL_DIR = GENOMES_ROOT / "refs" / "1000G_phased"
KG_POP = "EUR"   # superpopulation the primary percentile is matched to


def kg_vcf_source(chrom: str) -> str:
    """Local mirrored chromosome VCF if present, else the remote URL — same bcftools/plink2
    commands work unchanged against either. Requires the .tbi too: a partially-mirrored
    file (VCF without index) would fail every region query instead of falling back."""
    local = KG_VCF_LOCAL_DIR / f"chr{chrom}.vcf.gz"
    if local.exists() and Path(f"{local}.tbi").exists():
        return str(local)
    return KG_VCF_URL_TMPL.format(chrom=chrom)
# Genome-wide scores (millions of variants) can't use the analytic AF model above — instead we
# score the SAME PGS over the 1000G reference panel with plink2 and take an *empirical*
# percentile (count of reference individuals below the sample's score). PRS is additive, so we
# score per-chromosome and sum; the per-individual reference scores are cached per PGS
# (sample-independent → computed once, reused for the whole family). Best-effort: on any failure
# the genome-wide score still reports raw. The 3202-sample population map (sample → super-
# population) labels the EUR subset; it lives alongside the phased VCFs.
PRS_GENOMEWIDE_PERCENTILE = True
# The original 3202-sample ped+population file (20130606_g1k_3202_samples_ped_population.txt)
# has been removed from the 1000G FTP (404 as of 2026-07-02); its replacement in the same
# directory (1kGP.3202_samples.pedigree_info.txt) dropped the population columns entirely.
# Fall back to the phase-3 2504-unrelated-sample panel, which still carries population labels
# for all 5 superpopulations (~500-660 each) — the ~700 related samples added for the 3202
# high-coverage trio-phasing set are left unlabeled ('-') and are simply not used as the
# reference for a population percentile.
KG_PANEL_URL = ("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
                "integrated_call_samples_v3.20130502.ALL.panel")
KG_PANEL_FILE = PGS_DIR / "1kg_3202_population.txt"   # cached download

# --- Liftover references (hg19/GRCh37 -> GRCh38) ------------------------------
# For bringing an older-build callset over before `ingest` (see the README recipe). GATK
# LiftoverVcf needs a real target-build FASTA (not just a .dict) for ref-allele
# verification, with contig names matching the UCSC chain's target side ('chr1'). Rather
# than a second ~3GB hg38 download, setup rewrites our EXISTING REF_FASTA's contig headers
# to chr-prefixed into a cached copy, used only for liftover.
LIFTOVER_DIR = GENOMES_ROOT / "refs" / "liftover"
HG19TOHG38_CHAIN_URL = ("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/"
                        "hg19ToHg38.over.chain.gz")
HG19TOHG38_CHAIN = LIFTOVER_DIR / "hg19ToHg38.over.chain.gz"
HG38_CHR_FASTA = REFS_DIR / "GRCh38.chr_named.fa"   # chr-prefixed copy, LiftoverVcf-only

# --- Read-backed phasing (WhatsHap, pipeline/phasing.py) ---------------------
# Resolves in-cis vs in-trans for multi-het disease genes (the panel's "possible compound
# het — review" flag) by phasing the CRAM reads. WhatsHap is pip-installable, supplied at
# run time via `uv run --with`. Output depends on the genome, not ClinVar, so it lives in
# PHASE_DIR outside the dated snapshots. Runs where the CRAM is.
PHASE_DIR = GENOMES_ROOT / "phase"

# --- Share (pipeline/share.py) -------------------------------------------------
# Copies of a person's VCF or CRAM made for handing to someone else, with the working files'
# command lines, paths and internal ids taken out of the header.
SHARE_DIR = GENOMES_ROOT / "share"

# --- Force-call (all-sites genotyping at fixed positions, pipeline/forcecall.py) ---
# HaplotypeCaller emits only *variant* sites, so a marker absent from a sample's VCF is
# ambiguous: genuine homozygous-reference, or simply never looked at. The lookup engines
# (traits/HLA single-SNP, small PRS) need the genotype at *fixed* positions including
# hom-ref. This stage force-calls the union of those positions off the CRAM (GATK
# `--alleles` + EMIT_ALL_ACTIVE_SITES), turning an *assumed* hom-ref into an *evidenced*
# one and surfacing genuine no-calls (low coverage) instead of silently scoring them as
# reference. Genome-wide PRS stays on the plink2 path (millions of sites). Output depends
# on the genome, not ClinVar, so it lives in FORCECALL_DIR outside the dated snapshots;
# the sample-independent sites VCF is cached (hashed by its inputs) under FC_SITES_DIR.
FORCECALL_DIR = GENOMES_ROOT / "forcecalled"
FC_SITES_DIR = FORCECALL_DIR / "sites"
# A scored variant joins the force-call union only when its source list has at most this
# many variants — i.e. the curated/small lists, never the genome-wide demonstrator.
FC_MAX_LIST = PRS_PLINK_MIN

# --- AlphaMissense (deep missense pathogenicity, pipeline/alphamissense.py) ---
# DeepMind's precomputed pathogenicity for ~71M possible missense variants (GRCh38),
# independent of ClinVar. Large (~1GB); bgzipped + tabix-indexed by setup. Reports
# enrich themselves with it when present and degrade gracefully when absent.
ALPHAMISSENSE_FILE = REFS_DIR / "AlphaMissense_hg38.tsv.bgz"
ALPHAMISSENSE_URL = ("https://storage.googleapis.com/dm_alphamissense/"
                     "AlphaMissense_hg38.tsv.gz")

# --- gnomAD gene constraint (LOEUF/pLI, pipeline/gnomad_constraint.py) --------
# Per-gene loss-of-function intolerance from gnomAD v4.1 (~734k exomes). Strengthens the
# ACMG PVS1 criterion: a null variant only earns full PVS1 in a gene whose haploinsufficiency
# is plausible (LOF-constrained). `setup` downloads the per-transcript constraint metrics
# (~95MB) and reduces it to a small gene-keyed table (gene → pLI, LOEUF) keeping the
# MANE-select transcript per gene. Build-agnostic (keyed by gene symbol). Constraint is
# meaningful for haploinsufficient (dominant) genes; recessive LOF doesn't require it.
GNOMAD_CONSTRAINT_VERSION = "v4.1"
GNOMAD_CONSTRAINT_FILE = REFS_DIR / "gnomad_constraint_by_gene.tsv"   # gene\tpLI\tLOEUF
GNOMAD_CONSTRAINT_URL = ("https://storage.googleapis.com/gcp-public-data--gnomad/"
                         "release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv")
# A gene is "LOF-constrained" (haploinsufficiency plausible) at either threshold. pLI≥0.9 is
# the ClinGen-SVI operationalization; LOEUF<0.6 catches constrained genes pLI alone misses.
GNOMAD_PLI_CONSTRAINED = 0.9
GNOMAD_LOEUF_CONSTRAINED = 0.6

# --- SpliceAI (deep splice-disruption prediction, pipeline/spliceai.py) -------
# Illumina's precomputed delta scores for every possible SNV (and curated indels) in
# genic GRCh38 regions — a splice-impact signal orthogonal to AlphaMissense's missense
# pathogenicity, so it covers the variants AlphaMissense can't (splice-region, synonymous-
# but-splicing, intronic). Same precomputed-join pattern: a tabix point-query over the VCF,
# best-effort, reports degrade gracefully when absent. We default to the *masked* files
# (scores zeroed where a change is not splice-relevant — fewer false positives, the
# clinically recommended set). NOT auto-downloaded: the precomputed VCFs are distributed via
# Illumina BaseSpace behind a login, so there is no anonymous URL. Obtain them once and drop
# them in REFS_DIR (see SPLICEAI_OBTAIN below); `setup` tabix-indexes them if present.
SPLICEAI_SNV_FILE = REFS_DIR / "spliceai_scores.masked.snv.hg38.vcf.gz"
SPLICEAI_INDEL_FILE = REFS_DIR / "spliceai_scores.masked.indel.hg38.vcf.gz"
SPLICEAI_OBTAIN = (
    "SpliceAI precomputed scores are on Illumina BaseSpace (login required): "
    "https://basespace.illumina.com/s/otSPW8hnhaZR (Genome Annotations). Download "
    "the GRCh38 masked SNV + indel VCFs and place them in the refs dir as "
    "spliceai_scores.masked.{snv,indel}.hg38.vcf.gz, then re-run `setup`.")
# Illumina's recommended delta-score cut-offs: ≥0.8 high-precision, ≥0.5 recommended
# (likely splice-altering), ≥0.2 high-recall (possible). Below 0.2 → no predicted effect.
SPLICEAI_DS_HIGH = 0.8
SPLICEAI_DS_RECOMMENDED = 0.5
SPLICEAI_DS_LOW = 0.2

# --- Non-coding regulatory effect (GPU, pipeline/interpret.py) — I3 -----------
# The pipeline's only *on-demand GPU* in-silico signal: AlphaMissense/SpliceAI are joins over
# public precomputed tables, but non-coding regulatory effect isn't precomputed genome-wide
# anywhere usable, so we run a sequence-to-function model ourselves on a GPU box. It
# predicts thousands of regulatory tracks from a long context; the per-variant effect is the
# change in those tracks between ref and alt. Fills VEP+ClinVar's biggest blind spot
# (promoter/enhancer/UTR variants). The box is kept *stateless*: the Mac (which owns the
# reference) extracts the ref/alt windows and ships a compact job over SSH; the box runs the
# baked worker image (built from pipeline/interpret_worker.py) and returns scores. Scores depend
# only on (variant, model) — not on snapshot/sample — so each model's table is a genome-level
# cache that accumulates and is re-joined as the family's variant set grows.
#
# Two interchangeable backends; pick with INTERPRET_MODEL. Enformer is the default (lighter,
# 196 kb context); Borzoi is its RNA-seq successor (524 kb context, more tracks). Cache,
# calibration, targets and window length are all per-model (their deltas aren't comparable).
INTERPRET_MODEL = os.environ.get("INTERPRET_MODEL", "enformer")
# Each backend gets its own Docker image: enformer-pytorch pulls torch 2.12 (cu13) while Borzoi
# is proven on the base image's torch 2.3.1 (the newer torch breaks Borzoi's weight load), so
# they can't share one image. The worker script is shared — it imports each model's package
# lazily, only in that model's branch. Built on the box from pipeline/Dockerfile.{model}.
INTERPRET_MODELS = {
    "enformer": {"hf": "EleutherAI/enformer-official-rough", "seq_len": 196_608,
                 "image": "enformer-worker:latest",
                 "targets": ("https://raw.githubusercontent.com/calico/basenji/master/"
                             "manuscripts/cross2020/targets_human.txt")},
    "borzoi":   {"hf": "johahi/borzoi-replicate-0", "seq_len": 524_288,
                 "image": "borzoi-worker:latest",
                 "targets": ("https://raw.githubusercontent.com/calico/borzoi/main/"
                             "examples/targets_human.txt")},
}
_M = INTERPRET_MODELS[INTERPRET_MODEL]
INTERPRET_MODEL_HF = _M["hf"]
ENFORMER_SEQ_LEN = _M["seq_len"]        # the model's fixed input length (the Mac extracts windows of this)
INTERPRET_DIR = GENOMES_ROOT / "interpret"
ENFORMER_SCORES_FILE = INTERPRET_DIR / f"{INTERPRET_MODEL}_scores.tsv"   # per-model, by (chrom,pos,ref,alt)
# Remote GPU box: where the worker image lives and runs. SSH target + a scratch dir on the
# box for the per-run job/result files, and the persistent HF/model cache (downloaded once).
INTERPRET_SSH_HOST = os.environ.get("INTERPRET_SSH_HOST", "user@gpu-box.local")
INTERPRET_REMOTE_DIR = os.environ.get("INTERPRET_REMOTE_DIR", "~/enformer/jobs")
INTERPRET_CACHE_DIR = os.environ.get("INTERPRET_CACHE_DIR", "~/enformer")  # mounts to /cache
# The selected model's worker image (per-model — see INTERPRET_MODELS above).
ENFORMER_IMAGE = os.environ.get("ENFORMER_IMAGE", _M["image"])
ENFORMER_BATCH_SIZE = int(os.environ.get("ENFORMER_BATCH_SIZE", "4"))
# Human target metadata (index → assay description) for the selected model, so the report can say
# *which* regulatory signal moved (e.g. "CAGE:brain") rather than a bare track index. Fetched
# best-effort to the box cache on first run; the worker degrades to empty descriptions if absent.
# For Borzoi this file is load-bearing beyond labels: its `strand_pair` column is what lets the
# worker swap +/- track pairs when averaging the reverse-complement pass. Without that column the
# worker disables RC averaging rather than average stranded tracks against their opposite strand,
# so a targets URL that drops `strand_pair` silently costs the RC robustness pass.
ENFORMER_TARGETS_URL = _M["targets"]
# Which VEP consequences to send to Enformer: the non-coding/regulatory classes that the
# missense (AlphaMissense) and splice (SpliceAI) predictors can't score. Intergenic is
# excluded by default — too numerous and far from any regulatory target to be informative.
INTERPRET_CONSEQUENCES = (
    "regulatory_region_variant", "TF_binding_site_variant", "TFBS_ablation",
    "TFBS_amplification", "5_prime_UTR_variant", "3_prime_UTR_variant",
    "upstream_gene_variant", "downstream_gene_variant",
    "non_coding_transcript_exon_variant", "mature_miRNA_variant",
)
# Keep each run GPU-affordable and focused: cap the variant set, and prefer rare ones (common
# regulatory variants are unlikely to be the family's interesting findings). A carried variant
# at/under this gnomAD AF (or AF-absent) is eligible. ~0.5s/variant means low-hundreds is fine.
INTERPRET_MAX_VARIANTS = int(os.environ.get("INTERPRET_MAX_VARIANTS", "256"))
INTERPRET_MAX_AF = float(os.environ.get("INTERPRET_MAX_AF", "0.05"))
# Effect bands for delta_max (largest |alt-ref| track change at the central bins). Enformer
# deltas are *not* normalized to [0,1] like SpliceAI. These raw cut-offs are the *fallback*
# used only until a calibration exists (see below); once calibrated, banding is by percentile.
ENFORMER_DELTA_MODERATE = float(os.environ.get("ENFORMER_DELTA_MODERATE", "1.0"))
ENFORMER_DELTA_HIGH = float(os.environ.get("ENFORMER_DELTA_HIGH", "3.0"))

# --- Calibration: empirical percentile of Δmax vs a benign null (interpret calibrate) -------
# A raw Δmax is uninterpretable on its own; a *percentile against common (≈benign) variants of
# the same type* is. `interpret calibrate` scores a sample of common non-coding variants (high
# gnomAD AF — overwhelmingly benign for regulatory disruption) to build a null distribution of
# Δmax, then every scored variant is reported as its percentile against that null. Stratified by
# SNV vs indel because an indel's Δ scales with its length (pooling them would be wrong). Bands
# then mean "perturbs regulation more than N% of common variants of its kind" — the same idiom
# as the PRS percentile. The calibration is genome/model-level (sample-independent), cached here.
ENFORMER_CALIBRATION_FILE = INTERPRET_DIR / f"{INTERPRET_MODEL}_calibration.json"
INTERPRET_COMMON_AF = float(os.environ.get("INTERPRET_COMMON_AF", "0.2"))  # null = AF ≥ this
# Variants to score per stratum (SNV / indel) when building the null. Bigger = tighter tail
# estimate but slower (~1 s/variant); a few hundred is enough for 95th/99th-pctl triage.
INTERPRET_CALIBRATION_N = int(os.environ.get("INTERPRET_CALIBRATION_N", "400"))
# Percentile cut-offs for the calibrated bands (vs the benign null).
ENFORMER_PCTL_MODERATE = float(os.environ.get("ENFORMER_PCTL_MODERATE", "95"))
ENFORMER_PCTL_HIGH = float(os.environ.get("ENFORMER_PCTL_HIGH", "99"))

# Set only by `run.py demo`. When non-empty, every report written and every dashboard page
# built carries it as a callout directly under the title, so a page from the synthetic demo
# family can never be mistaken for a result about a real person.
DEMO_NOTICE = ""

# --- Notification hook (pipeline/notify.py) -----------------------------------
# Optional. A command run with the scan's one-line summary as its last argument (no shell).
# Unset = no notification. e.g. NOTIFY_CMD="ntfy publish genome"
NOTIFY_CMD = os.environ.get("NOTIFY_CMD", "")
NOTIFY_TIMEOUT = int(os.environ.get("NOTIFY_TIMEOUT", "30"))


# Haplogrep quality floor below which an mtDNA haplogroup call is reported as unresolved
# rather than as a result. A low-quality call lands on the tree root (rCRS, H2a2a1) — which
# is what "I found no informative variants" looks like, not a finding. A callset lifted from
# hg19 is the typical cause: every MT variant is offset 1-2 bp (the two builds use different
# mitochondrial reference sequences and the UCSC chain maps chrM naively), which produces a
# confident-looking H2a2a1 that is an artefact of the lift, not a result.
MT_MIN_QUALITY = float(os.environ.get("MT_MIN_QUALITY", "0.85"))

# Trio segregation / transmission-phased compound het (M2, pipeline/segregation.py).
# Runs off the PARENTS' CRAMs (the child needs only a VCF), so it lives outside the snapshots.
SEGREGATION_DIR = GENOMES_ROOT / "segregation"

# --- Trio de-novo / segregation (M2, pipeline/denovo.py) ---------------------
# Runs off the child's normalized VCF + the PARENTS' CRAMs, so it lives outside the dated
# snapshots like phase/force-call. Runs where the parents' CRAMs are.
DENOVO_DIR = GENOMES_ROOT / "denovo"

# Child-side candidate filters. A true de-novo SNV is heterozygous, balanced, and at
# ordinary depth; each FP mode skews one of those — mismapped/collapsed repeats pile up
# depth, post-zygotic mosaicism drops the allele fraction, and a low-GQ call is simply
# not evidence. Hom-alt is excluded outright: two independent de-novo hits at one site
# is vanishingly rare, so `1/1` means the parents were miscalled, not that the child is new.
DN_MIN_GQ = int(os.environ.get("DN_MIN_GQ", "20"))
DN_MIN_DP = int(os.environ.get("DN_MIN_DP", "10"))
# Upper depth bound as a multiple of the child's median het depth (collapsed-paralog guard).
DN_MAX_DP_FACTOR = float(os.environ.get("DN_MAX_DP_FACTOR", "2.5"))
# 0.40, not 0.25. A floor that low lets through a confirmed set many times the 40-150 SNVs a
# genome is expected to carry, with the survivors clustered just above the floor rather than
# at the 0.5 a germline het sits at. That shape is artefact leakage, not mosaicism. When
# changing it, check the confirmed SNV count against that prior.
DN_MIN_VAF = float(os.environ.get("DN_MIN_VAF", "0.40"))
# Indels carry their own artefact load (reference bias, slippage, alignment ambiguity) and get
# their own floor so they can be tightened without touching SNVs.
DN_MIN_VAF_INDEL = float(os.environ.get("DN_MIN_VAF_INDEL", "0.40"))
DN_MAX_VAF = float(os.environ.get("DN_MAX_VAF", "0.75"))
# Short insertions/deletions of a single repeated base in a homopolymer run are polymerase
# slippage, not mutation, and without this filter typically the largest artefact class after
# clustering. True de-novo events are ~90% SNV; an
# unfiltered set that is mostly indel has the ratio inverted.
DN_DROP_HOMOPOLYMER_INDELS = os.environ.get(
    "DN_DROP_HOMOPOLYMER_INDELS", "1") not in ("0", "false", "")
# Autosomes only by default: a male child's chrX/chrY are hemizygous (no meaningful "het"),
# and MT is maternally inherited with heteroplasmy — both need their own logic, not this one.
DN_AUTOSOMES_ONLY = os.environ.get("DN_AUTOSOMES_ONLY", "1") not in ("0", "false", "")

# Parent-side confirmation, read off the CRAM (NOT the parent's VCF). This is the whole
# point of the stage: "absent from the parent's variant-only VCF" conflates true hom-ref
# with never-called, and across a cross-caller/cross-build comparison the latter dominates.
DN_MIN_PARENT_DP = int(os.environ.get("DN_MIN_PARENT_DP", "10"))
# Alt-supporting reads tolerated in a parent before the site is called inherited. Not 0:
# a single alt read at 40x is indistinguishable from sequencing error.
DN_MAX_PARENT_ALT = int(os.environ.get("DN_MAX_PARENT_ALT", "1"))
# Second, CALLER-INDEPENDENT parental check: a raw pileup (no local assembly, no genotype
# model). Force-calling the parents with HaplotypeCaller is not an independent test at sites
# the child's caller found and HC systematically cannot — HC returns a confident, *wrong*
# hom-ref there. The pileup catches inherited sites that the force-call pass alone would
# have called de-novo.
DN_PILEUP_MIN_ALT = int(os.environ.get("DN_PILEUP_MIN_ALT", "3"))
DN_PILEUP_MIN_MQ = int(os.environ.get("DN_PILEUP_MIN_MQ", "20"))
DN_PILEUP_MIN_BQ = int(os.environ.get("DN_PILEUP_MIN_BQ", "20"))
# Real de-novo events are independent point mutations scattered across the genome; several
# landing within a kilobase of each other is not biology, it is one badly-behaved locus
# (mismapping, a collapsed paralog, a misassembled region) emitting a burst of calls. Scattered
# at random, essentially none would sit within 1 kb of another, so a large clustered fraction
# is the signature of this — often the single largest residual artefact class.
DN_CLUSTER_WINDOW = int(os.environ.get("DN_CLUSTER_WINDOW", "1000"))

# Published germline de-novo SNV rate per generation — the sanity scale for the output.
# Not a filter: a count far outside it means the pipeline is leaking, and the report says so.
DN_EXPECTED_MIN = int(os.environ.get("DN_EXPECTED_MIN", "40"))
DN_EXPECTED_MAX = int(os.environ.get("DN_EXPECTED_MAX", "150"))

# --- Report-producing directories (render + serve, this bug class) -----------
# Every stage directory that can hold rendered .md reports for a person (or, for
# INTERPRET_DIR, per-snapshot). Single source of truth for:
#   - `run.py cmd_render`: which roots to walk for markdown -> HTML conversion.
#   - `pipeline/serve.py _ALLOWED`: which subdirs of GENOMES_ROOT are web-servable.
# Every report lives under REPORTS_DIR (pipeline/reportpaths.py lays it out: md/, md-history/,
# html/). Stage directories hold working data only, so nothing else is ever a report dir —
# and only REPORTS_DIR/html is web-servable (see serve.py).
REPORT_DIRS = (REPORTS_DIR,)

# --- Classification / diff thresholds ----------------------------------------
# ClinVar classification TERMS we treat as "actionable". A record's CLNSIG is a compound of
# terms (`Pathogenic/Pathogenic,_low_penetrance|risk_factor`), so these are matched per
# term by pipeline/clinsig.py — never compare a whole CLNSIG string against this set.
ACTIONABLE_TERMS = {"Pathogenic", "Likely_pathogenic"}
# Significance values considered "uncertain" — a transition out of these into an
# actionable bucket is the highest-priority DELTA event.
UNCERTAIN_SIG = {"Uncertain_significance", "Conflicting_classifications_of_pathogenicity",
                 "not_provided", ""}
# gnomAD allele-frequency boundary separating "rare" from "common".
RARE_AF_THRESHOLD = 0.001
# Allele balance: a germline het carries the alt allele on ~50% of reads, so a het call
# whose alt fraction (AD) is below this is flagged "low allele fraction — review" in the
# reports. Flag, never drop: at ~45x a true het falls below 25% about 1 time in 2,500, and
# those misses concentrate at low-depth sites. Chosen so that the usual artefacts behind
# reported P/LP calls — homopolymer-slippage frameshifts, mismapped-paralog calls — fall
# well below it and true heterozygous carriers well above. Paralog clusters whose alt fraction sits near 50%, or right
# at the line, are NOT caught by this.
LOW_VAF_THRESHOLD = float(os.environ.get("LOW_VAF_THRESHOLD", "0.25"))

# ClinVar review-status ranking (star rating), used to detect confidence upgrades.
REVIEW_STATUS_RANK = {
    "no_assertion_provided": 0,
    "no_assertion_criteria_provided": 0,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "reviewed_by_expert_panel": 3,
    "practice_guideline": 4,
}

# --- Delta triage (multi-factor) ---------------------------------------------
# A changed variant's triage score is a weighted sum of five auditable factors,
# each normalized to roughly [-1, 1]. All weights and cut-offs live here so the
# logic stays adjustable (and tunable by the retrospective validation harness).
#   score = Σ weight[f] * factor[f]   →   tier by TRIAGE_CUTOFFS
TRIAGE_WEIGHTS = {
    "acmg": 0.35,       # reasoned ACMG tier + ClinVar divergence (panel genes only)
    "model": 0.20,      # in-silico distance from threshold, gated by consequence class
    "freq": 0.20,       # population rarity (gnomAD AF)
    "stability": 0.10,  # how steady the call has been across snapshots (Phase 2)
    "gene": 0.15,       # panel disease gene / LOF-constrained gene context
}
# Tier cut-offs on the composite score (descending). Calibrated to the factor ranges:
# a coding variant is gated out of the model factor, so the realistic max for the canonical
# "new pathogenic, rare, on-panel" finding is ~0.8 — cut-offs sit below that so each tier is
# reachable. `actionable_now` additionally requires real clinical context (a panel gene or a
# ClinVar-actionable call) — a high score from model+freq alone on an unknown gene is
# "monitor", not "act". These are the Feature-3 (validation harness) tuning surface.
TRIAGE_CUTOFFS = {"actionable_now": 0.60, "monitor": 0.35, "log_only": 0.15}
# Reasoned-ACMG-tier → factor contribution.
TRIAGE_ACMG_SCORE = {"Pathogenic": 1.0, "Likely pathogenic": 0.7, "VUS": 0.3,
                     "Likely benign": 0.1, "Benign": 0.0}
# Flat bonus when the reasoned call diverges from ClinVar (a candidate reclassification).
TRIAGE_DIVERGENCE_BONUS = 0.3
# How many recent snapshots define the stability window (Phase 2 stability factor).
TRIAGE_STABILITY_WINDOW = 5

# --- Incremental scan --------------------------------------------------------
# When a sample's genome + VEP cache are unchanged but ClinVar advanced, re-annotate ONLY the
# carried variants at ClinVar-changed coordinates (the rest are copied forward) instead of a
# full VEP pass. Correct because ClinVar is attached `type=exact`. Set to 0 to fall back to a
# full re-annotation whenever ClinVar changes.
SCAN_CLINVAR_DELTA = os.environ.get("SCAN_CLINVAR_DELTA", "1") not in ("0", "false", "")
