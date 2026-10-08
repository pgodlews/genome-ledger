# Hardware Profiles, Scope & Framing

What the pipeline is for, which input files unlock which features, what hardware each
workload needs, and the legal and ethical framing it ships under.

---

## 1. Core Philosophy & Value Proposition

The pipeline is *not* merely "another variant caller" — alignment and variant calling are largely commodity tools. Its unique value lies in **orchestration and time-aware longitudinal analysis**:

1. **Time-Aware Reclassification Engine**: Pairs dated ClinVar snapshots with a continuous DELTA triage system (`snapshot → annotate → delta`). As clinical knowledge evolves over years, the pipeline automatically flags newly pathogenic or downgraded variants across the family without re-running upstream compute.
2. **Pedigree-Scale & Family-Aware**: Designed around multi-generational trios/quads, compound heterozygosity phasing, de novo mutation verification, and carrier screening.
3. **Local Sovereignty & Privacy-First**: Runs entirely on hardware you own (from a laptop to a homelab server), keeping sensitive genomic data off third-party commercial clouds.
4. **Decade-Scale Stewardship**: Structured to preserve irreplaceable raw data and analytical contracts across 20+ years, surviving hardware ISA and software stack transitions.

---

## 2. Genomic File Hierarchy: What File Unlocks What

A common point of confusion in genomics is knowing which file types are required for which analyses. The pipeline is designed to be **modular**; you do not need every file type to use the framework.

| File Type | Format / Ext | Typical Size (30x WGS) | What Pipeline Capabilities It Unlocks | Required Hardware Profile |
|---|---|---|---|---|
| **Raw Reads** | FASTQ (`.fq.gz`) | ~60–90 GB / sample | • Primary alignment to reference genomes<br>• Variant calling from scratch<br>• Re-alignment to new reference builds (e.g. GRCh38 ➔ T2T-CHM13)<br>• Novel structural variant discovery & repeat expansions<br>• Discovery of unmapped reads (viral, bacterial, microbiome) | **Profile C** (CUDA GPU Box recommended) or large CPU cluster |
| **Aligned Reads** | CRAM / BAM (`.cram`, `.crai`) | ~15–25 GB / sample | • Read depth & coverage validation<br>• Visual pileup inspection (Samtools, IGV)<br>• *De novo* mutation validation via parent-child allele balance (`denovo.py`)<br>• Structural variant read-support confirmation (`sv.py`)<br>• Callable loci analysis & gap filling | **Profile B** (Workstation) or **Profile A** |
| **Called Variants** | VCF / gVCF (`.vcf.gz`, `.tbi`) | ~150–800 MB / sample | • **Core Downstream Pipeline (~90% of features):**<br>  - VEP functional consequence annotation (`annotate.py`)<br>  - ClinVar dated snapshot ingestion & DELTA reclassification engine (`incremental.py`, `diff.py`)<br>  - ACMG 5-tier classification (`acmg.py`)<br>  - Disease panel screening & carrier risk (`panel.py`)<br>  - Polygenic risk scores for a 44-trait panel (`prs.py`)<br>  - Pharmacogenomics CPIC dosing guidelines (`pgx.py`)<br>  - Single-SNP traits & HLA tag SNPs (`traits.py`, `hla.py`) | **Profile A** (Any modern laptop or desktop, CPU only) |
| **Analytical Warehouse** | DuckDB (`.duckdb`) | ~50–200 MB / family | • Sub-second SQL queries across whole family pedigree (`run.py query`)<br>• Longitudinal diffs between ClinVar releases<br>• Instant report generation (Markdown & HTML)<br>• Automated carrier compatibility & Mendelian inheritance checks | **Profile A** (Lightweight CPU) |

> [!TIP]
> **External User Takeaway**: If you already have your VCF from a sequencing provider or clinical test, you **do not** need raw FASTQ or CRAM files to run the entire downstream interpretation, ACMG classification, ClinVar DELTA reclassification, and PRS reports. You only need Profile A hardware!

---

## 3. Hardware Requirements

The framework separates its workloads into **clear compute profiles**; pick the lightest one that covers what you want to do.

### Profile A: Downstream Interpretation & Reclassification (Lightweight / Laptop)
*For running annotation, ClinVar snapshots, ACMG scoring, PRS, and family reports from existing VCFs.*
- **CPU**: 4–8 cores (Apple Silicon M1/M2/M3/M4, or Intel/AMD x86_64).
- **RAM**: 8 GB minimum (16 GB recommended).
- **Disk**: 20–50 GB SSD storage (for DuckDB, reference FASTA, and ClinVar release archives).
- **GPU**: **None required (Zero GPU)**.
- **OS**: macOS or Linux.

### Profile B: Standard Family Processing & De Novo Segregation (Workstation)
*For running local VEP annotation caches, bcftools normalization, and CRAM pileup inspection for de novo trios.*
- **CPU**: 8–16 cores.
- **RAM**: 32 GB recommended.
- **Disk**: 200 GB+ fast NVMe (accommodates VEP local cache files and multiple family CRAM files).
- **GPU**: Not required.
- **OS**: Linux or macOS.

### Profile C: High-Throughput Primary Processing (CUDA GPU Box)
*For raw FASTQ alignment (BWA-MEM / fq2bam) and DeepVariant germline variant calling.*
- **Hardware Requirements**:
  - **GPU**: Dedicated NVIDIA GPU with CUDA support, minimum **16 GB VRAM** (e.g. RTX 3080 16GB, RTX 3090, RTX 4080/4090, RTX A4000/A5000, Tesla T4/L4/A10G).
  - **Host CPU**: Modern x86_64 CPU with ≥16 threads.
  - **Host RAM**: 32–64 GB system RAM.
  - **Scratch Space**: ≥500 GB fast NVMe disk for intermediate BAM and sharded temporary sort files.
- **Software Stack**:
  - NVIDIA Container Toolkit / Docker / Apptainer.
  - Accelerated calling engine: **NVIDIA Clara Parabricks** (`pbrun fq2bam`, `pbrun deepvariant`) or Google DeepVariant GPU container.
- **Performance Impact**:
  - A 30x whole genome FASTQ alignment and variant call drops from ~18–24 hours on CPU down to **~40–60 minutes** on a modern CUDA GPU.

A GPU box is used purely as a stateless compute node: point `INTERPRET_SSH_HOST` (and the
related settings in `pipeline/config.py`) at it, or leave them unset and stay on Profile A/B.

---

## 4. Legal, Licensing & Framing

- **License**: MIT (see `LICENSE`).
- **Research and Educational Use Only.** This software is not a certified medical device and
  is not intended or validated for clinical diagnosis, treatment decisions, or patient care.
  Genetic variant classifications and polygenic risk predictions are research-grade
  estimations. Any medically actionable finding must be verified independently by a
  board-certified clinical genetics laboratory and discussed with a qualified genetic
  counselor or physician. Provided AS-IS without warranty of any kind.
- **Regulatory posture**: an open-source research and bioinformatic automation tool (in the
  same space as GEMINI, OpenCRAVAT, or slivar). No diagnostic claims.

---

## 5. Consent & Incidental Findings Filter (S1)

Not everyone in a family wants to know everything their genome says. The S1 filter gates
categories of findings a person may not wish to learn (adult-onset untreatable conditions,
late-onset neurodegenerative risk such as *APOE*-e4 homozygosity, carrier status) behind an
explicit per-person, per-category opt-in recorded in `consent.tsv`
(`consent.example.tsv` shows the format; `pipeline/consent.py` is the gate).

These findings are **withheld by default**: a report says that something was held back and
how to opt in, and `run.py consent` records the grant per person. An operator who wants full
disclosure for everyone sets `PRS_S1_GATE=0`. The synthetic demo family is shown everything.

**Withdrawing is immediate.** `run.py consent <sample> <category> no` deletes that person's
reports whose content depended on the grant — the repeat-expansion report and, while the S1
gate is on, the polygenic report — in every format (markdown, HTML, JSON), and rebuilds their
dashboard page so no summary of them remains. Each of those reports records the consent
state it was written under, so a copy that survives some other way (a restore, a sync from
another machine) is regenerated rather than reused the next time the stage runs.

See also [`stewardship-and-succession.md`](stewardship-and-succession.md) for the family
charter prompts (right not to know, consent for minors, custody handoff).
