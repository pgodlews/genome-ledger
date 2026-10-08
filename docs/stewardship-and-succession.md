# Family Genomic Trust & Succession — a template

Most documentation is about *features*. This one is about **custody of an asset meant to
outlive the person who set it up**: a family's genomes, and the pipeline that reads them,
handed to the next generation and kept useful for decades.

The goal: the genomes and this pipeline become a durable family archive that **gains value**
as databases improve, the pedigree deepens and tools advance — preserved well enough that a
descendant decades from now can still read the data, re-run the analysis, add a
newly-sequenced relative, and learn more than anyone could today.

Copy this file into your own records and fill it in. The technical sections are
prescriptive. **Governance and consent (§6) are prompts, not answers — they are family
decisions, not engineering defaults.**

---

## 1. The asset and its threat model

**What is irreplaceable:** each person's **raw sequencing reads** (FASTQ). Once someone has
died, their genome can never be sampled again. Everything downstream — alignments, VCFs,
annotations, reports — can be *recovered* by re-running the pipeline on the reads. The reads
cannot be recovered from anything. **Preserving the raw reads is the single most important
job.**

Failure modes a generational handoff must survive:

| Threat | Mitigation (section) |
|---|---|
| Data loss / silent bit-rot | redundant backup + integrity scrubbing (§2) |
| Format obsolescence | open formats + archived specs + migration (§3) |
| **CPU architecture obsolescence** | source-first, data-first; containers are a convenience, not a guarantee (§3) |
| Environment rot (package managers, download URLs vanish) | pinned container + rebuildable source (§3) |
| Knowledge loss (the heir can't operate it) | successor runbook + local UI + scheduled runs (§4) |
| Custody loss (the keys die with one person) | decentralized mirrors + succession plan (§5) |
| Consent / privacy breach across people and generations | family charter (§6) |

---

## 2. Preservation (do this first — it protects the irreplaceable)

- **3-2-1 backup of the raw reads and the reference:** at least 3 copies, on at least 2
  kinds of media, at least 1 off-site. For example: a primary copy on a NAS, a second on
  cold or encrypted storage, and a third at a *different physical location* — a relative's
  home, a safe-deposit box, or an encrypted cloud archive. A 30x whole genome is roughly
  50–100 GB of compressed FASTQ, and the archive grows with every relative sequenced —
  budget for it.
- **Integrity scrubbing:** silent bit-rot is the quiet killer of decade-scale archives. Use a
  checksumming filesystem (ZFS/Btrfs) with scheduled scrubs, and/or generate `par2` parity
  files. The **sha256 the pipeline records at ingest** is the tripwire — verify it on a
  schedule, not only at ingest.
- **⚠️ Archive the reference *with* the data.** CRAM is reference-*compressed*: a CRAM file
  cannot be decoded without the **exact** reference FASTA it was written against. Store the
  GRCh38 primary assembly (and any other reference you use) in the cold archive next to the
  CRAMs, or a future heir inherits unreadable files. Raw FASTQ has no such dependency —
  another reason it is the true master copy.
- **What goes in the cold archive, in priority order:** ① raw FASTQ + sha256 manifest,
  ② the reference FASTA(s), ③ the per-sample VCFs, ④ the family graph and sample manifest,
  ⑤ the pinned pipeline (source + container, §3), ⑥ this document and the format specs (§3).
  ③–④ can be regenerated from ①–② and ⑤, but they are cheap to keep.

---

## 3. Reproducibility across decades — and across architectures

Today's processor architectures will not last forever. A setup that runs on an arm64 laptop
and an x86 GPU server now may face a dominant architecture in twenty years that neither was
built for. So preservation is **layered**, from most durable to least, and leans on the
durable layers.

1. **Data in open formats — the most durable layer.** FASTQ (plain text), VCF (plain text),
   CRAM (open published spec), Parquet/DuckDB (open columnar formats). These can be decoded
   from their *public specifications* by software written fresh on any future machine.
   **Archive the specifications themselves** (the SAM/CRAM/VCF spec PDFs, the Parquet spec)
   beside the data. The ultimate guarantee: the data is DNA letters, coordinates and quality
   scores — with the spec in hand, a competent person can recover a genome even if none of
   the original software survives.
2. **Source + pinned manifests — rebuildable on any architecture.** The pipeline is Python
   over open-source tools (htslib, bwa, GATK, VEP) that are portable C/Java/Perl and will
   recompile on, or have successors for, whatever architecture comes next. The git
   repository plus `uv.lock` is this layer. Keep it mirrored (§5).
3. **Container image — a convenience, not the guarantee.** A pinned image (Docker, see
   `docs/docker.md`; or Apptainer/Singularity, which is a single file with no daemon and
   suits archiving) gives bit-for-bit reproducibility for perhaps 10–15 years. But an image
   is **architecture-bound**. Mitigate, don't rely on it:
   - Build **multi-arch** (amd64 + arm64) so it runs on both current worlds.
   - Treat the image as a 15-year convenience; *rebuild from source* is the long-term contract.
   - Emulation (QEMU and similar) is a *possible* fallback for running an old image on a
     future machine — plausible, never guaranteed. Don't make it the plan.
- **Migration over emulation.** Active stewardship beats frozen artifacts. Every decade or
  so, do a **preservation migration**: confirm the data still reads with then-current tools,
  re-pin a fresh multi-arch container, and convert any format that is fading to its
  successor. A living archive that is touched every decade outlasts a sealed time capsule.

---

## 4. Operability handoff (an heir is not necessarily a bioinformatician)

- **Successor runbook** — `docs/runbook.md` is the starting point: what this is, where the
  data lives (and its backups), how to run a scan, **how to add a newly-sequenced relative**
  (ingest → normalize → re-scan the whole family — the recurring generational action), how to
  restore from cold backup, and how to rebuild the environment from source if the container
  won't run. Add your own specifics: where *your* data and backups are.
- **Local UI** — `run.py serve` gives heirs the dashboards without the command line.
- **Scheduled scans** — run `check-releases` on a timer (launchd, cron or a systemd timer)
  so the reclassification engine keeps working untended across custody changes.

---

## 5. Custody & governance (technical)

If the data and the repository live on one machine with one keyholder, that is a single
point of failure for both *access* and *survival*.

- **Decentralize the code:** mirror the git repository to each successor custodian's machine
  and to a remote. Git is already a distributed full-history copy — use that property.
- **Replicate the data store**, not only back it up — at least one successor's location
  should hold a live, independently usable copy.
- **Secrets & access succession:** document where keys, passwords and encryption material
  live and how access transfers on death or incapacity. Fold the genomic archive into the
  actual **estate plan**, with at least two named custodians so it can neither be lost nor
  monopolized.

---

## 6. Consent & ethics charter — *family decisions, fill these in*

A family archive spans people who cannot all consent — minors today, descendants not yet born
— and one person's genome reveals information about every blood relative. That makes a short
written charter worth more than any code. Below are prompts, **not answers**; decide them as
a family and record the decisions here. Check the law where you live: in the EU, for example,
genetic data is a GDPR special category, which shapes lawful storage and access.

- **Who may be sequenced or added?** (Adults by consent; minors — decided by whom, and
  revisited when they reach adulthood?) → _____
- **The never-use list.** Commit in writing to what this data is *never* used for — e.g.
  sharing with insurers, employers, direct-to-consumer sites, or law-enforcement databases.
  → _____
- **Right *not* to know.** Each person chooses whether to learn their own results; children
  opt in as adults rather than having predictive findings imposed. (The pipeline supports
  this per person: `run.py consent`, see `docs/hardware-and-scope.md` §5.) → _____
- **Disclosure of actionable findings.** How is a serious, actionable result (e.g. a BRCA
  variant) surfaced to the person it affects, and at what age? → _____
- **Stewardship roles.** Who are the named custodians, and how is the role passed on? → _____
- **Legal posture.** Lawful basis (consent), each person's rights over their own data within
  the family, retention. → _____

---

## 7. Why preservation pays off

Keeping this alive is not sentimental — the asset gains value over time:

- **Automatic reclassification.** The dated-snapshot + DELTA engine re-interprets the *same*
  genomes against improving databases every cycle. Heirs inherit a process that keeps finding
  new things — including in relatives no longer alive to be re-tested.
- **A deepening pedigree.** Each newly-sequenced relative turns pairs into trios and quads,
  unlocking de-novo detection, phasing and multi-generational segregation of any variant of
  interest. Analytical power grows with the family tree.
- **A private longitudinal record.** Over decades, linking genotypes to who actually developed
  what is insight no consumer test can offer — and it only exists if the data is preserved
  continuously.

---

## Preservation checklist (start here)

- [ ] Verify a real **off-machine** backup of the raw FASTQs exists, with sha256 verified.
- [ ] Add a **third, off-site** copy; enable scheduled **integrity scrubbing**.
- [ ] **Archive the reference FASTA(s)** in the cold copy beside the CRAMs.
- [ ] Build a **multi-arch** container image of the pipeline; keep the build recipe with it.
- [ ] Mirror the **git repository** to at least one successor's machine and a remote.
- [ ] Extend `docs/runbook.md` with your specifics (where the data and backups are).
- [ ] Hold the **family charter** conversation (§6) and record the decisions.
