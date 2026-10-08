# Operator Runbook

Plain-language operations guide for whoever runs this next. You do **not** need to be
a bioinformatician — follow the recipes. For *why* it's built this way, see
`docs/architecture.md` and `docs/stewardship-and-succession.md`.

## What this is, in one paragraph

A re-runnable pipeline that turns a family member's raw DNA reads into a health report,
and — the key point — **re-checks everyone's DNA against updated medical databases on a
schedule**, telling you what newly became significant since last time. The raw DNA is
never changed; reports are regenerated.

## Where everything lives

- **Code** (this repo): `~/Projects/genome-ledger` (or cloned from your git remote). Safe to clone anywhere.
- **Data** (never in git): `~/genomes` (the `GENOMES_ROOT`). Holds raw reads, reference,
  databases, and reports. **This is the irreplaceable part — see backups below.**
- Reports land in `~/genomes/reports/` — markdown under `md/`, the pages you open under `html/`.

## Daily reality: you don't do anything

A weekly `launchd` job runs `check-releases`; it only scans when a database actually
changed (or when the last scan did not finish, or left out a sample that has since been
normalized), then writes new reports. The same job re-scores the polygenic reports when the
trait panel has changed since they were written — add a score to
`panels/complex_traits.tsv` and everyone already scored is brought up to date on the next
run (`docs/prs-guide.md`, "Adding a score"). Open `~/genomes/reports/html/index.html`.
The DELTA report groups what changed into
**Actionable now / Monitor closely / Log only** — read the Actionable tier first.

Scans are **incremental**: unchanged samples are copied forward and, when only ClinVar
moved, only the affected variants are re-annotated — so a routine weekly scan is fast.
Once a month (and on any tool-version bump) it does a full re-annotation automatically.

To run a scan by hand:
```bash
cd ~/Projects/genome-ledger
uv run python run.py scan            # incremental (normal)
uv run python run.py scan --full     # force a full re-annotation of everyone
```

If you ever doubt the incremental path, `scan --full` reproduces a snapshot from scratch;
`uv run python run.py validate --replay <old-date> --truth <newer-date>` checks how well
past alerts held up in later ClinVar.

## Recipe: add a newly-sequenced family member

This is the main thing you'll ever do. Say the new person is "kid2", sequenced as files
`ACME123_1.fq.gz` / `ACME123_2.fq.gz`. Other common namings work too: `ACME123_R1.fastq.gz` /
`_R2`, or Illumina's one pair per lane (`ACME123_S1_L001_R1_001.fastq.gz`, `_L002_`, …), which
`align` streams in lane order. It reads the sequencer (MGI/DNBSEQ or Illumina) off the read
names; set `READ_PLATFORM` only if it guesses wrong (e.g. `ELEMENT` for AVITI reads).

1. Put the read files in `~/genomes/incoming/`.
2. Add a line to `~/genomes/incoming/manifest.tsv` (tab-separated):
   ```
   kid2	ACME123	son	male
   ```
3. Run the front-end, then fold them into the family:
   ```bash
   cd ~/Projects/genome-ledger
   uv run python run.py align kid2     # reads → CRAM   (hours)
   uv run python run.py call  kid2     # CRAM → VCF     (hours)
   uv run python run.py ingest ~/genomes/called/kid2.GRCh38.vcf.gz \
       --sample kid2 --relation son --sex male
   uv run python run.py normalize kid2
   uv run python run.py scan           # re-annotates the WHOLE family together
   ```
4. Add kid2's row to `~/genomes/family/family.tsv` with his `father` and `mother` ids
   (format: `family.example.tsv`). That file is the family graph; `pedigree.ped` is
   generated from it, so don't edit the `.ped` by hand.

Each step is **idempotent and resumable** — safe to re-run; finished steps are skipped.
The long jobs are best run via the unattended driver pattern (see "Running unattended").

## Recipe: read the reports

```bash
uv run python run.py serve      # then open http://127.0.0.1:8765/
```

or open `~/genomes/reports/html/index.html` straight from disk. Pick a person; their page
opens on a summary and the tree on the left lists every report. The **Snapshots** section
at the foot of the tree switches between database releases, **Viewing as** relabels the
family from one person's point of view, and Light / Dark / Auto sets the theme. If you want
the reports as files rather than pages, they are markdown under `reports/md/` and JSON
under `reports/json/` (`docs/report-formats.md`).

## Recipe: I already have a VCF (vendor, or the GPU box)

Skip align/call — start at `ingest` with the VCF, then `normalize`, then `scan`.
For GPU-accelerated calling on a separate box, see `docs/hardware-and-scope.md`.

If the VCF is GRCh37/hg19, `ingest` will refuse it (it checks the chr1 contig length).
Lift it to GRCh38 first — see the LiftoverVcf recipe in `README.md`. Sanity-check the
result before trusting it: `bcftools norm` reporting `mismatch_removed: 0` against the
GRCh38 FASTA means every lifted REF agrees with the reference.

**chrM does not survive a build liftover.** hg19 and GRCh38 use different mitochondrial
reference sequences, and the UCSC chain maps chrM naively — every MT variant comes across
offset by 1-2 bp, and `RECOVER_SWAPPED_REF_ALT` then swaps REF/ALT at the wrong base. The
result looks like valid data: haplogrep returns the tree-root haplogroup (H2a2a1) at low
quality, which will contradict the maternal line. `ancestry` now refuses any MT call
below `MT_MIN_QUALITY`, but the honest fix is to exclude chrM from a lifted callset and call
mtDNA from the reads (`mito`) instead. The check that catches this in one line: a child's MT
haplogroup **must** equal their mother's.

chrY, by contrast, lifts cleanly: a son's Y haplogroup off a lifted callset matches his
father's. So the problem is chrM specifically, not liftover in general. The two haplogroup checks together (Y must match the father, MT must match the
mother) are the cheapest parentage/sanity test the pipeline has; run them before trusting
anything else about a newly-ingested family member.

First run of `ancestry` downloads Yleaf's own ~938 MB hg38 reference before it can predict
anything, and will appear to hang or fail on a slow link. It is cached afterwards — if the
stage dies the first time, just re-run it.

A VCF-only sample is a *thinner* sample. `traits`, `hla`, `prs`, `ancestry` and all the
DuckDB reports work off its normalized VCF; `sv`, `repeats`, `callable`, `mito`, `phase`,
`hla-type`, `force-call` and CYP2D6-within-`pgx` need reads and cannot run. If the lab
will give you a BAM/CRAM, take it over a FASTQ — it is half the storage and
`samtools collate | samtools fastq` recovers reads for realignment.

## Recipe: force-calling — making "not in the VCF" mean something

**Why.** A VCF lists only the places where someone *differs* from the reference. So when a
report looks up a fixed position — a trait SNP, an HLA tag, one of the variants in a small
polygenic score — and finds nothing, there are two possibilities: the person matches the
reference there, or the sequencer never got a usable look at that spot. Without
force-calling the pipeline assumes the first. That is right most of the time and silently
wrong the rest: a low-coverage position is scored as "reference", and a no-call at a celiac
tag SNP reads as "haplotype absent".

`force-call` goes back to the reads and genotypes exactly those positions, whether or not
they vary. An *assumed* reference genotype becomes an *evidenced* one, and a position the
reads cannot resolve becomes an explicit no-call that the reports show as indeterminate
instead of as a negative.

**When to do it.**

- Once per person, after `align` (it needs their CRAM), before you rely on `traits`, `hla`
  or the small polygenic scores. It is optional: everything runs without it, with the
  assumption above stated in each report's caveats.
- Again when the list of positions changes — you added a trait, an HLA tag or a small
  score to the panel. The stage notices by itself: it re-calls when the site list or the
  CRAM differs from last time, and skips otherwise.
- Not for a VCF-only person. There are no reads to go back to.
- Not for genome-wide polygenic scores (more than 50,000 variants). Those are scored with
  `plink2`, which handles missing reference positions its own way.

**How.**

```bash
uv run python run.py force-call <sample>             # reads the CRAM; minutes, not hours
uv run python run.py traits <sample> --force         # re-run what uses it
uv run python run.py hla    <sample> --force
uv run python run.py prs    <sample> --force
```

The first run also builds the list of positions, which downloads the small scoring files
from the PGS Catalog; that list is shared by everyone and cached. The consumers need
`--force` only because their existing report would otherwise be reused.

**How to tell it worked.**

- A **Force-call QC** report appears on the person's page: how many positions were asked
  for, how many were called, and the depth. Expect essentially all of them called on a 30×
  genome; a noticeable share of no-calls means a coverage problem worth knowing about.
- The caveats in the `traits`, `hla` and `prs` reports change from "absent positions are
  taken as reference" to "read from the force-called VCF".
- A result that was "negative" may become "indeterminate". That is the stage doing its job.

**What it does not do.** It does not find new variants, and it does not change the FULL,
panel or ACMG reports, which come from the ordinary calls. The trio stage `denovo` uses the
same machinery on the parents for a different purpose and runs it for you.

## Recipe: trio de-novo variants

Needs a sequenced child plus **both parents' CRAMs** (the confirmation reads the parents'
reads, not their VCFs). Set the child's `father`/`mother` in `~/genomes/family/family.tsv`,
then:

```bash
uv run python run.py denovo <child>              # full: candidates + parental confirmation
uv run python run.py denovo <child> --no-confirm # candidates only (FP-dominated; for tuning)
```

Read the headline before the table. The stage compares its own output against the published
germline rate (40–150 de-novo SNVs per generation) and says plainly when the count is too
high to be a call set. If the child is VCF-only, expect it to be — child-side artefacts
cannot be distinguished from real de-novo events without the child's reads.

## Recipe: compound heterozygosity in a recessive gene

Two damaging variants in one recessive gene mean very different things depending on whether
they sit on the same chromosome copy. Two stages answer it, and they are complementary:

```bash
uv run python run.py phase <sample>          # from reads: any sample, but only ~1 kb apart
uv run python run.py segregation <child>     # from transmission: any distance, needs a trio
```

`phase` (WhatsHap) walks the reads, so it works for anyone with a CRAM but goes quiet when the
variants are further apart than a read-pair spans. `segregation` uses the fact that the child's
two copies came from different parents, which has no distance limit and needs no reads from the
child at all — only the parents' CRAMs, to confirm that a parent genuinely lacks a variant
rather than merely having no record of it. Where both stages have an opinion they should agree.

Verdicts: **in-trans** (one from each parent — both copies hit, the reviewable case),
**in-cis** (both from one parent — the other copy is intact), **unresolved**. Note that when
*both* parents carry a variant, transmission cannot tell the copies apart; that is a genuine
limitation, not a failure, and reads are the only recourse.

## Recipe: give a file to a doctor or researcher

The working VCFs and CRAMs record how they were made: every command, with paths that include
your home folder, the run dates, and internal or provider sample ids. Don't send those. Make a
copy for sharing instead:

```bash
uv run python run.py share kid2            # → ~/genomes/share/kid2/kid2.GRCh38.vcf.gz (+ .tbi)
uv run python run.py share kid2 --cram     # → ~/genomes/share/kid2/kid2.cram (+ .crai)
```

The copy's header is rebuilt from an allowlist and the sample is named after the person
(`--name` to choose another). What it cannot change: the file *is* that person's genome, and
the reads in a CRAM still carry the sequencer's run and flowcell ids in their names. A CRAM
also needs the same GRCh38 reference to be read — say which one when you send it.

## Running unattended (long jobs)

Calling a whole genome takes hours. Run it **detached** so it survives the terminal
closing, and keep the Mac awake:
```bash
nohup caffeinate -i -s bash your_steps.sh > ~/genomes/logs/run.log 2>&1 &
```
Set `NOTIFY_CMD` (see the README) to be told when a scan finishes.
⚠️ Closing the laptop lid still sleeps the machine (pauses, resumes on wake).

## Restore from backup / move to a new machine

1. Install prerequisites and clone the repo:
   ```bash
   git clone <your-repo-url>/genome-ledger.git && cd genome-ledger
   uv sync
   ```
2. Restore `~/genomes` from backup (NAS / off-site copy). **The raw `incoming/`
   FASTQs are the master copy** — everything else regenerates from them.
3. `uv run python run.py setup` rebuilds tools, the reference + bwa index, and the VEP
   cache. Then re-run `scan`.
4. If reproducing exactly, use the pinned container (see stewardship doc) rather than
   `setup`, since external download URLs rot over time.

## Backups — the one thing that matters

The raw reads of someone who has died can never be re-made. Keep **3 copies, 2 media,
1 off-site**, and verify checksums periodically (see `stewardship-and-succession.md`).
Archive the reference FASTA *with* the data — a CRAM is unreadable without it.

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `call` fails instantly, "CRAM version 3.1 is not supported" | samtools wrote a too-new CRAM for GATK. Pinned to 3.0 in `call.py`; if it recurs, transcode: `samtools view --output-fmt cram,version=3.0`. |
| `normalize` output is empty | the VCF had no `FILTER=PASS` records (raw GATK output is unfiltered). `call` hard-filters to add PASS; check that ran. |
| A scan stage failed and stopped | read that stage's log under `~/genomes/logs/`. Fix, then re-run the same command — it resumes. The retry is compared against the last scan that *completed*, so nothing that changed in between is lost. |
| The scan notification says a sample was "NOT scanned" | it has no normalized VCF. Run `normalize <sample>`; the next `check-releases` scans it without waiting for a new ClinVar release. A scan in which no sample is normalized fails outright. |
| A report says results were withheld | incidental findings are opt-in per person: `uv run python run.py consent <sample> adult_onset_untreatable yes`, then re-run that stage (`repeats`, `prs`) — it sees the change and regenerates without `--force`. |
| Someone withdraws that consent | `uv run python run.py consent <sample> adult_onset_untreatable no` deletes their repeat-expansion and polygenic reports in every format and rebuilds their dashboard page. Re-run the stages to get the withheld versions back. |
| Dashboards don't show reports that exist | the reports are in the old layout — see "Moving from the old layout" in `docs/report-formats.md`. |
| Pages look unstyled or use the wrong fonts | run `uv run python run.py render`; it rebuilds `reports/html/` including `assets/fonts/`. |
| Disk nearly full | local Time Machine snapshots may be hiding freed space: `tmutil deletelocalsnapshots /`. |
| Tests | `uv run python tests/test_units.py` — sanity-checks inheritance, scatter, triage, validation, and incremental-scan logic. |

## Command reference

See the README's command table, or `uv run python run.py --help`.
