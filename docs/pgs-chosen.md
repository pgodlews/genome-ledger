# Polygenic scores you choose

The [curated panel](prs-guide.md) scores one fixed score per trait, picked by
[objective criteria](polygenic-panel-design.md). This is the other way in: search the whole
[PGS Catalog](https://www.pgscatalog.org/) (about 7,000 scores), pick any score, decide per
person whether to compute it, and read the results in a report of their own.

It uses the same scoring engine as the panel (`pipeline/prs.py`). What it adds is the
search, the consent matrix, bundles for air-gapped machines, and a separate report.

## The short version

```bash
uv run python run.py pgs sync                       # once, on a connected machine
uv run python run.py pgs search "macular degeneration"
uv run python run.py pgs show PGS004606
uv run python run.py pgs consent PGS004606          # the family matrix; tick people
uv run python run.py pgs score PGS004606            # scores whoever consented
```

Or in a browser: `uv run python run.py serve`, then open `http://127.0.0.1:8765/pgs/`.

## Consent comes before any computing

Each person decides per **trait**. A score is downloaded and computed for someone only after
they have consented to every trait it maps to. Nothing is computed and then hidden, so a
result nobody agreed to never exists on disk. That keeps it out of backups, dashboards and
anything that reads the reports.

- Grants go in `family/consent.tsv` as category `pgs:<trait id>`, with all the other consent
  grants. A grant belongs to the person, so it covers all of their callsets.
- Consent is per trait and results are per score. Agreeing to a trait covers a newer score
  for it too, and every result still names the exact score that produced it.
- **Withdrawing consent deletes the data.** `pgs consent PGS… --revoke NAME`, an unticked box
  in the picker, or `run.py consent NAME pgs:<trait id> no` deletes every result for that
  trait, for every callset of the person, and rebuilds the report from what remains.
- The consent matrix shows this warning every time: medical genetics guidelines advise against
  predictive testing of children for adult-onset conditions with no treatment. There is no age
  check. That is advice, and the decision is the parent's.
- A trait the curated panel treats as severe with no established intervention (its S1 list)
  is marked ⚠ in the matrix. It still needs the same explicit consent as any other.

## Picking a score

`pgs search TERM` lists matching traits; with one match (or `--trait EFO_…`) it lists the
trait's scores. `pgs categories` lists the Catalog's own trait categories. `pgs show PGS…`
prints everything about one score, including every published evaluation.

The **published performance** column is deliberately conservative:

- It is the **median** over European evaluation cohorts, not the best one. Cohorts with an
  impossible case count, or fewer than 100 cases, are left out (they are still listed by
  `show`).
- Only evaluations **without age, sex or clinical covariates** count as the score's own
  figure. Ancestry principal components are allowed. A model that includes age reports an
  AUROC that age earns as much as the score does.
- When the only figures come from models with covariates, it says
  `only with age/sex or other covariates`. A covariate field recorded as `0` or `NA` says
  `covariates not stated`.

Some scores cannot be computed by an additive engine. These are weights published as plain
odds or hazard ratios, per-dosage or haplotype weight columns, and interaction terms. They
are marked `cannot score` and refused with the reason, instead of producing a wrong number.

## What it costs

| Score size | First time | Each further person |
|---|---|---|
| Small (≤ 50,000 variants) | seconds with a local 1000 Genomes copy (`run.py mirror-1kg`), minutes without | seconds |
| Genome-wide | 10–15 min on its own; far less per score in a batch (below) | minutes |

The reference work does not depend on the family. It is done once per score and reused.

Almost all of a genome-wide score's reference build is reading each 1000 Genomes chromosome
end to end, whatever positions are asked for. `pgs fetch --batch N` (and the panel's own `prs`
stage, which batches automatically) reads each chromosome **once per batch** and then scores
every member against that one panel. The result is the same as building each score alone; a
test builds both ways and compares them.

Batches are sized to the machine. A batch closes at `PRS_REF_BATCH` scores (default 40), or
earlier before its variant rows or distinct positions pass what a quarter of the RAM and half
the free disk allow (`PRS_BATCH_MAX_ROWS`, `PRS_BATCH_MAX_UNION` cap them further). Scores
built on the same sites, typically one publication's, share almost all their positions, so
listing ids in publication order makes batches cheap. Each parallel plink2 run gets its share
of the cores and memory rather than all of them. `PRS_REF_BATCH=1` turns batching off.

## Air-gapped machines: bundles

A machine with genomes on it never needs the network. Everything a score needs that does not
depend on the family is built elsewhere and carried in.

```bash
# On a connected machine. It needs no genomes: only public files are downloaded.
uv run python run.py pgs sync --to /media/usb
uv run python run.py pgs fetch PGS004606 PGS000040 --to /media/usb --batch 40 --prune

# On the air-gapped machine
uv run python run.py pgs import /media/usb
PGS_OFFLINE=1 uv run python run.py pgs score PGS004606
```

- The browsing and consent steps work offline too, against the carried catalog snapshot
  (about 12 MB). The only thing that needs to cross *outward* is the list of PGS ids you
  want. The picker's *Download wishlist* button writes it; `pgs fetch --from-file` reads it.
- A bundle holds the scoring file, the 1000 Genomes reference (frequencies for a small score;
  the reference distribution and the REF alleles at its sites for a genome-wide one) and a
  sha256 manifest. `import` checks every file and refuses a damaged bundle whole. A
  genome-wide bundle is tens to ~170 MB.
- `--prune` deletes the reference panel's genotype intermediates once a genome-wide score's
  bundle data exists. They can run to gigabytes per score and are only read to build it. Run
  `pgs fetch <ids> --prune` over scores already cached by the panel to reclaim that space.
- For a long list, `--shard I/N` takes every N-th id starting at I, so N machines can each
  build a slice. `PGS_OFFLINE=1` makes `pgs score` refuse to reach the network.
- Keep bundle stores private. Some scoring files' licences forbid redistribution.

## Where results go

| What | Where |
|---|---|
| one result per callset and score | `prs_chosen/<callset>/<PGS id>.json` |
| one report per callset | `reports/md/genome/<person>/chosen/prs_<callset>.md` (and the same path under `reports/html/`) |
| catalog snapshot | `tools/pgs_catalog/<release>/` |
| scoring files and reference data | `tools/pgs/`, shared with the panel |

The report is listed on the person's dashboard as *Polygenic scores you chose*. The curated
panel report is never changed by any of this.

## Keeping results current

Every result records what it was computed from: the scoring file, its harmonisation date, the
reference data, the genotypes and the engine version.

- `pgs status` lists every result as `current`, or names what changed.
- `pgs score --stale` recomputes the stale ones, for people whose consent still stands.
- `pgs check` asks the Catalog whether a scoring file was **re-harmonised**: the same model,
  with its positions re-mapped (`#HmPOS_date` in its header). If so, run
  `pgs fetch <id> --refresh`, then `pgs score --stale`. It also lists **newer scores** for
  the same trait. Those are not run automatically: switching models is a new choice,
  although consent for the trait already covers it.
- A PGS id is never reused for a different model. When authors revise a score, it gets a new
  id.

## A bundle store that follows the Catalog

The PGS Catalog publishes a release about once a month: usually 5–65 new scores, occasionally a
dump of a thousand or more. `pgs update --store DIR` keeps a store of bundles in step:

```bash
uv run python run.py pgs update --store /path/to/pgs-bundles      # e.g. from a nightly timer
```

- Each run costs one request to the Catalog. Nothing else happens until there is a new release.
- On a new release it takes a fresh catalog snapshot, builds bundles (in batches) for every
  computable score the store does not hold into `DIR/<release>/`, and remembers the scores the
  engine cannot model in `DIR/refused.tsv` so they are not retried. More than `--max-new` (300)
  new scores and it builds nothing, and says so: that is a job for a bigger machine.
- For scores the family was actually scored with, it checks whether the Catalog re-harmonised
  their files, rebuilds those bundles and recomputes the results that depended on them.
- It sends one notification when something changed. With `NOTIFY_QUIET_HOURS=00:00-07:00`
  (local time in `NOTIFY_TZ`, default Europe/London) a message raised in that window is held and
  goes out after it: with the next notification, or when `run.py notify-flush` runs from a
  morning timer.

Every release folder in the store is `pgs import`-able on its own.
