# AlphaGenome 2 integration — implementation plan

Written 2026-11-11 for a coding agent working in this repository. Read the whole document
before starting; phases are ordered so each one is shippable on its own.

## 1. Goal

Add Google DeepMind's AlphaGenome 2 (released 2026-11-11) as an optional, opt-in source of
non-coding variant-effect evidence, and use the one thing it does that Enformer and Borzoi
cannot: predict a variant's effect **in a named cell type or tissue**.

The work is finished when:

- `INTERPRET_MODEL=alphagenome2 run.py interpret <snapshot>` scores the family's carried rare
  non-coding variants through the AlphaGenome API with no GPU box, and every existing consumer
  of regulatory scores (triage, DELTA, provenance, limitations, dashboards) works unchanged.
- Each panel gene is scored in the tissue relevant to its disease.
- Three reports gain AlphaGenome evidence: a new recessive "second hit" candidate report, the
  de novo report, and a splice column where SpliceAI has no score.
- Nothing leaves the machine unless the operator has configured a key **and** every carrier of
  a variant has consented. The default install sends nothing.
- `uv run python tests/test_units.py` passes offline, with no API key.

## 2. What is known, and how

### Verified by reading the `alphagenome` 2.0.0 wheel from PyPI

These were read from source, not run. Phase 0 confirms them against the live service.

- The client is `alphagenome.models.v2.dna_client`. `create(api_key, model_version=None,
  timeout=None, address=None)` opens a gRPC channel to `gdmscience.googleapis.com:443` and
  authenticates with an `x-goog-api-key` header. There is no REST endpoint in the package.
- `DnaClient.score_variant(interval, variant, condition=None, variant_scorers=(), *,
  organism=HOMO_SAPIENS, cell_conditioning=None, cell_condition_name=None)` returns a list of
  `AnnData`, one per scorer. With no scorers given it uses the recommended set.
- `condition` is a `model_condition.ModelCondition` or a composite CURIE string. A condition
  holds at most one term from each of five ontologies: CLO, UBERON, CL, EFO, NTR.
  `model_condition.from_curies([...])` builds one. One condition per call.
- `DnaClient.supported_conditions()` returns a DataFrame of the reference conditions and, per
  output type, whether that output is observed or imputed for the condition.
- Interval width must be one of 16,384 / 131,072 / 524,288 / 1,048,576 bp.
- Limits in the client: 23 scorers per request, 256 conditioning cells, ISM interval ≤ 10 bp.
- `variant_scorers.RECOMMENDED_VARIANT_SCORERS` keys include `RNA_SEQ` (gene-level log fold
  change), `RNA_SEQ_ACTIVE`, `CAGE`, `DNASE`, `ATAC`, `CHIP_TF`, `CHIP_HISTONE`,
  `SPLICE_SITES`, `SPLICE_SITE_USAGE`, `SPLICE_JUNCTIONS`, `POLYADENYLATION`, `METHYLATION`,
  `CONTACT_MAPS`, and single-cell variants (`SC_ATAC`, `SC_GENE_EXPRESSION`).
- `variant_scorers.tidy_scores(...)` flattens results to a DataFrame with `raw_score` and,
  where the server supplies it, `quantile_score`.
- `ModelVersion` has only `ALL_FOLDS` (default) and `FOLD_1`. **There is no dated or numbered
  server model version to pin.** Section 5.7 deals with that.
- `genome.Variant(chromosome, position, reference_bases, alternate_bases)` is 1-based.
- The package also contains `alphagenome.atlas` (a client for precomputed AlphaGenome 1 scores)
  and `alphagenome.models.v1`. Neither is needed for this plan.
- The package is Apache-2.0 and pulls in pandas, anndata, grpcio, matplotlib, pyarrow and
  scipy. This repository currently depends on only `duckdb` and `requests`.

### From same-day press coverage only (the DeepMind post was not readable when this was written)

Treat these as unconfirmed until Phase 0.

- Free API key for non-commercial research; commercial use is a private preview.
- Model weights are not released. One outlet says they follow "upon final publication".
- Terms: outputs for non-commercial research only, not to be used to train other models, not
  for clinical decisions.
- Rate limits and quotas: no figure found anywhere.

### How the repository handles regulatory scores today

- `pipeline/interpret.py` selects carried rare non-coding variants from DuckDB
  (`_select_variants`), builds ref/alt windows, ships a job over SSH to a GPU box
  (`_run_remote`), and merges results into a genome-level cache (`_merge_scores`).
- The cache is `config.ENFORMER_SCORES_FILE` = `interpret/<model>_scores.tsv`, columns
  `chrom pos ref alt delta_max l2_center top_track top_desc`, keyed by the first four.
- `pipeline/enformer.py` is the read side: `available()`, `lookup()`, `percentile()`,
  `band()`, `label()`, `is_regulatory_hit()`. Bands come from a benign-null calibration
  (`interpret/<model>_calibration.json`, SNV and indel strata) built by `interpret --calibrate`.
- Consumers call only that read API: `triage.model_lookups` and `triage._f_model`,
  `diff.py`, `provenance.py` (`record_history`, `record_provenance`), `limitations.py`.
- `acmg.py` does not use regulatory scores at all. Keep it that way.
- `config.INTERPRET_MODELS` maps a model name to `hf`, `seq_len`, `image`, `targets`.
  `config.set_genomes_root()` re-derives every data path and must be updated with any new one.
- Chromosomes are stored without a `chr` prefix (`"17"`), GRCh38.
- `pipeline/consent.py` gates disclosure per (person, category), default not consented.
- Tests live in one file, `tests/test_units.py`: plain functions using `check(...)`, each
  registered in the list at the bottom of the file.

## 3. Constraints that apply to every phase

1. **Opt-in egress.** No network call to the AlphaGenome service without both an API key and
   consent (section 5.2). The demo (`run.py demo`) never calls it.
2. **Send the minimum.** A request carries a variant, an interval and an ontology term.
   Never put a sample id, person name, family id or carrier list in a request, a request
   metadata header, or a log line that leaves the machine.
3. **Evidence weight stays at triage level.** AlphaGenome scores feed `triage._f_model`
   through the existing band, and appear in reports. They do not become ACMG criteria and
   never set or change an ACMG tier.
4. **Report wording.** Every AlphaGenome-derived statement is a model prediction and a
   candidate, never a finding. Keep the research-use caveat that existing reports carry.
5. **Optional dependency.** `alphagenome` goes in `[project.optional-dependencies]` as an
   extra named `alphagenome`, pinned `>=2.0,<3`, imported lazily inside the functions that
   need it. A base install and the Docker demo must work without it. Update `uv.lock`.
6. **No real API output in git.** Test fixtures use invented numbers. Real scores live under
   `GENOMES_ROOT/interpret/`, which `.gitignore` already covers.
7. **Tests are offline.** All API access goes through one thin wrapper that tests replace
   with a fake. No test may need a key or a network.
8. **Enformer and Borzoi keep working** exactly as now, including their caches and tests.

## 4. Phase 0 — confirm the API before writing pipeline code

This needs an API key and network access to `gdmscience.googleapis.com:443`. If either is
missing, stop and ask the operator; do not guess the answers below.

Write a throwaway script outside the repository tree and record the answers in section 9 of
this file. Use a well-known public variant, not a family variant.

1. Does `score_variant` accept chromosome `"17"` or does it require `"chr17"`? Which assembly
   does it assume, and is there a parameter for it?
2. Print `supported_conditions()`. How many conditions are there, and which CURIEs cover
   liver/hepatocyte, airway epithelium, cardiomyocyte, colon epithelium, breast epithelium,
   neuron, kidney, skeletal muscle, pancreas, retina and whole blood?
3. For one condition, which scorers return rows, and which return `quantile_score`? Is
   `RNA_SEQ` available for every condition or only imputed for some?
4. What does `tidy_scores` return for a gene-level scorer: one row per gene per track? Which
   columns identify the gene (`gene_id`, `gene_name`) and the track?
5. How long does one call take at 1 Mb, and what error comes back when the quota is hit
   (gRPC status code and message)? Is there a documented daily or per-minute limit?
6. Is the same variant scored twice bit-identical? If not, how large is the jitter? Section
   5.7 needs a tolerance.
7. What is the largest indel `score_variant` accepts? Does a 5 kb or 50 kb deletion work as a
   `Variant`, or must it go through `predict_sequence` with a hand-built alternate sequence?
8. Read the current API terms page and confirm the three terms listed in section 2. If the
   terms forbid storing outputs locally or redistributing code that calls the API, stop and
   report.

Decisions that depend on Phase 0: the summary metric in 5.3 (assumes `RNA_SEQ` gene-level
scores exist for most conditions), Phase 5 (assumes splice scorers return usable rows), and
Phase 6 (assumes large deletions can be scored at all).

## 5. Phase 1 — `alphagenome2` as an `interpret` backend

### 5.1 Configuration

- Add an `alphagenome2` entry to `config.INTERPRET_MODELS` and a `transport` key to every
  entry: `"ssh"` for enformer and borzoi, `"api"` for alphagenome2. `seq_len` is 1,048,576.
- New environment settings: `ALPHAGENOME_API_KEY` (read from the environment only, never
  written to disk, logs, reports or provenance), `ALPHAGENOME_MAX_CALLS` (hard cap on API
  calls per run, default 500), `ALPHAGENOME_PANEL_ONLY` (default on: only variants whose gene
  is in `panels/monogenic.tsv`).
- New paths under `INTERPRET_DIR`, also added to `set_genomes_root()`:
  `alphagenome2_scores_long.tsv`, `alphagenome2_egress.tsv`, `alphagenome2_epoch.json`.

### 5.2 Privacy gate

- Add a consent category `external_scoring` in `pipeline/consent.py`, default not consented,
  settable with the existing `run.py consent` command.
- A variant is eligible to be sent only if **every** carrier in the snapshot has consented.
  Relatives share variants, so one person's consent does not cover a variant another person
  also carries.
- Add `interpret --dry-run`: print the count and list of variants that would be sent, the
  conditions, and the number of API calls, then exit without opening a channel.
- Append every request actually sent to `alphagenome2_egress.tsv`: timestamp, chrom, pos,
  ref, alt, condition, scorers. This is the operator's record of what left the machine.
- With no key set, the stage exits with one clear message naming the variable and the docs
  section. With a key but no consenting carriers, it says so and sends nothing.

### 5.3 Scoring

- Split `interpret.interpret_report` so window building and `_run_remote` run only for the
  `ssh` transport. For `api`, skip `samtools`, the FASTA check and SSH entirely.
- New module `pipeline/alphagenome_client.py`, the only place that imports `alphagenome`.
  It exposes roughly `score(variants, condition_for) -> rows` and hides the client, retry and
  tidy logic. Tests substitute this module's client factory.
- For each variant, the interval is the 1 Mb window centred on the variant position. Skip a
  variant whose window runs off the contig, as the existing code does.
- Call sequentially, with exponential backoff on quota and unavailable errors, stopping
  cleanly at `ALPHAGENOME_MAX_CALLS`. A stopped run must still merge what it scored: flush
  results to the cache in batches so an interrupted run loses at most one batch.
- Request an explicit scorer list, not the default set: `RNA_SEQ`, `CAGE`, `DNASE`,
  `CHIP_TF`, `SPLICE_SITES`, `SPLICE_SITE_USAGE`, `SPLICE_JUNCTIONS`, `POLYADENYLATION`.
  Adjust after Phase 0 if some are unavailable.
- Write every returned row to `alphagenome2_scores_long.tsv`: `chrom pos ref alt condition
  scorer gene track raw_score quantile_score epoch`. Key is the first seven columns. This is
  the full record; later phases read it.
- Also write one summary row per variant to the standard `alphagenome2_scores.tsv` in the
  existing eight-column format, so `enformer.lookup()` and every consumer keep working:
  - `delta_max` = absolute `RNA_SEQ` log fold change for the variant's own annotated gene in
    its resolved condition (section 6). If the gene has no row, fall back to the largest
    absolute `DNASE` score at the variant.
  - `l2_center` = 0.0. It is unused by consumers; keep the column for format compatibility.
  - `top_track` = 0, `top_desc` = `"<scorer>:<condition label>"`, e.g. `RNA_SEQ:hepatocyte`.
    `enformer.label()` prints the part before the colon, so the scorer name must come first.
- The summary is a deliberate simplification: one gene, one tissue, one number. Do not take
  a maximum across scorers; their units are not comparable.

### 5.4 Calibration

`interpret --calibrate` must work for the API transport using the same summary metric, with
each null variant scored in the condition its own gene resolves to. Keep the SNV and indel
strata. At the default 400 per stratum this is 800 calls, so make it resumable and have it
respect `ALPHAGENOME_MAX_CALLS` across runs. Null variants are common variants (AF ≥ 0.2) and
are subject to the same consent rule.

### 5.5 Report and CLI text

`interpret._report` hard-codes "run on the GPU box", "predicted regulatory tracks" and "114 kb
receptive field". Make that text per-model. For alphagenome2 the report names the condition
used for each variant, states that scores came from an external API on a given date, and
gives the number of variants withheld for lack of consent. Update the `interpret` help string
in `run.py`.

### 5.6 Provenance

`provenance.record_provenance` already labels the source with `config.INTERPRET_MODEL`. Add
the condition and the epoch id (5.7) to the alphagenome2 provenance value. Add `alphagenome2`
to the source comments in `schema.py`. `_conf_enformer` applies unchanged.

### 5.7 Detecting a silent model change

Because the service exposes no version, keep a canary:

- A fixed list of about five public sentinel variants in a tracked file
  (`panels/alphagenome_sentinels.tsv`), each with a fixed condition.
- At the start of every run, score the sentinels. Compare with the values stored in
  `alphagenome2_epoch.json`, using the tolerance measured in Phase 0.
- If they match, reuse the current epoch id. If they differ, start a new epoch: record the
  date, the installed package version and the new sentinel values, and log a warning.
- Every cache row carries its epoch. Rows from an older epoch stay in the long table but are
  treated as unscored by selection, so they are rescored over subsequent runs.

### 5.8 Tests

- Fake client returning fixed rows: summary row derivation, fallback when the gene is absent,
  long-table merge and re-merge, resume after a simulated quota error mid-run.
- Consent: a variant with one non-consenting carrier is never passed to the client.
- `--dry-run` and the no-key path never construct a client.
- Epoch: unchanged sentinels reuse the epoch; changed sentinels start a new one and
  invalidate old rows for selection.
- The existing `test_enformer` and the report-discovery tests still pass.

## 6. Phase 2 — score each gene in its relevant tissue

- New tracked panel file `panels/gene_tissue.tsv`: `gene`, `condition` (CURIE from
  `supported_conditions()`), `label` (short human name), `rationale`. Cover every gene in
  `panels/monogenic.tsv`. Use the CURIEs found in Phase 0; where a gene's disease tissue has
  no supported condition, pick the closest and say so in `rationale`.
- Resolution rule, in one function used by scoring, calibration and reports: the gene's row
  if present, otherwise a single default condition set in config.
- A gene may list a second condition on another row (for example a cardiac and a skeletal
  muscle condition). Each is a separate API call and a separate long-table row; the summary
  row uses the first listed.
- A test asserts that every panel gene resolves to a condition, and that every CURIE in the
  file parses with `model_condition.from_curies`. Guard that second check so it is skipped
  when the optional package is not installed.
- Per-person and family reports that show a regulatory score now show the tissue beside it.

## 7. Later phases

Each depends on Phase 1 and 2 and on the long table. Do them in this order.

### Phase 3 — second-hit candidates in recessive genes

The clinically useful case: a person carries one pathogenic variant in a recessive gene and
has a compatible phenotype, but no second coding hit.

- New stage and CLI command `second-hit <snapshot>`, new module `pipeline/secondhit.py`,
  report under the family reports directory following `docs/report-formats.md`.
- Index cases: for each sample, each gene marked `AR` in `panels/monogenic.tsv` where the
  sample has exactly one heterozygous variant matching `clinsig`'s actionable expression.
- Candidates: that sample's other heterozygous variants in the same gene with gnomAD AF ≤ 0.01
  or absent, of any consequence class, including intronic. Note that
  `config.INTERPRET_CONSEQUENCES` excludes intronic variants, so this stage needs its own
  selection query. Cap candidates per gene.
- Score candidates in the gene's condition through the Phase 1 client, with the same consent
  rule and call cap. Rank by the expression score and by the splice scorers separately.
- Phase status for each candidate against the known pathogenic variant: trans, cis or
  unknown. Use parental origin as `queries/compound_het.sql` does when parents are loaded.
  `pipeline/phasing.py` currently only considers genes with two or more damaging variants
  (`_candidate_genes`); extend it to accept an explicit pair so read-backed phasing can be
  requested for a candidate. A cis candidate is reported but ranked below unknown.
- Wording: "candidate second variant, model-predicted, unconfirmed". The report must say that
  absence of a candidate does not exclude a second hit.

### Phase 4 — de novo candidates

- Find how `pipeline/denovo.py` persists confirmed candidates and add them as a priority
  class in the interpret selection, ahead of the rare-first ordering, when non-coding.
- Add a regulatory column to the de novo report using `enformer.label()`, with the tissue.
- No change to de novo calling or confirmation logic.

### Phase 5 — splice evidence where SpliceAI is silent

- New read function over the long table returning the splice scorers' values for a variant
  in its gene's condition.
- Where `spliceai.lookup()` returns nothing for a variant with a splice-region or intronic
  consequence, reports show the AlphaGenome splice value in a separate, clearly labelled
  column. It is display-only in this phase: it does not enter `triage._f_model` and does not
  touch `acmg.py`, because there is no validated cut-off for it yet.
- Add the non-validated status to `limitations.py`.

### Phase 6 — structural variants (only if Phase 0 item 7 shows it is possible)

- For deletions and duplications from `pipeline/sv.py` that overlap or lie within 100 kb of a
  gene in `panels/sv_genes.tsv`, and that fit inside one 1 Mb window, score the effect on that
  gene's expression in its condition.
- If only `predict_sequence` works, the alternate sequence must be built locally from the
  reference FASTA. Note that this sends sequence, not coordinates; record it in the egress
  log as such.
- Report as an extra column in the SV report. If Phase 0 shows large events are not
  supported, write that down in section 9 and skip this phase.

### Phase 7 — DELTA on model change

- When a new epoch (5.7) changes a variant's band between two snapshots, DELTA reports it as
  a model-driven change, distinct from a ClinVar or gnomAD change.
- `variant_history.model_band` already stores the band per snapshot; add the epoch id
  alongside it so the two causes can be told apart.
- Model-driven changes rank below database reclassifications and never reach the top triage
  tier on their own.

## 8. Documentation to update as each phase lands

- `README.md`: the `interpret` row in the command table; the "Private by design" bullet must
  state plainly that the alphagenome2 backend is the one feature that sends data out, that it
  is off by default, and what is sent.
- `docs/architecture.md`: the interpret node in both diagrams and the stage table.
- `docs/hardware-and-scope.md`: `interpret` no longer strictly needs a GPU; licence and terms
  summary for the API (non-commercial, each operator needs their own key).
- `docs/docker.md`: `interpret` with the API backend can run in the container if the extra is
  installed and the key is passed; say how.
- `docs/runbook.md`: recipe for getting a key, granting `external_scoring` consent, running
  `--dry-run`, reading the egress log, and what to do when a new epoch is announced.
- `docs/report-formats.md`: the second-hit report and the new columns.
- `docs/stewardship-and-succession.md`: a hosted model can be withdrawn; scores already cached
  remain usable, and Enformer/Borzoi remain the fully local path.

## 9. Phase 0 findings

Fill this in before starting Phase 1.

Measured 2026-11-12 against the live service with `alphagenome` 2.0.0, using public
variants only (the scoring-tutorial variant chr2:219243701 C>T, and the five sentinels).

| Question | Answer |
|---|---|
| Chromosome naming and assembly | `chr` prefix required: bare `"2"` fails with `INVALID_ARGUMENT: Unrecognized chromosome '2' for organism`. GRCh38; no assembly parameter. The client adds `chr` and strips it on the way back. |
| Conditions available for the panel's tissues | `supported_conditions()` returns 2,799 conditions. RNA-seq is `observed` for 306, `imputed_supported` for 2,183, `imputed_exploratory` for 310. Every panel gene is mapped in `panels/gene_tissue.tsv` (85 CURIEs, all present in the list and all parse with `from_curies`); the rationale column records the RNA status and any proxy choice. PMP22's Schwann-cell condition is `imputed_exploratory`. |
| Scorers that return rows, and which have quantiles | `RNA_SEQ`, `CAGE`, `DNASE`, `CHIP_TF`, `SPLICE_SITES`, `SPLICE_SITE_USAGE`, `SPLICE_JUNCTIONS` return tidy rows. `POLYADENYLATION` returns an AnnData but no tidy rows, so it is not requested. No scorer returns `quantile_score`: v2 gives `raw_score` only (quantiles exist in AlphaGenome 1 only), which is why bands need `interpret --calibrate`. |
| Shape of gene-level tidy output | One row per gene per track. Gene: `gene_id` (Ensembl), `gene_name`; track: `track_name`, plus `imputation_mode`, `biosample_curie`, `biosample_name`. `RNA_SEQ` returned 73 genes × 5 tracks for one variant in one condition; splice junction rows add `junction_Start`/`junction_End`. |
| Latency per call, quota, quota error | ~1.5 s per 1 Mb call with 8 scorers. No quota was hit. No published limit: DeepMind's community answer says quotas are heuristic and change with demand. The client backs off on `RESOURCE_EXHAUSTED`, `UNAVAILABLE` and `DEADLINE_EXCEEDED`, then stops cleanly and keeps what it scored. |
| Determinism and tolerance | Two identical calls were bit-identical (max difference 0). The epoch canary uses rel 1e-6 / abs 1e-7. |
| Largest scorable indel | 5 kb and 50 kb deletions and duplications are accepted by `score_variant` when passed as full reference/alternate bases (~1.6–2.0 s each). Larger events were not tried; Phase 6 caps at 50 kb. This sends sequence, not just coordinates. |
| Terms confirmed | Terms and Output Terms read 2026-11-12: outputs are for non-commercial use, must not train models that predict variant effects, and must not be used for clinical purposes or relied on as medical advice. Storing outputs locally is not forbidden. Publishing or distributing outputs requires conspicuous notice that they are subject to the Output Terms, so every AlphaGenome-bearing report carries that link. |

### Implementation notes (deviations from the plan)

- The long-table key also includes `gene_id`, `junction_start` and `junction_end`.
  Without them, splice-junction rows for one gene and track would overwrite each other.
- The summary metric is the largest absolute own-gene `RNA_SEQ` value across that
  condition's RNA tracks (one gene, one tissue, one scorer), with the `DNASE` fallback.
- Uncalibrated alphagenome2 scores band as `low` and add no triage weight. Raw scores have
  no validated cut-off, so bands appear only after `interpret --calibrate` for the current
  epoch. A calibration from an older epoch is ignored.
- Phase 6 consent: SVs are not in DuckDB, so carriers among relatives cannot be listed.
  SV scoring therefore needs `external_scoring` from every sequenced person in
  `samples.tsv`, which fails closed.
- Phase 3 phenotype compatibility is not recorded anywhere in the pipeline. The second-hit
  report says so and leaves it to human review.

## 10. Out of scope

- Running AlphaGenome 2 locally. The weights are not available.
- A local AlphaGenome 1 backend. It is possible in principle with the published v1 weights,
  but it has no cell-type conditioning and is a separate piece of work.
- The AlphaGenome Atlas of precomputed v1 scores.
- Custom single-cell conditioning (`cell_conditioning=`). The family has no single-cell data.
- Any change to ACMG scoring, polygenic scores, pharmacogenomics or HLA.
