# Polygenic risk scores — what is included, and how to add one

A practical companion to two longer documents:
[`polygenic-panel-design.md`](polygenic-panel-design.md) explains *why* these traits were
chosen, and [`prs-architecture.md`](prs-architecture.md) explains *how* scoring works.
This page answers two questions: which scores does `run.py prs` report, and what do I do to
add another.

## The panel

`run.py prs <sample>` scores every row of
[`panels/complex_traits.tsv`](../panels/complex_traits.tsv): 44 traits, one pinned
[PGS Catalog](https://www.pgscatalog.org/) score each. 22 are small scores handled
by the exact per-variant engine; the other 22 are genome-wide and go through
`plink2`. The table below mirrors that file (a test fails if the two drift apart).

| Trait | PGS ID | Variants | Engine | Score AUROC/C | Actionability | Consent-gated |
|-------|--------|---------:|--------|--------------:|---------------|:-------------:|
| angina pectoris | [PGS012540](https://www.pgscatalog.org/score/PGS012540/) | 969,769 | plink2 | 0.64 (AUROC) | Preventive treatment |  |
| asthma | [PGS002727](https://www.pgscatalog.org/score/PGS002727/) | 985,837 | plink2 | 0.66 (AUROC) | Preventive treatment |  |
| coronary artery disorder | [PGS003446](https://www.pgscatalog.org/score/PGS003446/) | 538,084 | plink2 | 0.82 (AUROC) | Preventive treatment |  |
| gout | [PGS004931](https://www.pgscatalog.org/score/PGS004931/) | 1,138 | exact | 0.73 (AUROC) | Preventive treatment |  |
| heart failure | [PGS012544](https://www.pgscatalog.org/score/PGS012544/) | 969,507 | plink2 | 0.62 (AUROC) | Preventive treatment |  |
| hypertensive disorder | [PGS004934](https://www.pgscatalog.org/score/PGS004934/) | 15,056 | exact | 0.72 (AUROC) | Preventive treatment |  |
| Ischemic stroke | [PGS012545](https://www.pgscatalog.org/score/PGS012545/) | 969,253 | plink2 | 0.6 (AUROC) | Preventive treatment |  |
| low density lipoprotein cholesterol measurement | [PGS000875](https://www.pgscatalog.org/score/PGS000875/) | 36 | exact | 0.67 (AUROC) | Preventive treatment |  |
| peripheral arterial disease | [PGS012546](https://www.pgscatalog.org/score/PGS012546/) | 969,771 | plink2 | 0.62 (AUROC) | Preventive treatment |  |
| stroke disorder | [PGS000039](https://www.pgscatalog.org/score/PGS000039/) | 3,225,583 | plink2 | 0.582 (AUROC) | Preventive treatment |  |
| systolic blood pressure | [PGS002807](https://www.pgscatalog.org/score/PGS002807/) | 1,084,157 | plink2 | 0.69 (C-index) | Preventive treatment |  |
| triglyceride measurement | [PGS004937](https://www.pgscatalog.org/score/PGS004937/) | 8,337 | exact | 0.68 (AUROC) | Preventive treatment |  |
| urate measurement | [PGS000126](https://www.pgscatalog.org/score/PGS000126/) | 114 | exact | 0.67 (AUROC) | Preventive treatment |  |
| venous thromboembolism | [PGS000043](https://www.pgscatalog.org/score/PGS000043/) | 297 | exact | 0.67 (C-index) | Preventive treatment |  |
| vitamin D level | [PGS000882](https://www.pgscatalog.org/score/PGS000882/) | 1,094,650 | plink2 | 0.63 (AUROC) | Preventive treatment |  |
| abdominal aortic aneurysm | [PGS003429](https://www.pgscatalog.org/score/PGS003429/) | 831,447 | plink2 | 0.708 (AUROC) | Screening pathway |  |
| atrial fibrillation | [PGS004929](https://www.pgscatalog.org/score/PGS004929/) | 178,814 | plink2 | 0.71 (AUROC) | Screening pathway |  |
| basal cell carcinoma | [PGS000730](https://www.pgscatalog.org/score/PGS000730/) | 47 | exact | 0.624 (AUROC) | Screening pathway |  |
| brain aneurysm | [PGS003408](https://www.pgscatalog.org/score/PGS003408/) | 6,671,269 | plink2 | 0.77 (C-index) | Screening pathway |  |
| breast carcinoma | [PGS000004](https://www.pgscatalog.org/score/PGS000004/) | 313 | exact | 0.673 (C-index) | Screening pathway |  |
| celiac disease | [PGS000040](https://www.pgscatalog.org/score/PGS000040/) | 228 | exact | 0.9 (AUROC) | Screening pathway |  |
| chronic kidney disease | [PGS012543](https://www.pgscatalog.org/score/PGS012543/) | 969,768 | plink2 | 0.58 (AUROC) | Screening pathway |  |
| colorectal cancer | [PGS000149](https://www.pgscatalog.org/score/PGS000149/) | 41 | exact | 0.57 (C-index) | Screening pathway |  |
| cutaneous melanoma | [PGS019946](https://www.pgscatalog.org/score/PGS019946/) | 7,280,842 | plink2 | 0.643 (C-index) | Screening pathway |  |
| endometrial carcinoma | [PGS002735](https://www.pgscatalog.org/score/PGS002735/) | 19 | exact | 0.56 (AUROC) | Screening pathway |  |
| endometriosis | [PGS003447](https://www.pgscatalog.org/score/PGS003447/) | 14 | exact | 0.6 (AUROC) | Screening pathway |  |
| glaucoma | [PGS000137](https://www.pgscatalog.org/score/PGS000137/) | 2,673 | exact | 0.71 (AUROC) | Screening pathway |  |
| hypothyroidism | [PGS004935](https://www.pgscatalog.org/score/PGS004935/) | 6,127 | exact | 0.7 (AUROC) | Screening pathway |  |
| inflammatory bowel disease | [PGS000017](https://www.pgscatalog.org/score/PGS000017/) | 6,907,112 | plink2 | 0.69 (AUROC) | Screening pathway |  |
| juvenile idiopathic arthritis | [PGS000114](https://www.pgscatalog.org/score/PGS000114/) | 26 | exact | 0.657 (AUROC) | Screening pathway |  |
| keratoconus | [PGS002292](https://www.pgscatalog.org/score/PGS002292/) | 36 | exact | 0.756 (AUROC) | Screening pathway |  |
| lung carcinoma | [PGS004691](https://www.pgscatalog.org/score/PGS004691/) | 950,353 | plink2 | 0.618 (AUROC) | Screening pathway | yes |
| ovarian carcinoma | [PGS002250](https://www.pgscatalog.org/score/PGS002250/) | 27,240 | exact | 0.624 (AUROC) | Screening pathway |  |
| prostate carcinoma | [PGS003415](https://www.pgscatalog.org/score/PGS003415/) | 268 | exact | 0.79 (C-index) | Screening pathway |  |
| systemic lupus erythematosus | [PGS000328](https://www.pgscatalog.org/score/PGS000328/) | 57 | exact | 0.83 (AUROC) | Screening pathway |  |
| type 1 diabetes mellitus | [PGS000024](https://www.pgscatalog.org/score/PGS000024/) | 85 | exact | 0.921 (AUROC) | Screening pathway |  |
| ulcerative colitis | [PGS002066](https://www.pgscatalog.org/score/PGS002066/) | 566,637 | plink2 | 0.57 (AUROC) | Screening pathway |  |
| body mass index | [PGS000716](https://www.pgscatalog.org/score/PGS000716/) | 295 | exact | 0.739 (AUROC) | Lifestyle |  |
| high density lipoprotein cholesterol measurement | [PGS004932](https://www.pgscatalog.org/score/PGS004932/) | 17,151 | exact | 0.6 (AUROC) | Lifestyle |  |
| type 2 diabetes mellitus | [PGS003103](https://www.pgscatalog.org/score/PGS003103/) | 945,820 | plink2 | 0.75 (AUROC) | Lifestyle |  |
| exocrine pancreatic carcinoma | [PGS004693](https://www.pgscatalog.org/score/PGS004693/) | 6,351,686 | plink2 | 0.543 (AUROC) | No intervention | yes |
| multiple sclerosis | [PGS002726](https://www.pgscatalog.org/score/PGS002726/) | 476,399 | plink2 | 0.73 (AUROC) | No intervention | yes |
| renal carcinoma | [PGS004690](https://www.pgscatalog.org/score/PGS004690/) | 6,351,669 | plink2 | 0.564 (AUROC) | No intervention | yes |
| urinary bladder carcinoma | [PGS004687](https://www.pgscatalog.org/score/PGS004687/) | 1,077,775 | plink2 | 0.558 (AUROC) | No intervention | yes |

- **Engine** follows from the variant count: above `PRS_PLINK_MIN` (50,000) a score is
  genome-wide and scored with `plink2`; at or below it, by the exact per-variant path.
- **Score AUROC/C** is the published discrimination of the score on its own (no age, sex or
  other covariates) in a European evaluation cohort of at least 1,000 people. It bounds how
  much any percentile can mean.
- **Actionability** decides where a trait appears in the report, ahead of the percentile.
- **Consent-gated** traits are severe and have no established intervention. They are
  withheld unless the person has opted in (`run.py consent`), or the operator turns the gate
  off for everyone (`PRS_S1_GATE=0`). See [`hardware-and-scope.md` §5](hardware-and-scope.md).

Percentiles are against the **European** 1000 Genomes reference (`KG_POP` in
`pipeline/config.py`). Most scores are European-derived and miscalibrate in other
ancestries; the report says so, and so should anyone reading a percentile.

## Adding a score

A score is one row in `panels/complex_traits.tsv`. Nothing else in the code names
individual scores, so adding a row is the whole change — but pick the score carefully first.

### 1. Choose the score

On the PGS Catalog page for the candidate, check:

- **A harmonized GRCh38 scoring file exists.** The pipeline downloads
  `<PGS id>_hmPOS_GRCh38.txt.gz` from the Catalog's `Harmonized/` folder and needs the
  columns `hm_chr`, `hm_pos`, `effect_allele`, `effect_weight`, and either `other_allele`
  or `hm_inferOtherAllele`. A file without usable rows fails loudly rather than scoring 0.
- **A performance metric for the score alone.** The panel's own bar is an AUROC or C-index
  with *no covariates*, in a European cohort of ≥ 1,000. Metrics that include age, sex or
  clinical covariates overstate what the genotype contributes.
- **The trait mapping is right.** Confirm the EFO term by reading it. Automated trait
  matching once mapped "carpal tunnel syndrome" to a body-height score with a perfectly
  plausible AUROC; nothing downstream can catch that.
- **The licence.** It is recorded, not enforced. Personal, non-commercial use satisfies
  every licence in the Catalog; redeploying commercially does not.
- **The size.** Above 50,000 variants the score is genome-wide: it needs `plink2` and its
  own 1000 Genomes reference run (step 3).

### 2. Add the row

Append a tab-separated line to `panels/complex_traits.tsv`:

| Column | What to put |
|---|---|
| `efo_id` | The ontology term, e.g. `MONDO_0005147` |
| `trait` | The label shown in the report |
| `pgs_id` | The Catalog id, e.g. `PGS000024` |
| `n_variants` | Variant count from the Catalog page (decides the engine and the tier cap) |
| `metric`, `metric_kind`, `eval_n` | The score-alone figure, `AUROC` or `C-index`, and the evaluation cohort size |
| `licence` | A short tag, e.g. `pgs-catalog-default`, `CC-BY-4.0` |
| `actionability` | `treat`, `screen`, `life` or `none` — or `review` if you have not decided |
| `s1_gate` | `yes` for a severe trait with no intervention, otherwise `no` |

A second score for a trait already in the panel is allowed; give it a distinguishable
`trait` label so the two rows can be told apart in the report.

### 3. Run it

```bash
uv run python run.py prs <sample>
```

No `--force`: every report records the panel it was scored against, so `prs` sees that the
panel changed and re-scores. You also do not have to run it for each person. The scheduled
`check-releases` job does the same comparison and re-scores **every callset that already
has a PRS report** — after the ClinVar scan if one is due, and on its own if not — then
rebuilds the dashboards and sends one notification. A callset that has never been scored is
left alone; run `prs` for it once and it is kept current from then on.
`PRS_AUTO_RESCORE=0` turns the scheduled re-scoring off.

What happens the first time a new score is seen:

- The scoring file is downloaded once into `$GENOMES_ROOT/tools/pgs/` and cached.
- **Small score:** European allele frequencies at its sites are queried from 1000 Genomes
  and cached. Quick.
- **Genome-wide score:** the whole 1000 Genomes panel is scored for it, chromosome by
  chromosome, and cached under `tools/pgs/<PGS id>_1kg/`. This is slow for the first person
  and reused for everyone after. `run.py mirror-1kg` (about 35 GB) makes it far faster than
  querying the remote copy. `--max-variants 50000` skips genome-wide scores for a quick run.

Scores are isolated from each other: if the new one fails, the report lists it under
"Scores that failed" and every other trait is unaffected.

### 4. Refresh the force-called sites (small scores only)

If you use `force-call`, its site list is built from the same panel, so a new small score's
positions are not in an existing force-called VCF. Re-run it — it notices the list changed
and re-calls without `--force` — then re-run the report. Here `--force` is needed: the
panel has not changed since the report was written, only the genotypes behind it.

```bash
uv run python run.py force-call <sample>
uv run python run.py prs <sample> --force
```

Until then those positions fall back to the ordinary calls, where an absent site is assumed
homozygous reference. The runbook's force-calling recipe explains the stage.

### 5. Check the result

- **Variants (used / total)** should be close to the total. A low count means the score's
  sites are largely missing from the sample.
- A percentile of `—` on a genome-wide score means the reference run did not complete; the
  raw score is still valid. Re-run; partial reference results are never cached.
- A percentile pinned at the extreme for every family member points at the comparison, not
  the family. Read the two bias traps in `prs-architecture.md` before trusting it.

## Regenerating the panel

`panels/build_complex_traits.py` rebuilds the whole file from the current PGS Catalog
release by the criteria in the design document:

```bash
uv run python panels/build_complex_traits.py
```

It **overwrites** `complex_traits.tsv`, so hand-added rows are lost, and the best score for
a trait can change between Catalog releases. Review the diff before committing it. To make a
curation decision survive regeneration, put it in the script: `DROP` removes a trait,
`ACTIONABILITY` assigns its section, `S1_GATE` gates it. A score that does not meet the
criteria cannot be kept by the generator at all — if you add one by hand, re-add it after
regenerating.

## Removing or replacing a score

Delete or edit the row. Any edit to the panel makes the existing reports stale, so they
are re-scored the same way as for a new score (`prs <sample>`, or the next
`check-releases`). Cached files for a removed score stay in
`tools/pgs/` until deleted by hand; they are harmless.
