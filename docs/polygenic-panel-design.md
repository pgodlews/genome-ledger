# Polygenic Panel — Design & Trait Selection

How the complex-disease (polygenic) report picks its traits, what it does differently from
the simple lead-SNP scores common in consumer genomics, and how the scores it uses are
licensed.

Companion to [`prs-guide.md`](prs-guide.md) (the list of scores and how to add one),
[`prs-architecture.md`](prs-architecture.md) (how scoring works) and the `prs` pipeline
module. This doc covers **which traits, chosen how, and on whose evidence**.

## 1. The lead-SNP approach, and where it falls short

A widely used way to produce a "polygenic" result for a consumer report is the
**single-GWAS lead-SNP additive score**: take one landmark genome-wide association study,
take its genome-wide-significant lead SNPs, weight each by the published per-allele effect
size, sum, and place the result in a band against some comparison group.

The method is legible and easy to reproduce from public GWAS summary tables, which is why
it is common. It is also about a decade behind the state of the art.

### Where it is weak

| Weakness | Consequence |
|---|---|
| Lead SNPs only, one GWAS | Discards genome-wide signal; no LD-aware weighting (no LDpred/PRS-CS) |
| Convenience reference group | A percentile against a self-selected user base is not a percentile against a population cohort |
| Coarse bands | A handful of bins hides most of the distribution, and their labels tend to understate the extremes |
| No absolute risk | A percentile alone can't inform a screening or treatment decision |
| No coverage/QC | Rarely states what fraction of the score's variants were actually callable |
| Siloed | Polygenic, monogenic and PGx results never combine |
| No inclusion floor | A trait can be reported on a handful of loci with no validated score behind it |

That last row is the biggest quality problem. A percentile computed from half a dozen loci
is noise presented as a finding, and nothing in the method itself stops it being reported.

## 2. Seven concrete improvements

1. **Validated scores, not lead-SNP sums.** Use PGS Catalog harmonized scoring files with
   published performance metrics, which is what the `prs` stage does.
2. **Curated trait→score mapping, never fuzzy-matched.** See §4 — this is a correctness
   requirement, not a nicety.
3. **Continuous ancestry-matched percentile with an interval**, against 1000G, not bands
   against a convenience sample. Already built; keep the EUR-matched *and* overall columns.
4. **Absolute risk, not just percentile.** Convert percentile → lifetime/10-year risk using
   population incidence (GBD / national registries) plus the score's effect size per SD.
   This is the single largest differentiator: it is what turns "97th percentile" into a
   decision. Requires per-trait incidence tables and per-SD effect sizes — real work,
   and the honest blocker is that many Catalog `OR` entries don't state whether they are
   per-SD or per-quantile, so this needs curation per trait, not a bulk import.
5. **Report coverage and QC per score.** What fraction of the score's variants were callable
   in this sample, how many were imputed as hom-ref, and what that does to the percentile.
   The `force-call` stage already produces the evidence needed.
6. **Combine the evidence layers.** A pathogenic monogenic variant *and* a high polygenic
   score is a different finding than either alone; polygenic background is known to modify
   monogenic penetrance. The pipeline holds monogenic, PGx and pedigree data in one place —
   a report built from one result type at a time cannot do this.
7. **Use the family as a calibration check.** Parent–offspring PRS correlation should be
   ≈0.5 for an additive score. Trios that deviate are a bug signal — a free, continuous
   correctness test that is unavailable when only one person is sequenced.

## 3. Trait inclusion criteria

The panel is **generated from these criteria**, not hand-picked and not copied. A trait
enters `panels/complex_traits.tsv` when all hold:

- ≥1 PGS Catalog score with a **PGS-alone** performance metric — AUROC or C-index — with the
  `Covariates Included in the Model` field **empty**
- evaluated in an ancestry-matched cohort (currently `Broad Ancestry Category = European`)
  of **≥1,000 individuals**
- metric strictly between 0.5 and 1.0
- trait passes curation: not a proxy phenotype, ontology artefact, or redundant subtype

Applied to the current PGS Catalog release this yields **70 candidate traits**, curated
down to **44**. That curation step removed ontology artefacts the criteria alone let through
(`female` and `male` as "traits", `comparative body size at age 10, self-reported`), six
redundant breast-cancer receptor subtypes collapsed into one, and three juvenile idiopathic
arthritis sub-phenotypes.

**Covariate strictness matters.** Relaxing it inflates everything: type 1 diabetes appears
at AUROC 0.96 — but that evaluation includes autoantibodies and family history as
covariates. The PGS-alone figure is 0.90, still the best in the panel, and it is the honest
one.

### What the criteria bring in

Applying the criteria surfaces high-evidence traits that a lead-SNP panel has no natural
way to rank:

| Trait | PGS-alone | Why it matters |
|---|---|---|
| **Breast carcinoma** | C-index 0.673 | `PGS000004` (313-variant) is among the best-validated scores in existence and feeds established clinical models |
| **Celiac disease** | AUROC 0.900 | HLA-driven, near-diagnostic negative predictive value, trivially actionable |
| **Venous thromboembolism** | C-index 0.670 | Directly actionable (prophylaxis, contraceptive choice) |
| **LDL / HDL / triglycerides / urate** | 0.60–0.68 | Quantitative, measurable, treatable — and checkable against an actual blood test |
| **Blood pressure / hypertensive disorder** | 0.69–0.72 | Largest modifiable cardiovascular contributor |
| **Ovarian, endometrial carcinoma; endometriosis** | 0.56–0.62 | Female-specific traits with validated scores |
| **Body mass index** | AUROC 0.739 | Strong score; frames every metabolic result |

### What the criteria keep out

- **Traits with no validated score.** If PGS Catalog holds no score for a trait with a
  PGS-alone metric in a matched cohort, the trait is not reported — however familiar it is.
- **Traits resting on a handful of loci.** A few genome-wide-significant hits do not make a
  usable predictor.
- **Nuisance and self-report phenotypes**, ontology artefacts, and redundant subtypes
  (removed by the curation lists in `panels/build_complex_traits.py`).

Sometimes the right move is to *re-base* rather than drop: where a disease has no good score
of its own but a closely related measurable trait does (bone density for osteoporosis, a
venous-thromboembolism score for deep vein thrombosis), the validated score is the one to
use.

### Gate, don't drop

Non-actionable severe conditions (multiple sclerosis, renal/bladder/pancreatic carcinoma)
stay in the panel but carry `s1_gate = yes`, deferring to the S1 incidental-findings consent
filter in [`hardware-and-scope.md §5`](hardware-and-scope.md). Withheld by default; each person
opts in (`run.py consent`), or an operator turns the gate off with `PRS_S1_GATE=0`. Evidence quality and disclosure policy are separate
axes and must not be conflated.

## 4. Why mapping must be curated

An automated trait→EFO mapper using fuzzy string matching produced these, silently:

```
"Carpal tunnel syndrome"  ->  OBA_VT0001253 (body height)   AUROC 0.861
"Colon polyp"             ->  OBA_VT0001253 (body height)   AUROC 0.861
```

Both inherited a *short-stature-in-females* score and a plausible-looking metric. Nothing
downstream would have caught it — the numbers look fine, the report renders, the percentile
is real. It is simply the wrong trait.

**Therefore:** `panels/complex_traits.tsv` pins an explicit `efo_id` **and** `pgs_id` per
row, both human-reviewed. Fuzzy matching may *propose* rows; it must never *populate* them.
Score selection is pinned by ID and versioned alongside the Catalog release, so a report is
reproducible and a score swap shows up as a reviewable diff.

## 5. Provenance and licensing

**The trait list is generated, not hand-picked.** `panels/complex_traits.tsv` is produced by
the §3 criteria running over PGS Catalog metadata; `panels/build_complex_traits.py` is the
whole procedure and needs nothing but the Catalog to run. No third-party trait list is an
input. Two house rules follow from that:

- Don't name commercial vendors in code, panel files, or user-facing output.
- Don't describe the tool as compatible with, equivalent to, or a replacement for any
  commercial product.

**Per-score licensing: record it, don't filter on it.** Of 6,972 PGS Catalog scores, 6,879
carry the default terms; 31 are CC BY-NC-ND and ~9 are academic- or research-only.
Treating those as disqualifying would be wrong, and filtering on them is actively harmful:

- **NonCommercial is satisfied by definition.** Analysing your own genome is not commercial use.
- **NoDerivatives restricts *distribution*.** `prs.py` fetches scoring files at runtime and
  never redistributes them, so the trigger never fires — and a score is arguably not a
  derivative work of a weight table regardless.
- The one restriction that does bite is **downstream**: a fork operating commercially would
  need to clear those ~40 scores. That is the operator's compliance question, not a reason
  to remove capability from every personal user.

Filtering would cost real coverage: it would exclude one panel trait, **inflammatory bowel
disease** (`PGS000017`, 6.9M variants, AUROC 0.69), whose licence reads *"freely available
to the academic community for research use"* with a contact address **for commercial
purposes** — i.e. explicitly permitting the use this project makes of it.

So the panel carries a `licence` column (`pgs-catalog-default` / `CC-BY-4.0` / `CC0` /
`CC-BY-NC-ND` / `academic-noncommercial` / `research-only`) and `prs.py` should surface it
in the report. It is provenance, not a gate. Current panel: 42 default, 1 CC-BY-4.0,
1 academic-noncommercial.

**Attribution.** Cite the PGS Catalog and GWAS Catalog per their terms, record the release
version used, and keep fetching scoring files at runtime rather than vendoring them — that
sidesteps redistribution entirely and is what `prs.py` already does.

**Regulatory posture** is unchanged from `hardware-and-scope.md §4`: research/educational
use, no clinical claims, no diagnostic framing. Adding absolute-risk output (§2.4) pushes
closer to the line that EU MDR Rule 11 / IVDR care about, so it must ship with explicit
"not for medical decision-making" framing and confirmation-by-clinician language.

## 6. Seed panel

`panels/complex_traits.tsv` — 44 traits, criteria-derived, columns:

| Column | Meaning |
|---|---|
| `efo_id`, `trait` | pinned ontology term, human-reviewed |
| `pgs_id`, `n_variants` | pinned score |
| `metric`, `metric_kind`, `eval_n` | PGS-alone AUROC/C-index and evaluation cohort size |
| `licence` | score's terms of use — provenance for commercial redeployment, not a gate |
| `actionability` | `screen` / `treat` / `life` / `none` — drives report prominence |
| `s1_gate` | defer to the incidental-findings consent filter |

Distribution: 22 `screen`, 15 `treat`, 3 `life`, 4 `none`. Every row is actionable or
explicitly gated — which is the point.

**Not yet done:** `actionability` is a first-pass judgement and needs review against actual
guidelines; absolute-risk conversion (§2.4) needs per-trait incidence and per-SD effect
sizes; surfacing the `licence` column in rendered reports is specified here but not
implemented.
