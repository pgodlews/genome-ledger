# PRS subsystem — architecture

How the polygenic-risk-score stage works, written for someone who reads code, not chromatograms.
All of it lives in [`pipeline/prs.py`](../pipeline/prs.py); everything else is config,
a panel file, and two external data sources.

## The 60-second version, in IT terms

A PGS (polygenic score) is a **linear model with no intercept**, published as a flat file:

| Genomics term | What it actually is |
|---|---|
| Scoring file (`PGS000018_hmPOS_GRCh38.txt.gz`) | A TSV of `(feature_key, expected_value, weight)` — one row per genome position |
| Variant / position `chr9:22125504` | The primary key. Two of them per file column: `hm_chr`, `hm_pos` |
| Effect allele | Which of the two possible values at that key counts |
| Dosage | The feature value: `0`, `1` or `2` copies of the effect allele |
| Raw score | `Σ dosage × weight` — a plain dot product |
| VCF | A **sparse delta** of the sample against the reference genome. Rows exist only where the sample differs |
| Hom-ref | "Key absent from the delta" ⇒ the sample has the reference value, i.e. two copies of it |
| Percentile | The rank of this dot product inside a population of 3202 reference individuals |
| 1000 Genomes (1000G) | That reference population, published as remote VCFs with per-population frequencies |
| GRCh38 / hg38 | The coordinate system. Everything here is GRCh38; mixing builds is a silent key-mismatch bug |

Two properties of the input drive almost every design decision in this module:

1. **The VCF is sparse.** Absence means hom-ref, not "unknown". Where the *effect allele is the
   reference allele*, an absent row means dosage **2**, not 0 — the single most common source of
   bias in this code (see `_ref_effect_supplement`).
2. **Score sizes span five orders of magnitude.** 14 variants (endometriosis) to ~7.3M (cutaneous melanoma). One engine
   cannot serve both, so there are two, split at `PRS_PLINK_MIN = 50_000`.

## Component map

```mermaid
flowchart TB
  subgraph entry["Entry point"]
    CLI["run.py prs SAMPLE --force"]
    REP["prs_report()"]
  end

  subgraph inputs["Local inputs"]
    PANEL["panels/complex_traits.tsv<br/>44 traits: pgs_id, actionability,<br/>s1_gate, licence, AUROC"]
    CFG["pipeline/config.py<br/>PRS_PLINK_MIN, PRS_MAX_VARIANTS,<br/>PRS_S1_GATE, KG_POP=EUR"]
    VCF["sample VCF<br/><i>call.sample_vcf()</i>"]
    FC["force-called VCF<br/><i>forcecall.py</i> — optional,<br/>preferred for small scores"]
    FA["reference FASTA<br/><i>util.faidx_bases</i>"]
  end

  subgraph remote["External data — cached on first use"]
    PGSC["PGS Catalog FTP<br/>harmonized GRCh38 scoring files"]
    KG["1000 Genomes phased VCFs<br/>remote HTTP or local 34.6GB mirror<br/><i>config.kg_vcf_source()</i>"]
    MAP["1000G population map<br/>sample id to superpopulation"]
  end

  subgraph engines["Two scoring engines"]
    EX["exact path — under 50k variants<br/>bcftools query per position<br/><i>_score_one</i>"]
    PL["plink2 path — genome-wide<br/>plink2 --score over a pgen<br/><i>_score_plink2</i>"]
  end

  subgraph calib["Two percentile methods"]
    AN["analytic<br/>Normal(sum 2fw, sum 2f(1-f)w^2)<br/>HWE + CLT<br/><i>_calibrate</i>"]
    EMP["empirical<br/>score all 3202 individuals,<br/>count those below<br/><i>_reference_scores</i>"]
  end

  subgraph caches["Caches — ~/genomes/tools"]
    C1["pgs/PGSxxxx_hmPOS_GRCh38.txt.gz"]
    C2["pgs/PGSxxxx_1kg_af.tsv"]
    C3["pgs/PGSxxxx_1kg/scores_all.tsv<br/>+ per-chrom reference pgen"]
    C4["plink2/SAMPLE.pgen<br/>stamped on VCF identity"]
  end

  OUT["reports/md/genome/PERSON/prs_SAMPLE.md<br/>grouped by actionability<br/>(+ HTML and JSON renderings)"]

  CLI --> REP
  PANEL --> REP
  CFG --> REP
  REP --> EX
  REP --> PL
  VCF --> EX
  FC --> EX
  FA --> EX
  VCF --> C4 --> PL
  PGSC --> C1 --> EX
  C1 --> PL
  EX --> AN
  PL --> EMP
  KG --> C2 --> AN
  KG --> C3 --> EMP
  MAP --> EMP
  AN --> OUT
  EMP --> OUT

  classDef ext fill:#2d3f52,stroke:#5b7ea8,color:#e8eef5
  classDef cache fill:#3d3524,stroke:#8a7440,color:#f2ece0
  class PGSC,KG,MAP ext
  class C1,C2,C3,C4 cache
```

## End-to-end sequence

One `prs` run scores every row of the panel in a loop. Each row is isolated: a 404 or a truncated
download costs that one trait, not the run (`prs_report` in [`prs.py`](../pipeline/prs.py)).

```mermaid
sequenceDiagram
  autonumber
  actor U as CLI
  participant P as prs_report
  participant PF as panel + filters
  participant SF as _scoring_file
  participant SC as _score_entry
  participant BC as bcftools / plink2
  participant KG as 1000G
  participant R as report writer

  U->>P: run.py prs SAMPLE001 --force
  P->>P: require_tools bcftools, samtools
  P->>P: VCF exists? report exists, stamps current, and not --force?
  P->>P: forcecall.available(sample)?<br/>if yes, small scores read the force-called VCF
  P->>PF: _panel_scores() + _panel_filters()
  PF-->>P: N rows, minus S1-gated, minus over PRS_MAX_VARIANTS

  loop for each trait row (failures isolated)
    P->>SC: _score_entry(sample, row, gt_vcf, out_dir)
    SC->>SF: scoring file for PGS id
    SF->>SF: cache hit? else download from PGS Catalog FTP
    SF-->>SC: list of (chrom, pos, effect, other, weight)
    Note over SC: zero parsed rows raises —<br/>a clean-looking 0.0 over 0 variants<br/>is the worst possible report row

    alt variants <= PRS_PLINK_MIN (50k)
      SC->>BC: bcftools query -R positions, PASS-only
      BC-->>SC: genotypes at present positions
      SC->>BC: samtools faidx for absent positions
      BC-->>SC: reference base ⇒ hom-ref dosage
      SC->>KG: effect-allele frequencies at these positions
      KG-->>SC: AF_EUR / AF (cached per PGS)
      SC->>SC: analytic percentile, same variant subset both sides
    else variants > 50k (genome-wide)
      SC->>BC: build/reuse sample pgen, then plink2 --score
      BC-->>SC: raw SCORE1_SUM
      SC->>KG: extract 1000G genotypes per chromosome
      KG-->>SC: reference pgen (cached per PGS, family-wide)
      SC->>BC: plink2 --score over 3202 individuals
      BC-->>SC: per-individual reference scores
      SC->>SC: add _ref_effect_supplement, then rank sample among them
    end
    SC-->>P: row + raw + percentile + method
  end

  P->>R: tables grouped by actionability, caveats, failures
  R-->>U: reports/md/genome/PERSON/prs_SAMPLE.md
```

## Engine routing and the two bias traps

The interesting part is not the dot product — it's making the sample's score and the reference
distribution cover **the same variant set**. Both branches get this wrong in opposite directions
if you are careless, and both corrections are in the code with the incident written into the comment.

```mermaid
flowchart TD
  START["variants from scoring file"] --> N{"len(variants)<br/>&gt; PRS_PLINK_MIN<br/>(50,000)?"}

  N -->|no — exact path| E1["bcftools query, PASS or '.' only"]
  E1 --> E2["position present in VCF?"]
  E2 -->|yes| E3["dosage from GT"]
  E2 -->|no| E4["faidx reference base<br/>⇒ dosage of ref allele x2"]
  E3 --> E5["contrib[key] = dosage x weight"]
  E4 --> E5
  E5 --> E6["1000G AF per position<br/>cached per PGS"]
  E6 --> E7["percentile from<br/>Normal(mu, sigma^2) under HWE+CLT"]
  E7 --> E8["TRAP: sum the sample score and the<br/>population moments over the SAME subset —<br/>only positions with both a genotype and an AF"]

  N -->|yes — plink2 path| G1["pgen from VCF<br/>--var-filter --autosome --rm-dup<br/>stamped on VCF identity + recipe"]
  G1 --> G2["plink2 --score, no-mean-imputation"]
  G2 --> G3["reference: extract 1000G genotypes<br/>per chromosome, build pgen, score all 3202"]
  G3 --> G4["TRAP: score both sides over the FULL<br/>variant set, autosomes 1-22 only"]
  G4 --> G5["plink2 cannot see hom-ref sites at all.<br/>Where effect allele = REF allele,<br/>absent means TWO copies, not zero"]
  G5 --> G6["_ref_effect_supplement adds 2 x weight back,<br/>reading REF from the reference pvar"]
  G6 --> G7["empirical percentile =<br/>share of 3202 scoring below the sample"]

  E8 --> DONE["row: raw, used/total, EUR %ile, band"]
  G7 --> DONE
```

Why the traps matter, concretely: restricting the reference scoring to the sample's own covered
sites is textbook **ascertainment bias**. A variant-only VCF's covered sites are exactly the sites
where the sample is non-reference, so the sample carried alt alleles at every scored site by
construction while the reference individuals carried them at population frequency. Every sample
lands in the top percentile — and two unrelated people both in the top 1% of one score has a
prior probability of about 1 in 10,000, which is how the bug shows itself. That is what the
`_ref_effect_supplement` design replaced.

## Caches and their invalidation

Everything expensive is cached; most caches are **sample-independent**, which is why the second
family member scores in minutes.

| Cache | Path | Keyed on | Invalidation | Shared across |
|---|---|---|---|---|
| Scoring file | `tools/pgs/PGSxxxx_hmPOS_GRCh38.txt.gz` | PGS id | manual delete | everyone |
| 1000G allele freqs | `tools/pgs/PGSxxxx_1kg_af.tsv` | PGS id | manual delete; **never written partial** | everyone |
| Reference pgen + scores | `tools/pgs/PGSxxxx_1kg/` | PGS id | manual delete; **never written partial** | whole family |
| Sample pgen | `tools/plink2/SAMPLE.pgen` | `file_token(VCF)` + `_PGEN_RECIPE` | automatic — re-call or re-ingest rebuilds it | one sample |
| Report | `reports/md/genome/PERSON/prs_SAMPLE.md` | panel digest + variant cap (`prs/SAMPLE/prs_SAMPLE.panel-stamp`) and S1 consent state (`.consent-stamp`) | automatic — a panel edit or a consent change re-scores it, on the next `prs` or the next `check-releases`; `--force` for anything else | one sample |

The "never written partial" rule appears three times in this module and is load-bearing. A reference
distribution missing one chromosome is not a degraded percentile, it is a **wrong** one: the sample
is still scored genome-wide, so every reference individual sits artificially low and the sample's
percentile is inflated — permanently, family-wide, with no invalidation path short of deleting the
file by hand. Hence: partial result ⇒ log a warning, return it for this run, refuse to cache.

```mermaid
stateDiagram-v2
  [*] --> Check
  Check --> Reuse: pgen exists and stamp matches VCF token + recipe version
  Check --> Build: missing, or VCF changed, or _PGEN_RECIPE bumped
  Build --> Rename: plink2 --make-pgen into a per-process temp prefix
  Rename --> Stamp: os.replace parts into place, .pgen last
  Stamp --> Reuse: write_stamp is written LAST, so a crash mid-build leaves it stale
  Reuse --> [*]
```

## Report shape

The output is ordered by **actionability, not by score** — a 60th-percentile result with a screening
pathway outranks a 95th for something nothing can be done about
(`_ACTION_SECTIONS` in [`prs.py`](../pipeline/prs.py)). Sections: `treat` → `screen` → `life` → `none` → `review`,
each sorted by percentile within itself. Every row carries the score's own published AUROC, because
that bounds what any percentile can mean: a 0.55 score separates cases from controls barely better
than a coin flip, whatever percentile it returns.

Two filters run before scoring (`_panel_filters` in [`prs.py`](../pipeline/prs.py)):

- **`PRS_S1_GATE`** — withholds traits flagged severe *and* non-actionable from anyone who has not
  opted in (`run.py consent`). On by default; `PRS_S1_GATE=0` or `prs --no-s1-gate` turns it off for
  everyone. The report is stamped with the state it was written under
  (`prs/<sample>/prs_<sample>.consent-stamp`): a cached report is reused only while that state
  still holds, and withdrawing consent deletes it in every format.
- **`PRS_MAX_VARIANTS`** — a tier cap. Half the panel is genome-wide and each such score needs its
  own cached 1000G reference run, so the cap lets the ~22 cheap scores land in minutes while the
  genome-wide half backfills separately.

## Where the files are

| Concern | File |
|---|---|
| All scoring logic | [`pipeline/prs.py`](../pipeline/prs.py) |
| CLI wiring | the `prs` and `force-call` commands in [`run.py`](../run.py) |
| Tunables, paths, 1000G source selection | the `PRS_*`, `PGS_*` and `KG_*` settings in [`pipeline/config.py`](../pipeline/config.py) |
| Trait panel (data) | [`panels/complex_traits.tsv`](../panels/complex_traits.tsv) |
| Panel generator + selection criteria | [`panels/build_complex_traits.py`](../panels/build_complex_traits.py), [`docs/polygenic-panel-design.md`](polygenic-panel-design.md) |
| Force-call integration | `_prs_sites` in [`pipeline/forcecall.py`](../pipeline/forcecall.py) |
| Where the report is filed and listed | `reportpaths.genome_report("prs", …)`; `GENOME_REPORTS` in [`pipeline/render.py`](../pipeline/render.py) |
| Which scores, and adding one | [`docs/prs-guide.md`](prs-guide.md) |
| Tests | [`tests/test_units.py`](../tests/test_units.py) |

External dependencies at runtime: `bcftools`, `samtools` (both required up front), `plink2`
(downloaded by `setup`, only needed for genome-wide scores), and network access to the PGS Catalog
and 1000G FTP unless both are already mirrored locally.
