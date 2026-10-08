# Demo Multi-Generational Family VCF Fixtures

Synthetic, lightweight, and fully functional whole-genome VCF fixtures for a 3-generation family cohort (Grandparents, Parents/Daughter & Marry-in Partner, and Grandchildren).

Designed for:
- Zero-data demonstration of Genome Ledger's multi-generational reports, reproductive risk, and workflows (`run.py demo`).
- Fast, self-contained unit and regression testing without requiring gigabyte-scale sequencing data.
- Demonstrating custody handoff, marry-in founder genetics (partner with no parental history), and longitudinal stewardship.
- GRCh38 assembly compliance.

## Included Samples (3 Generations, 6 Members)

| Generation | Sample ID | Role | Sex | File |
|---|---|---|---|---|
| **Gen 1** | `father` | Grandfather | Male | `father.vcf.gz` |
| **Gen 1** | `mother` | Grandmother | Female | `mother.vcf.gz` |
| **Gen 2** | `daughter` (`child`) | Daughter (Mother in Gen 2) | Female | `daughter.vcf.gz` |
| **Gen 2** | `partner` | Partner (Father in Gen 2 — Marry-in Founder) | Male | `partner.vcf.gz` |
| **Gen 3** | `grandson` | Grandson | Male | `grandson.vcf.gz` |
| **Gen 3** | `granddaughter` | Granddaughter | Female | `granddaughter.vcf.gz` |

Accompanying metadata:
- `manifest.tsv`: Ingest manifest for `python run.py ingest`.
- `family.tsv`: The family graph (one row per person); `run.py demo` copies it in and generates `pedigree.ped` from it.

---

## Seeded Variants & Expected Report Outcomes

### 1. Monogenic Panel (`panels/monogenic.tsv`)

#### A. Familial Hypercholesterolemia (`APOB` Arg3527Gln, rs5742904, `chr2:21002409 G>A`)
- **Inheritance**: Autosomal Dominant (AD)
- **Clinical context**: Elevated LDL cholesterol, non-fatal, highly actionable with routine statins/lifestyle.
- **Genotypes**:
  - `father`: `0/1` (At-Risk)
  - `mother`: `0/0` (Normal)
  - `daughter`: `0/1` (At-Risk, inherited from father)
  - `partner`: `0/0` (Normal marry-in founder)
  - `grandson`: `0/1` (At-Risk, inherited from daughter)
  - `granddaughter`: `0/0` (Normal, spared)
- **Report Impact**:
  - Highlights **3-Generation 50% Mendelian AD Transmission**: Passed from Grandfather ➔ Daughter ➔ Grandson, while Granddaughter inherits the normal allele.

#### B. Hereditary Hemochromatosis (`HFE` C282Y, rs1800562, `chr6:26093141 G>A`)
- **Inheritance**: Autosomal Recessive (AR)
- **Genotypes**:
  - `father`: `0/1` (Carrier)
  - `mother`: `0/1` (Carrier)
  - `daughter`: `1/1` (Affected)
  - `partner`: `0/1` (Carrier marry-in founder)
  - `grandson`: `1/1` (Affected)
  - `granddaughter`: `0/1` (Carrier)
- **Report Impact**:
  - Gen 1: Both grandparents carry the allele ➔ Daughter is affected (`1/1`).
  - Gen 2: Daughter (`1/1`) pairs with Carrier partner (`0/1`) ➔ Flags **50% child risk** in `panel_family.md`.
  - Gen 3: Grandson is affected (`1/1`), Granddaughter is carrier (`0/1`).

#### C. Cystic Fibrosis (`CFTR` G542X, rs113993959, `chr7:117548628 G>A`)
- **Inheritance**: Autosomal Recessive (AR)
- **Genotypes**:
  - `mother`: `0/1` (Carrier)
  - `daughter`: `0/1` (Carrier)
  - `partner`: `0/0` (Normal marry-in founder)
  - `grandson`: `0/1` (Carrier)
  - `granddaughter`: `0/0` (Normal)
- **Report Impact**:
  - Shows carrier screening peace of mind: Daughter is a carrier, but because partner is non-carrier, their children are safe from the disease.

### 2. Single-SNP Traits (`traits/traits.tsv`)
- **Lactose digestion (adult)**: rs4988235 on `chr2:135851076 G>A`
  - `father` `G/A`, `mother` `G/G`, `daughter` `G/A`, `partner` `G/A`, `grandson` `G/G` (intolerant), `granddaughter` `G/A` (persistent).
- **Eye colour (blue vs brown)**: rs12913832 on `chr15:28120472 A>G`
  - `father` `A/G` (hazel), `mother` `G/G` (blue), `daughter` `G/G` (blue), `partner` `A/G` (hazel), `grandson` `A/G` (hazel), `granddaughter` `G/G` (blue).
- **Muscle fibre type (ACTN3)**: rs1815739 on `chr11:66560624 C>T`
  - Shows transmission of sprint (`C`) vs endurance (`T`) alleles across generations.
- **Alcohol flush reaction**: rs671 on `chr12:111803962 G>A`
  - Maternally transmitted ALDH2 deficiency from `mother` ➔ `daughter` ➔ `granddaughter`.
- **Earwax type & body odor**: rs17822931 on `chr16:48224287 C>T`
  - Wet earwax (`C/T`) segregation.

### 3. De Novo Mutation Candidate
- Private Heterozygous Mutation in `grandson`: `chr1:10000000 A>T` (`0/1`, DP=36, GQ=60), absent in `daughter` and `partner` (`0/0`).

### HLA tag SNP

`rs2187668` (`chr6:32638107 C>T`, tags the celiac DQ2.5 haplotype) is heterozygous in
`father`, `daughter` (`child`) and `grandson`, and absent in the others — so the `hla`
report, which the demo runs for real, shows a transmitted haplotype.

## Illustrative results (`../demo/toy_results.json`)

Repeat expansions, pharmacogenomics, haplogroups and polygenic scores cannot be computed
from these VCFs. `run.py demo` takes hand-written sample values for them from
`../demo/toy_results.json` and passes them through the real report code (see
`pipeline/demo.py`). They follow inheritance within the family but are not derived from the
genotypes here.
