-- Compound heterozygous: child carries TWO different het variants in the same gene, one
-- from each parent. Params: ${snapshot} ${child} ${father} ${mother}
--
-- WHAT THIS ESTABLISHES, AND WHAT IT DOES NOT
-- Each side now requires the variant to be present in one parent and ABSENT from the other.
-- Without that exclusion (as this query originally stood) a variant carried by both parents
-- satisfied "from_dad" and "from_mom" simultaneously, so two variants that both came from
-- the same parent — in cis, on one chromosome copy, clinically unremarkable — were reported
-- as a compound heterozygote.
--
-- Absence from a parent's rows still means "not called", not "reference": a variant-only
-- VCF omits hom-ref sites, uncovered sites and filter-failed calls alike. So this is a
-- CANDIDATE list. `run.py phase <child>` resolves cis/trans from reads directly and is what
-- should decide the question; `run.py segregation` force-calls the parents at these sites.
--
-- Restricted to ClinVar-actionable variants: any (Likely) pathogenic TERM of the compound
-- CLNSIG value, the same expression pipeline/clinsig.py uses (`sql_actionable`).
-- Without that filter the gene-level pair join is a cartesian product over every het the
-- child carries: 274 million rows as this query originally stood, 57 million even with the
-- parental-origin exclusions above. Neither is an answer to anything.
WITH v AS (SELECT * FROM variants WHERE snapshot_id = '${snapshot}'
           AND regexp_matches(clinvar_sig,
               '(^|[/|])(Likely_pathogenic|Pathogenic)(,_low_penetrance)?($|[/|])')),
child AS (SELECT * FROM v WHERE sample_id='${child}' AND zygosity='HET'),
from_dad AS (
  SELECT c.* FROM child c
  JOIN v f ON f.sample_id='${father}'
    AND f.chrom=c.chrom AND f.pos=c.pos AND f.ref=c.ref AND f.alt=c.alt
  WHERE NOT EXISTS (SELECT 1 FROM v m WHERE m.sample_id='${mother}'
    AND m.chrom=c.chrom AND m.pos=c.pos AND m.ref=c.ref AND m.alt=c.alt)),
from_mom AS (
  SELECT c.* FROM child c
  JOIN v m ON m.sample_id='${mother}'
    AND m.chrom=c.chrom AND m.pos=c.pos AND m.ref=c.ref AND m.alt=c.alt
  WHERE NOT EXISTS (SELECT 1 FROM v f WHERE f.sample_id='${father}'
    AND f.chrom=c.chrom AND f.pos=c.pos AND f.ref=c.ref AND f.alt=c.alt))
SELECT d.gene,
       d.chrom||':'||d.pos||' '||d.ref||'>'||d.alt AS paternal_variant,
       d.clinvar_sig AS paternal_sig,
       m.chrom||':'||m.pos||' '||m.ref||'>'||m.alt AS maternal_variant,
       m.clinvar_sig AS maternal_sig
FROM from_dad d JOIN from_mom m ON d.gene=m.gene
WHERE NOT (d.chrom=m.chrom AND d.pos=m.pos AND d.ref=m.ref AND d.alt=m.alt)
  AND d.gene IS NOT NULL
ORDER BY d.gene;
