-- X-linked recessive: variant on chrX in a son, carried heterozygously by the mother.
-- Params: ${snapshot} ${son} ${mother}
--
-- The mother-is-HET condition is now enforced. The header always described it, but the
-- query only joined on the mother carrying the variant at all, so a mother homozygous for
-- it — a different situation entirely — matched the same way.
--
-- The son's zygosity is reported, not filtered on: a hemizygous male X call surfaces as HET
-- or HOM depending on the caller's ploidy model, so neither value can be required without
-- silently dropping true hits. Read son_zyg rather than assuming it. Pseudoautosomal
-- regions are genuinely diploid in males and are NOT excluded here.
-- Restricted to ClinVar-actionable variants: any (Likely) pathogenic TERM of the compound
-- CLNSIG value, the same expression pipeline/clinsig.py uses (`sql_actionable`).
-- Unfiltered this returns every chrX variant the son shares with a het mother — tens of
-- thousands of rows, which answers nothing.
WITH v AS (SELECT * FROM variants WHERE snapshot_id = '${snapshot}'
           AND regexp_matches(clinvar_sig,
               '(^|[/|])(Likely_pathogenic|Pathogenic)(,_low_penetrance)?($|[/|])'))
SELECT s.gene, s.chrom, s.pos, s.ref, s.alt, s.clinvar_sig, s.consequence, s.gnomad_af,
       s.zygosity AS son_zyg, m.zygosity AS mother_zyg
FROM v s
JOIN v m ON m.sample_id='${mother}' AND m.chrom=s.chrom AND m.pos=s.pos
        AND m.ref=s.ref AND m.alt=s.alt AND m.zygosity='HET'
WHERE s.sample_id='${son}' AND s.chrom IN ('X','chrX')
ORDER BY (s.clinvar_sig IS NULL), s.gene;
