-- Segregation: actionable variants and which family members carry them, so you can
-- check co-segregation with an affected phenotype by eye against pedigree.ped.
-- Params: ${snapshot}
SELECT gene, chrom, pos, ref, alt, clinvar_sig,
       string_agg(sample_id || '(' || zygosity || ')', ', ' ORDER BY sample_id) AS carriers,
       count(*) AS n_carriers
FROM variants
WHERE snapshot_id = '${snapshot}'
  AND clinvar_sig IN ('Pathogenic','Likely_pathogenic','Pathogenic/Likely_pathogenic')
GROUP BY gene, chrom, pos, ref, alt, clinvar_sig
ORDER BY n_carriers DESC, gene;
