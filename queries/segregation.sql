-- Segregation: actionable variants and which family members carry them, so you can
-- check co-segregation with an affected phenotype by eye against pedigree.ped.
-- Any (Likely) pathogenic TERM of a compound CLNSIG counts, the same expression
-- pipeline/clinsig.py uses (`sql_actionable`). Params: ${snapshot}
SELECT gene, chrom, pos, ref, alt, clinvar_sig,
       string_agg(sample_id || '(' || coalesce(zygosity, '?') || ')', ', ' ORDER BY sample_id) AS carriers,
       count(*) AS n_callsets
FROM variants
WHERE snapshot_id = '${snapshot}'
  AND regexp_matches(clinvar_sig,
      '(^|[/|])(Likely_pathogenic|Pathogenic)(,_low_penetrance)?($|[/|])')
GROUP BY gene, chrom, pos, ref, alt, clinvar_sig
ORDER BY n_callsets DESC, gene;
