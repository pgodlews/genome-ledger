-- Autosomal recessive: child homozygous-alt, both parents heterozygous carriers.
-- Params: ${snapshot} ${child} ${father} ${mother}
WITH v AS (SELECT * FROM variants WHERE snapshot_id = '${snapshot}')
SELECT c.gene, c.chrom, c.pos, c.ref, c.alt, c.clinvar_sig, c.consequence, c.gnomad_af
FROM v c
JOIN v f ON f.sample_id='${father}' AND f.chrom=c.chrom AND f.pos=c.pos
        AND f.ref=c.ref AND f.alt=c.alt AND f.zygosity='HET'
JOIN v m ON m.sample_id='${mother}' AND m.chrom=c.chrom AND m.pos=c.pos
        AND m.ref=c.ref AND m.alt=c.alt AND m.zygosity='HET'
WHERE c.sample_id='${child}' AND c.zygosity='HOM'
ORDER BY (c.clinvar_sig IS NULL), c.gene;
