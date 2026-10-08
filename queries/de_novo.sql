-- CANDIDATE de novo: present in the child, ABSENT FROM THE PARENTS' CALLS.
-- Params: ${snapshot} ${child} ${father} ${mother}
--
-- THIS IS NOT A DE-NOVO CALL SET. "Absent from a parent's rows" is not "the parent is
-- reference": a variant-only VCF omits hom-ref sites, uncovered sites and filter-failed
-- calls alike, and the three are indistinguishable here. Most rows this returns are
-- inherited variants the parent's callset simply did not emit.
--
-- `run.py denovo <child>` is the real implementation: it force-calls both parents at each
-- candidate site off their CRAMs and then pileups the raw reads, so "parent is reference"
-- rests on evidence rather than on absence. It also applies child-side quality filters
-- (VAF floor, homopolymer-slip rejection) that this query has no way to express. Expect it
-- to reject the large majority of what you see below — the published germline rate is
-- 40-150 de-novo SNVs per generation, and this query typically returns thousands.
WITH v AS (SELECT * FROM variants WHERE snapshot_id = '${snapshot}')
SELECT c.gene, c.chrom, c.pos, c.ref, c.alt,
       c.clinvar_sig, c.consequence, c.zygosity, c.gnomad_af
FROM v c
WHERE c.sample_id = '${child}'
  AND NOT EXISTS (SELECT 1 FROM v p WHERE p.sample_id = '${father}'
        AND p.chrom=c.chrom AND p.pos=c.pos AND p.ref=c.ref AND p.alt=c.alt)
  AND NOT EXISTS (SELECT 1 FROM v m WHERE m.sample_id = '${mother}'
        AND m.chrom=c.chrom AND m.pos=c.pos AND m.ref=c.ref AND m.alt=c.alt)
ORDER BY (c.clinvar_sig IS NULL), c.gene;
