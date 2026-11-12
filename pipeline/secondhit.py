"""Research-only recessive second-variant candidates; never changes ACMG."""
from __future__ import annotations

from . import alphagenome_client as ag, clinsig, config, enformer, pedigree, reportpaths
from .load import connect
from .panel import load_panel
from .util import write_report


def select(snapshot, cap=50):
    con = connect()
    ar = {g for d in load_panel() if d.inh == 'AR' for g in d.genes}
    rows = con.execute("SELECT DISTINCT sample_id,gene,chrom,pos,ref,alt,consequence,clinvar_sig,gnomad_af "
                       "FROM variants WHERE snapshot_id=? AND zygosity='HET'", [snapshot]).fetchall()
    con.close()
    grouped = {}
    for sample,gene,chrom,pos,ref,alt,csq,sig,af in rows:
        if gene in ar:
            grouped.setdefault((sample,gene), {})[(chrom,pos,ref,alt)] = {
                'sample':sample,'gene':gene,'chrom':chrom,'pos':pos,'ref':ref,'alt':alt,
                'consequence':csq,'sig':sig,'af':af}
    candidates = []
    for vs in grouped.values():
        hits = [v for v in vs.values() if clinsig.is_actionable(v['sig'])]
        if len(hits) != 1:
            continue
        known = hits[0]
        others = [v for v in vs.values() if ag.key(v) != ag.key(known)
                  and (v['af'] is None or v['af'] <= .01)]
        others.sort(key=lambda v: (v['af'] is not None,v['af'] or 0,ag.key(v)))
        candidates += [{**v,'known':known} for v in others[:cap]]
    return candidates


def parental_phase(snapshot, candidate):
    """Variant-only parental origin is suggestive, never read-confirmed phase."""
    fam = pedigree.load(sequenced=pedigree.sequenced_ids())
    trio = fam.trio(candidate['sample'])
    if not trio:
        return 'unknown'
    con = connect()
    def origin(v):
        carriers = {fam.person_of(r[0]) for r in con.execute(
            "SELECT DISTINCT sample_id FROM variants WHERE snapshot_id=? AND chrom=? AND pos=? "
            "AND ref=? AND alt=? AND zygosity IN ('HET','HOM')", [snapshot,*ag.key(v)]).fetchall()}
        parents = set(trio[1:]) & carriers
        return next(iter(parents)) if len(parents) == 1 else None
    try:
        a,b=origin(candidate),origin(candidate['known'])
        return ('trans (parental-origin candidate)' if a != b else 'cis (parental-origin candidate)') if a and b else 'unknown'
    finally:
        con.close()


def second_hit(snapshot, *, dry_run=False, phase=False):
    if config.INTERPRET_MODEL != 'alphagenome2':
        raise SystemExit('second-hit requires INTERPRET_MODEL=alphagenome2')
    candidates = select(snapshot)
    unique = {(*ag.key(v),v['gene']):v for v in candidates}
    import fcntl
    config.INTERPRET_DIR.mkdir(parents=True,exist_ok=True)
    with (config.INTERPRET_DIR/'alphagenome2.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        _,withheld = ag.score_selected(snapshot,list(unique.values()),dry_run=dry_run)
    if dry_run:
        return config.INTERPRET_DIR
    for v in candidates:
        v['phase'] = parental_phase(snapshot,v)
        if phase:
            from .phasing import phase_pair
            v['phase'] = phase_pair(v['sample'],v['gene'],v['known'],v)
        rows=ag.evidence(*ag.key(v),v['gene'])
        expr=ag.summary([r for r in rows if r['scorer']=='RNA_SEQ'],v['gene'])
        v['expression']=expr[0] if expr else None
        v['splice']=ag.splice_label(*ag.key(v),v['gene'],v['consequence'])
    candidates.sort(key=lambda v:(v['phase'].startswith('cis'),v['expression'] is None,-(v['expression'] or 0)))
    lines=[f'# Recessive second-variant candidates — {snapshot}', '', ag.TERMS, '',
           'Candidate second variant, model-predicted, unconfirmed. Phenotype compatibility is not '
           'recorded in this pipeline and requires independent review. Absence of a candidate does '
           'not exclude a second hit. Up to 50 rare candidates per gene and sample are considered.',
           'Parental-origin phase is provisional: missing variant calls do not prove reference. '
           'Use --phase for read-backed pairs. Expression and splice evidence are reported separately.',
           f'{withheld} variants withheld for carrier consent.', '',
           '| Sample | Gene / tissue | Known variant | Candidate | Phase | Absolute expression prediction | Splice (unvalidated) |',
           '|---|---|---|---|---|---|---|']
    for v in candidates:
        lines.append(f'| {v["sample"]} | {v["gene"]} / {ag.conditions_for(v["gene"])[0]["label"]} | '
                     f'{ag.fmt(v["known"])} | {ag.fmt(v)} | {v["phase"]} | '
                     f'{v["expression"] if v["expression"] is not None else "unscored"} | {v["splice"]} |')
    dest=reportpaths.family_report('second_hit',snapshot)
    write_report(dest,lines)
    return dest
