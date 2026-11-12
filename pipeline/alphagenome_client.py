"""Opt-in AlphaGenome 2 transport and cache. No SDK imports or network at import time.

Only _create_client imports the optional SDK. All remote requests pass through Session;
callers must first resolve EVERY snapshot carrier, including non-panel annotations.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import math
import os
import time
from datetime import datetime, timezone
from pathlib import Path

from . import config, consent, enformer, reportpaths
from .util import log, read_tsv, write_report

SCORERS = ('RNA_SEQ', 'CAGE', 'DNASE', 'CHIP_TF', 'SPLICE_SITES',
           'SPLICE_SITE_USAGE', 'SPLICE_JUNCTIONS')
# POLYADENYLATION returned an AnnData but no tidy rows in the live v2 probe.
LONG_COLS = ('chrom pos ref alt condition scorer gene track raw_score quantile_score epoch '
             'gene_id junction_start junction_end imputation_mode date').split()
EGRESS_COLS = ['timestamp', 'chrom', 'pos', 'ref', 'alt', 'condition', 'scorers', 'payload']
ROOT = Path(__file__).resolve().parents[1]
WINDOW_HALF = 524_288
SV_MAX_BP = 50_000          # largest event size verified against the live service in Phase 0
SV_NEAR_BP = 100_000        # Phase 6: events within this distance of a panel gene are scored
# GRCh38 primary-assembly lengths; unknown/alternate contigs fail closed.
LENGTHS = dict(zip([str(i) for i in range(1, 23)] + ['X','Y'], [
    248956422,242193529,198295559,190214555,181538259,170805979,159345973,
    145138636,138394717,133797422,135086622,133275309,114364328,107043718,
    101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,
    156040895,57227415]))
TERMS = ('Model-predicted candidates only; non-commercial research, not clinical decisions. '
         'Outputs and derived summaries are subject to '
         'https://deepmind.google.com/science/alphagenome/output-terms. '
         'Summary scores are derived absolute effects; splice scores are unvalidated and display-only.')


def now():
    return datetime.now(timezone.utc).isoformat()


def key(v):
    return (str(v['chrom']).removeprefix('chr'), int(v['pos']), v['ref'], v['alt'])


def fmt(v):
    """Human-readable variant, the form the other reports use: `17:43100000 ->A`."""
    chrom, pos, ref, alt = key(v)
    return f'{chrom}:{pos} {ref}>{alt}'


def _atomic(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'{path.name}.{os.getpid()}.part')
    try:
        tmp.write_text(text)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def _write_tsv(path, rows, columns):
    import io
    f = io.StringIO(newline='')
    w = csv.DictWriter(f, columns, delimiter='\t', extrasaction='ignore')
    w.writeheader(); w.writerows(rows)
    _atomic(path, f.getvalue())


def conditions_for(gene):
    rows = [r for r in read_tsv(ROOT / 'panels/gene_tissue.tsv') if r['gene'] == gene]
    return rows or [{'gene': gene, 'condition': config.ALPHAGENOME_DEFAULT_CONDITION,
                     'label': 'blood CD14-positive monocyte (fallback)',
                     'rationale': 'No panel mapping; reference fallback, not disease-specific.'}]


def epoch():
    try:
        return json.loads(config.ALPHAGENOME_EPOCH_FILE.read_text())
    except FileNotFoundError:
        return {}


# The long table reaches ~1,100 rows per scored variant (CHIP_TF alone is ~900), so a family
# run is millions of rows. Parse it once per process, index it by variant, and re-read only
# when the file or the epoch changes. Later rows win for a duplicate key.
_cache = {'sig': None, 'rows': [], 'by_key': {}}


def current_rows():
    path = config.ALPHAGENOME_LONG_FILE
    eid = epoch().get('id')
    try:
        st = path.stat()
    except FileNotFoundError:
        return []
    sig = (str(path), st.st_mtime_ns, st.st_size, eid)
    if _cache['sig'] != sig:
        latest = {}
        if eid:
            for r in read_tsv(path):
                # A crash mid-append can leave one short last line: skip it.
                if r.get('epoch') == eid and r.get('raw_score') not in (None, ''):
                    latest[_row_key(r)] = r
        rows = list(latest.values())
        by_key = {}
        for r in rows:
            by_key.setdefault(key(r), []).append(r)
        _cache.update(sig=sig, rows=rows, by_key=by_key)
    return _cache['rows']


def summary(rows, gene):
    """Maximum absolute value within own-gene RNA tracks; only then DNase fallback."""
    own = [r for r in rows if r['scorer'] == 'RNA_SEQ' and r['gene'] == gene]
    values = own or [r for r in rows if r['scorer'] == 'DNASE']
    finite = [r for r in values if math.isfinite(float(r['raw_score']))]
    if not finite:
        return None
    r = max(finite, key=lambda r: abs(float(r['raw_score'])))
    return abs(float(r['raw_score'])), r['scorer']


def evidence(chrom, pos, ref, alt, gene=None):
    target = (str(chrom).removeprefix('chr'), int(pos), ref, alt)
    current_rows()
    rows = list(_cache['by_key'].get(target, []))
    if gene:
        condition = conditions_for(gene)[0]['condition']
        rows = [r for r in rows if r['condition'] == condition]
    return rows


def splice_label(chrom, pos, ref, alt, gene, consequence):
    from . import spliceai
    if not any(c in (consequence or '') for c in ('splice', 'intron')):
        return '—'
    if spliceai.lookup(chrom, pos, ref, alt):
        return '—'
    rows = evidence(chrom, pos, ref, alt, gene)
    labels = []
    for scorer in SCORERS[-3:]:
        vals = [abs(float(r['raw_score'])) for r in rows if r['scorer'] == scorer
                and r['gene'] in ('', gene) and math.isfinite(float(r['raw_score']))]
        if vals:
            labels.append(f'{scorer} {max(vals):.3g}')
    return ('; '.join(labels) + ' (model prediction, unvalidated; ' +
            conditions_for(gene)[0]['label'] + ')') if labels else '—'


def _create_client(api_key):
    # This is the sole SDK import boundary, replaced by a fake in offline tests.
    try:
        from alphagenome.data import genome
        from alphagenome.models.v2 import dna_client, variant_scorers
    except ImportError:
        raise SystemExit('Install the optional client: uv sync --extra alphagenome') from None
    remote = dna_client.create(api_key, timeout=60)

    class Adapter:
        def score(self, variant, condition):
            chrom, pos, ref, alt = key(variant)
            if 'ref_seq' in variant:
                # Structural variant: full reference and alternate bases, built locally from
                # REF_FASTA (Phase 0: 5 kb and 50 kb DEL/DUP are accepted this way).
                v = genome.Variant('chr' + chrom, pos, variant['ref_seq'], variant['alt_seq'])
                interval = v.reference_interval.resize(1_048_576)
            else:
                v = genome.Variant('chr' + chrom, pos, '' if ref == '-' else ref,
                                   '' if alt == '-' else alt)
                interval = genome.Interval('chr' + chrom, pos - 1 - 524288, pos - 1 + 524288)
            results = remote.score_variant(interval, v, condition=condition,
                variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS[s] for s in SCORERS])
            return variant_scorers.tidy_scores(results).to_dict('records')
    return Adapter()


class Stopped(Exception):
    """Budget or transient service failure; already persisted results remain usable."""


class Session:
    def __init__(self, api_key, max_calls=None):
        self.client = _create_client(api_key)
        self.limit = config.ALPHAGENOME_MAX_CALLS if max_calls is None else max_calls
        self.calls = 0

    def score(self, v, condition):
        for attempt in range(4):
            if self.calls >= self.limit:
                raise Stopped('AlphaGenome call cap reached; rerun to resume.')
            # Persist intent immediately before dispatch, including retry attempts. A crash
            # here can conservatively overcount one attempted request, never conceal egress.
            path = config.ALPHAGENOME_EGRESS_FILE
            path.parent.mkdir(parents=True, exist_ok=True)
            exists = path.exists() and path.stat().st_size > 0
            with path.open('a', newline='') as f:
                w = csv.writer(f, delimiter='\t')
                if not exists:
                    w.writerow(EGRESS_COLS)
                payload = (f'reference sequence {len(v["ref_seq"])} bp + alternate {len(v["alt_seq"])} bp'
                           if 'ref_seq' in v else 'coordinates')
                w.writerow([now(), *key(v), condition, ','.join(SCORERS), payload])
                f.flush(); os.fsync(f.fileno())
            self.calls += 1
            try:
                return self.client.score({k: v[k] for k in ('chrom','pos','ref','alt','ref_seq','alt_seq')
                                          if k in v}, condition)
            except Exception as exc:
                code = str(exc.code()) if hasattr(exc, 'code') else type(exc).__name__
                if code.rsplit('.', 1)[-1] not in ('RESOURCE_EXHAUSTED', 'UNAVAILABLE', 'DEADLINE_EXCEEDED'):
                    # Never print raw RPC errors: metadata may contain credentials.
                    raise Stopped(f'AlphaGenome request failed ({code}); cached results retained.') from None
                if attempt == 3:
                    raise Stopped(f'AlphaGenome unavailable ({code}); rerun to resume.') from None
                time.sleep(2 ** attempt)


def _normalize(v, condition, rows, eid):
    result = []
    def clean(value):
        if value is None or isinstance(value, float) and math.isnan(value):
            return ''
        return str(value)
    for r in rows:
        value = r.get('raw_score')
        if value is None or not math.isfinite(float(value)):
            continue
        scorer = clean(r.get('output_type')).rsplit('.', 1)[-1]
        if scorer not in SCORERS:
            continue
        result.append(dict(zip(('chrom','pos','ref','alt'), key(v))) | {
            'condition': condition, 'scorer': scorer, 'gene': clean(r.get('gene_name')),
            'track': clean(r.get('track_name')), 'raw_score': float(value),
            'quantile_score': '', 'epoch': eid, 'gene_id': clean(r.get('gene_id')),
            'junction_start': clean(r.get('junction_Start')), 'junction_end': clean(r.get('junction_End')),
            'imputation_mode': clean(r.get('imputation_mode')), 'date': now()})
    return result


def _row_key(r):
    return (*key(r), *(r.get(k, '') for k in ('condition','scorer','gene','track','epoch',
                                              'gene_id','junction_start','junction_end')))


def merge(rows):
    """Append scored rows. Rewriting the whole table per variant was quadratic in run size;
    appending is safe because a (variant, condition, epoch) is only scored once, and
    current_rows() collapses any duplicate key to its latest row."""
    path = config.ALPHAGENOME_LONG_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    new = not path.exists() or path.stat().st_size == 0
    torn = False
    if not new:
        with path.open('rb') as f:
            f.seek(-1, os.SEEK_END)
            torn = f.read(1) != b'\n'
    with path.open('a', newline='') as f:
        if torn:
            f.write('\n')   # a crash left a partial line: never glue new rows onto it
        w = csv.DictWriter(f, LONG_COLS, delimiter='\t', extrasaction='ignore')
        if new:
            w.writeheader()
        w.writerows(rows)
        f.flush(); os.fsync(f.fileno())


def _sentinel_values(rows):
    # Compare all finite raw outputs in stable track/gene/junction order, not just a max.
    return {json.dumps(_row_key(r)[4:8] + _row_key(r)[9:]): float(r['raw_score']) for r in rows}


def establish_epoch(session):
    old = epoch()
    values = {}
    for v in read_tsv(ROOT / 'panels/alphagenome_sentinels.tsv'):
        raw = session.score(v, v['condition'])
        normalized = _normalize(v, v['condition'], raw, '')
        if not normalized:
            raise Stopped('Sentinel returned no scores; refusing to reuse an unverified epoch.')
        values[json.dumps(key(v))] = _sentinel_values(normalized)
    previous = old.get('sentinels', {})
    matches = values.keys() == previous.keys() and all(
        values[k].keys() == previous[k].keys() and all(
            math.isclose(v, previous[k][s], rel_tol=1e-6, abs_tol=1e-7)
            for s, v in values[k].items()) for k in values)
    if matches:
        return old['id']
    eid = hashlib.sha256((now() + json.dumps(values, sort_keys=True)).encode()).hexdigest()[:16]
    try:
        version = importlib.metadata.version('alphagenome')
    except importlib.metadata.PackageNotFoundError:
        version = 'test-client'
    # Invalidate summary/calibration before publishing the new epoch. Long rows are retained.
    _write_tsv(config.ENFORMER_SCORES_FILE, [],
               ['chrom','pos','ref','alt','delta_max','l2_center','top_track','top_desc'])
    config.ENFORMER_CALIBRATION_FILE.unlink(missing_ok=True)
    _atomic(config.ALPHAGENOME_EPOCH_FILE, json.dumps({
        'id': eid, 'date': now(), 'package_version': version, 'sentinels': values,
        'previous_epoch': old.get('id'), 'rtol': 1e-6, 'atol': 1e-7}))
    enformer._table = enformer._calibration = None
    if old:
        log.warning('AlphaGenome sentinel change: new epoch %s; prior scores need rescoring.', eid)
    return eid


def eligible(snapshot, variants, *, panel_only=None):
    """Re-query all carriers before egress, independent of rare/consequence/gene filters."""
    import duckdb
    from .panel import load_panel
    panel = {g for d in load_panel() for g in d.genes}
    con = duckdb.connect(str(config.DUCKDB_FILE), read_only=True)
    allowed, withheld = [], 0
    try:
        for v in variants:
            if panel_only and v.get('gene') not in panel:
                continue
            chrom, pos, ref, alt = key(v)
            if (chrom not in LENGTHS or pos - 1 - 524288 < 0 or
                pos - 1 + 524288 > LENGTHS[chrom] or
                any(set(x) - set('ACGTN-') for x in (ref, alt)) or
                max(len(ref), len(alt)) > 50001):
                continue
            carriers = con.execute('SELECT DISTINCT sample_id FROM variants WHERE snapshot_id=? '
                "AND chrom=? AND pos=? AND ref=? AND alt=? AND zygosity IN ('HET','HOM')",
                [snapshot, chrom, pos, ref, alt]).fetchall()
            if not carriers or not all(consent.has_consented(c[0], consent.CATEGORY_EXTERNAL_SCORING)
                                       for c in carriers):
                withheld += 1
                continue
            allowed.append(v)
    finally:
        con.close()
    return allowed, withheld


def score_selected(snapshot, variants, *, dry_run=False, panel_only=None):
    # The calibration null is common variants genome-wide; a panel filter would leave almost
    # none of them, so run_snapshot passes panel_only=False for --calibrate.
    if panel_only is None:
        panel_only = config.ALPHAGENOME_PANEL_ONLY
    allowed, withheld = eligible(snapshot, variants, panel_only=panel_only)
    requests = [(v, c) for v in allowed for c in conditions_for(v.get('gene'))]
    if dry_run:
        print(f'AlphaGenome dry run: {len(allowed)} eligible variants; {withheld} withheld for consent.')
        for v, c in requests:
            print(f'{fmt(v)}\t{c["condition"]}\t{c["label"]}')
        print(f'Up to {len(requests) + (5 if requests else 0)} initial calls including public sentinels; '
              f'cap {config.ALPHAGENOME_MAX_CALLS} counts retries. Cached current-epoch requests are skipped '
              'after live sentinel verification. No channel opened.')
        return allowed, withheld
    api_key = os.environ.get('ALPHAGENOME_API_KEY', '').strip()
    if not api_key:
        raise SystemExit('Set ALPHAGENOME_API_KEY; see docs/runbook.md, AlphaGenome 2.')
    if not requests:
        log.info('AlphaGenome: no eligible consenting carriers; nothing sent.')
        return allowed, withheld
    if config.ALPHAGENOME_MAX_CALLS <= 5:
        raise SystemExit('ALPHAGENOME_MAX_CALLS must exceed the five sentinel calls to make progress.')
    try:
        session = Session(api_key)
        eid = establish_epoch(session)
        existing = current_rows()
        completed = {(*key(r), r['condition']) for r in existing}
        for v, c in requests:
            if (*key(v), c['condition']) in completed:
                continue
            rows = _normalize(v, c['condition'], session.score(v, c['condition']), eid)
            if rows:
                merge(rows)  # Flush every successful variant/condition; never lose a batch.
                completed.add((*key(v), c['condition']))
    except Stopped as exc:
        log.warning('%s', exc)
    rebuild_summary(allowed)
    return allowed, withheld


def rebuild_summary(variants):
    # Include prior current-epoch summaries; new epoch initialization has already emptied them.
    cols = ['chrom','pos','ref','alt','delta_max','l2_center','top_track','top_desc']
    cache = {key(r): r for r in read_tsv(config.ENFORMER_SCORES_FILE)}
    rows = current_rows()
    grouped = {}
    for r in rows:
        grouped.setdefault((*key(r), r['condition']), []).append(r)
    for v in variants:
        c = conditions_for(v.get('gene'))[0]
        result = summary(grouped.get((*key(v), c['condition']), []), v.get('gene'))
        if result:
            cache[key(v)] = dict(zip(cols, [*key(v), result[0], 0, 0, f'{result[1]}:{c["label"]}']))
    _write_tsv(config.ENFORMER_SCORES_FILE, cache.values(), cols)
    enformer._table = None


def run_snapshot(snapshot, *, calibrate=False, dry_run=False, force=False):
    from . import interpret
    if not dry_run and not os.environ.get("ALPHAGENOME_API_KEY", "").strip():
        raise SystemExit("Set ALPHAGENOME_API_KEY; see docs/runbook.md, AlphaGenome 2.")
    variants = interpret._select_common(snapshot) if calibrate else interpret._select_variants(snapshot, set())
    # One writer owns epoch/cache transactions; flock is released even after interruption.
    import fcntl
    config.INTERPRET_DIR.mkdir(parents=True, exist_ok=True)
    with (config.INTERPRET_DIR / 'alphagenome2.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        allowed, withheld = score_selected(snapshot, variants, dry_run=dry_run,
                                           panel_only=False if calibrate else None)
        if dry_run:
            return config.INTERPRET_DIR
        scored = []
        for v in allowed:
            en = enformer.lookup(*key(v))
            if en:
                scored.append({**v, 'delta_max': en[0], 'top_desc': en[3]})
        if calibrate:
            data = {'model': 'alphagenome2', 'epoch': epoch().get('id'),
                    'af_threshold': config.INTERPRET_COMMON_AF}
            for stratum, snv in [('snv',True),('indel',False)]:
                data[stratum] = sorted(v['delta_max'] for v in scored
                    if enformer.is_snv(v['ref'],v['alt']) == snv)
                data['n_'+stratum] = len(data[stratum])
            _atomic(config.ENFORMER_CALIBRATION_FILE, json.dumps(data))
            enformer._calibration = None
            return config.ENFORMER_CALIBRATION_FILE
        dest = reportpaths.family_report('interpret_alphagenome2', snapshot)
        lines = [f'# AlphaGenome 2 model predictions — {snapshot}', '', TERMS, '',
                 f'External API; report date {now()}; epoch {epoch().get("id", "none")}.',
                 f'{len(scored)} scored; {withheld} withheld for lack of carrier consent.',
                 '1,048,576 bp GRCh38 context; absolute own-gene RNA log fold change '
                 '(maximum across same-condition RNA tracks), or DNase fallback. Scorer units differ.',
                 'Uncalibrated values do not establish a regulatory hit. Run interpret --calibrate.', '',
                 '| Variant | Gene | Model prediction / tissue | AlphaGenome splice (unvalidated) |',
                 '|---|---|---|---|']
        for v in sorted(scored, key=lambda v: -v['delta_max']):
            en = enformer.lookup(*key(v))
            lines.append(f'| {fmt(v)} | {v.get("gene", "")} | {enformer.label(en,v["ref"],v["alt"])} | '
                         f'{splice_label(*key(v),v.get("gene"),v.get("consequence"))} |')
        write_report(dest, lines)
        return dest


# --- Phase 6: structural variants ---------------------------------------------------------

def _sv_key(sv):
    """Long-table key for an SV: symbolic alleles, so no sequence is stored in the cache."""
    t = sv['svtype']
    return {'chrom': str(sv['chrom']), 'pos': int(sv['start']),
            'ref': f'<{t}>', 'alt': f'<{t}:{int(sv["end"])}>'}


def sv_scorable(sv):
    """DEL/DUP no larger than the size verified in Phase 0, whose 1 Mb window fits the contig."""
    chrom = str(sv['chrom']).removeprefix('chr')
    if sv['svtype'] not in ('DEL', 'DUP') or chrom not in LENGTHS:
        return False
    if int(sv['end']) - int(sv['start']) + 1 > SV_MAX_BP:
        return False
    mid = (int(sv['start']) + int(sv['end'])) // 2
    return mid - WINDOW_HALF >= 0 and mid + WINDOW_HALF <= LENGTHS[chrom]


def _sv_sequences(sv, fetch):
    """(ref, alt) bases from the local reference. Manta POS is the anchor base before the
    event and END its last base, so POS..END covers anchor + affected bases."""
    seq = fetch(str(sv['chrom']), int(sv['start']), int(sv['end']))
    if len(seq) < 2 or set(seq) - set('ACGTN'):
        return None
    return (seq, seq[0]) if sv['svtype'] == 'DEL' else (seq, seq + seq[1:])


def sv_consent_ok(sample):
    """SVs are not in DuckDB, so carriers among relatives cannot be enumerated. Fail closed:
    every sequenced person (and this sample) must have granted external_scoring."""
    from .pedigree import sequenced_ids
    people = sequenced_ids() | {sample}
    return all(consent.has_consented(s, consent.CATEGORY_EXTERNAL_SCORING) for s in people)


def _sv_label(rows, gene, label):
    own = [r for r in rows if r['scorer'] == 'RNA_SEQ' and r['gene'] == gene
           and math.isfinite(float(r['raw_score']))]
    if not own:
        return None
    best = max(own, key=lambda r: abs(float(r['raw_score'])))
    return f'{float(best["raw_score"]):+.3f} RNA log-fold change ({label}; model prediction)'


def sv_predictions(sample, pairs, *, fetch=None):
    """pairs: [(sv, gene_name)]. Returns ({(chrom,start,end,svtype,gene): label}, note).

    Scores missing pairs only when a key is set and consent holds; otherwise shows whatever
    the current epoch already has cached. The request carries reference sequence, not just
    coordinates; the egress log records that."""
    if fetch is None:
        from .interpret import _faidx as fetch
    pairs = [(sv, g) for sv, g in pairs if sv_scorable(sv)]
    if not pairs:
        return {}, 'no deletion/duplication within the scorable size near a panel gene'
    note = ''
    api_key = os.environ.get('ALPHAGENOME_API_KEY', '').strip()
    if not api_key:
        note = 'ALPHAGENOME_API_KEY not set; showing cached predictions only'
    elif not sv_consent_ok(sample):
        note = ('not every sequenced family member has granted external_scoring; '
                'nothing sent, cached predictions only')
    else:
        import fcntl
        config.INTERPRET_DIR.mkdir(parents=True, exist_ok=True)
        with (config.INTERPRET_DIR / 'alphagenome2.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            try:
                session = Session(api_key)
                eid = establish_epoch(session)
                done = {(*key(r), r['condition']) for r in current_rows()}
                for sv, gene in pairs:
                    k = _sv_key(sv)
                    cond = conditions_for(gene)[0]['condition']
                    if (*key(k), cond) in done:
                        continue
                    seqs = _sv_sequences(sv, fetch)
                    if not seqs:
                        continue
                    v = {**k, 'ref_seq': seqs[0], 'alt_seq': seqs[1]}
                    rows = _normalize(k, cond, session.score(v, cond), eid)
                    if rows:
                        merge(rows)
                        done.add((*key(k), cond))
            except Stopped as exc:
                log.warning('%s', exc)
                note = str(exc)
    grouped = {}
    for r in current_rows():
        grouped.setdefault((*key(r), r['condition']), []).append(r)
    out = {}
    for sv, gene in pairs:
        c = conditions_for(gene)[0]
        lab = _sv_label(grouped.get((*key(_sv_key(sv)), c['condition']), []), gene, c['label'])
        if lab:
            out[(str(sv['chrom']), int(sv['start']), int(sv['end']), sv['svtype'], gene)] = lab
    return out, note
