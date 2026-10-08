"""Phase 4: incremental scan — skip re-annotation when nothing a sample's rows depend on changed.

A sample's `variants` rows for a snapshot are a deterministic function of exactly three inputs:
the normalized VCF (the genome), the snapshot's ClinVar VCF (attached `type=exact`), and the VEP
cache version. The model lookups (AlphaMissense/SpliceAI/Enformer/constraint) are joined at
report time, never stored — so they don't enter this decision. When all three digests match what
produced the previous snapshot's rows, the new snapshot's rows are *byte-identical*: we copy them
forward (a cheap re-stamp) instead of re-running VEP.

A `scan_inputs` ledger records the digests per (sample, snapshot). The reproducibility ground
truth is preserved by a **periodic full anchor**: the first scan of each calendar month, any VEP
cache-version bump, and an explicit `--full` all force a from-scratch re-annotation. The `validate`
harness is the drift detector between anchors.

Two incremental paths (see `incremental_kind`):
  * copy_forward    — all three inputs unchanged → re-stamp prev's rows (no VEP).
  * clinvar_delta   — genome + cache unchanged, only ClinVar advanced (the normal weekly case)
    → copy everything forward, then re-run VEP on ONLY the carried variants at ClinVar-changed
    coordinates and merge them in. Correct because ClinVar is `type=exact`: a variant's ClinVar
    columns can only change if the ClinVar record at its exact key changed, and the Δ set is
    computed by parsing both ClinVar VCFs identically (so change-detection can't miss a change).
    The functional columns come from the byte-identical copy-forward; the ClinVar columns from a
    real VEP pass over the tiny subset → identical to a full run.
"""

from __future__ import annotations

from pathlib import Path

from . import config
from .annotate import annotate_sites, clinvar_vcf
from .load import TABLE, connect, load_subset, transaction
from .normalize import normalized_path
from .snapshot import load_manifest
from .util import log, run, sha256

# Inputs that determine a sample's variants rows (everything else is joined at report time).
_INPUT_SOURCES = ("normalized_vcf", "clinvar", "vep_cache")
_ANCHOR_SAMPLE = "__anchor__"  # sentinel row in scan_inputs marking a full-anchor snapshot

# The full variants column list, for the copy-forward re-stamp (snapshot_id is replaced).
_COLS = ("sample_id, chrom, pos, ref, alt, gene, consequence, clinvar_sig, clinvar_revstat, "
         "clinvar_alleleid, clinvar_disease, gnomad_af, sift, polyphen, zygosity")


# Bumped when the digest SCHEME changes, so old-format ledger rows compare unequal and the
# affected samples re-annotate once rather than being silently trusted under a new meaning.
_DIGEST_SCHEME = "sha256"


def source_digests(sample: str, snapshot_id: str) -> dict[str, str]:
    """Content digests for the three inputs that determine a sample's rows.

    The normalized VCF is hashed, not fingerprinted by (size, mtime). The old proxy argued
    that a content change at identical size and whole-second mtime "is not physically
    realizable" — true of an edit in place, but not of the operations this pipeline actually
    performs on these files. `rsync -a` and `cp -p` preserve mtime by design, so a file
    restored from backup, copied between hosts, or re-staged from the NAS presents the
    previous fingerprint with different bytes; a server migration or restore may copy files
    with timestamps preserved for exactly that reason. Being wrong here costs either a missed
    re-annotation or a needless one, and hashing measured 0.7 s against a ~90 min
    scan — <0.01% — so there is no reason to approximate.

    ClinVar and the VEP cache come from the pinned manifest, so their release identifiers are
    already exact.
    """
    meta = load_manifest(snapshot_id)
    return {
        "normalized_vcf": f"{_DIGEST_SCHEME}:{sha256(normalized_path(sample))}",
        "clinvar": meta.get("clinvar_date", "?"),
        "vep_cache": meta.get("vep_cache_version", "?"),
    }


def _ledger(con, sample: str, snapshot_id: str) -> dict[str, str]:
    return {s: d for s, d in con.execute(
        "SELECT source, digest FROM scan_inputs WHERE sample_id = ? AND snapshot_id = ?",
        [sample, snapshot_id]).fetchall()}


def record_inputs(sample: str, snapshot_id: str, digests: dict[str, str],
                  rows_loaded: int) -> None:
    con = connect()
    with transaction(con):   # a half-written ledger row makes the next scan misjudge the kind
        con.execute("DELETE FROM scan_inputs WHERE sample_id = ? AND snapshot_id = ?",
                    [sample, snapshot_id])
        con.executemany(
            "INSERT INTO scan_inputs (sample_id, snapshot_id, source, digest, rows_loaded) "
            "VALUES (?,?,?,?,?)",
            [[sample, snapshot_id, src, digests[src], rows_loaded] for src in _INPUT_SOURCES])
    con.close()


def incremental_kind(sample: str, prev: str | None, curr: str) -> str:
    """Per-sample decision within an incremental scan: 'copy_forward' (nothing changed),
    'clinvar_delta' (only ClinVar advanced — re-annotate just the affected carried variants),
    or 'full' (genome or VEP cache changed, or no comparable prev → re-annotate everything)."""
    if prev is None:
        return "full"
    con = connect()
    prior = _ledger(con, sample, prev)
    con.close()
    if not prior:
        return "full"
    now = source_digests(sample, curr)
    if all(prior.get(k) == now.get(k) for k in _INPUT_SOURCES):
        return "copy_forward"
    if (config.SCAN_CLINVAR_DELTA
            and prior.get("normalized_vcf") == now.get("normalized_vcf")
            and prior.get("vep_cache") == now.get("vep_cache")):
        return "clinvar_delta"  # genome + cache unchanged, only ClinVar differs
    return "full"


def _clinvar_select(path: Path) -> str:
    """Parse a ClinVar VCF's (chrom,pos,ref,alt) + the four CLINVAR INFO fields. Used identically
    on both snapshots' ClinVar, so change-detection is independent of how VEP formats the values
    — a difference here is a real ClinVar change."""
    return f"""
        SELECT chrom, pos, ref, alt,
               regexp_extract(info, 'CLNSIG=([^;]*)', 1)      AS clnsig,
               regexp_extract(info, 'CLNREVSTAT=([^;]*)', 1)  AS clnrevstat,
               regexp_extract(info, 'CLNDN=([^;]*)', 1)       AS clndn,
               regexp_extract(info, 'ALLELEID=([^;]*)', 1)    AS alleleid
        FROM read_csv('{path}', delim='\\t', header=false, comment='#', compression='gzip',
                      auto_detect=false,
                      columns={{'chrom':'VARCHAR','pos':'BIGINT','id':'VARCHAR','ref':'VARCHAR',
                                'alt':'VARCHAR','qual':'VARCHAR','filter':'VARCHAR',
                                'info':'VARCHAR'}})"""


# Memo of the ClinVar-level Δ coordinate set per (prev, curr) — sample-independent, so it is
# computed once per scan rather than re-parsing ClinVar for every family member.
_DELTA_CACHE: dict[tuple[str, str], list[tuple]] = {}


def clinvar_delta_coords(prev: str, curr: str) -> list[tuple]:
    """(chrom,pos,ref,alt) where ClinVar's CLNSIG/CLNREVSTAT/CLNDN/ALLELEID changed (incl.
    records added or removed) between the two snapshots. Conservative: any difference is
    included, so no real change is ever missed."""
    key = (prev, curr)
    if key in _DELTA_CACHE:
        return _DELTA_CACHE[key]
    con = connect()
    rows = con.execute(f"""
        WITH cp AS ({_clinvar_select(clinvar_vcf(prev))}),
             cc AS ({_clinvar_select(clinvar_vcf(curr))})
        SELECT coalesce(cc.chrom, cp.chrom), coalesce(cc.pos, cp.pos),
               coalesce(cc.ref, cp.ref),     coalesce(cc.alt, cp.alt)
        FROM cc FULL OUTER JOIN cp USING (chrom, pos, ref, alt)
        WHERE cc.clnsig     IS DISTINCT FROM cp.clnsig
           OR cc.clnrevstat IS DISTINCT FROM cp.clnrevstat
           OR cc.clndn      IS DISTINCT FROM cp.clndn
           OR cc.alleleid   IS DISTINCT FROM cp.alleleid
    """).fetchall()
    con.close()
    _DELTA_CACHE[key] = rows
    log.info("ClinVar Δ %s → %s: %d changed coordinate(s).", prev, curr, len(rows))
    return rows


def delta_match_rows(coords: list[tuple]) -> list[list]:
    """Expand ClinVar Δ coordinates into (chrom, pos, ref, alt) match rows for the join below.

    The two sides speak different allele languages. ClinVar is a VCF: an indel carries an
    anchor base and sits at the anchor's coordinate ("1 100 G GA"). The variants table holds
    VEP's Uploaded_variation representation, which drops the anchor and shifts one base right
    ("1_101_-/A"). An exact four-column join therefore matches SNVs and matches *no* indel at
    all — about a fifth of a family database's rows are dash-notation, against a few thousand
    in anchored form, and indels are 31% of ClinVar's own pathogenic/likely-pathogenic records.
    Every incremental scan silently skipped that third of the delta; the monthly full anchor
    was the only thing that ever caught up, so a P/LP reclassification on a frameshift could
    sit unreported for up to five weekly scans.

    Same-length alleles (SNV/MNV) keep the exact allele match — VEP leaves those in place.
    Indels degrade to a position match at both the anchored and the shifted coordinate, with
    the allele test dropped. Over-matching is the safe direction here: a surplus coordinate
    costs one extra variant in the VEP subset and re-annotation is idempotent, whereas a miss
    is a ClinVar change that never reaches the report.
    """
    rows: list[list] = []
    for chrom, pos, ref, alt in coords:
        if ref and alt and len(ref) == len(alt):
            rows.append([chrom, pos, ref, alt])
        else:
            rows.append([chrom, pos, None, None])
            rows.append([chrom, pos + 1, None, None])
    return rows


def _carried_changed(sample: str, prev: str, coords: list[tuple]) -> list[tuple]:
    """Intersect the ClinVar Δ with the variants this sample actually carries."""
    match = delta_match_rows(coords)
    n_indel = sum(1 for _c, _p, ref, _a in match if ref is None) // 2
    if n_indel:
        log.info("ClinVar Δ: %d indel coordinate(s) matched by position (anchor±1) — VEP and "
                 "ClinVar disagree on indel representation.", n_indel)
    con = connect()
    con.execute("CREATE TEMP TABLE _delta(chrom VARCHAR, pos BIGINT, ref VARCHAR, alt VARCHAR)")
    con.executemany("INSERT INTO _delta VALUES (?,?,?,?)", match)
    # Equi-join on (chrom,pos) with the allele test as a residual predicate, rather than an
    # OR across four columns: this stays a hash join against a 100M-row table.
    rows = con.execute(
        f"SELECT DISTINCT v.chrom, v.pos, v.ref, v.alt FROM {TABLE} v "
        f"JOIN _delta d ON v.chrom = d.chrom AND v.pos = d.pos "
        f"WHERE (d.ref IS NULL OR (v.ref = d.ref AND v.alt = d.alt)) "
        f"AND v.snapshot_id = ? AND v.sample_id = ?", [prev, sample]).fetchall()
    con.close()
    return rows


def site_regions(keys: list[tuple]) -> list[tuple]:
    """(chrom, pos) lookups in the sample's VCF for variants keyed the way the database
    keys them.

    The other half of the representation gap `delta_match_rows` crosses, in the opposite
    direction. `keys` are VEP coordinates; the normalized VCF is anchored. An insertion the
    database holds as "1:101 -/A" is the VCF record "1 100 G GA", which occupies position
    100 only — a lookup at 101 returns nothing, so the insertion was matched against the
    ClinVar Δ and then silently left out of the re-annotation, keeping its copied-forward
    classification until the next full anchor. A deletion survived only because its
    anchored record happens to span the shifted position.

    Every unanchored key therefore also asks for the base before it. Over-asking is safe for
    the same reason it is in `delta_match_rows`: a neighbouring record that comes along is
    re-annotated to the row it already has.
    """
    out = set()
    for chrom, pos, ref, alt in keys:
        out.add((chrom, pos))
        if ref == "-" or alt == "-":
            out.add((chrom, pos - 1))
    return sorted(out)


def _extract_sites(sample: str, snapshot: str, keys: list[tuple]) -> Path:
    """Pull the carried, ClinVar-changed variants out of the sample's (indexed) normalized VCF
    into a small VCF, by region — fast via the tabix index."""
    work = config.ANNOTATED_DIR / snapshot
    work.mkdir(parents=True, exist_ok=True)
    regions = work / f"{sample}.clinvar_delta.regions.tsv"
    regions.write_text("".join(f"{c}\t{p}\n" for c, p in site_regions(keys)))
    out = work / f"{sample}.clinvar_delta.vcf.gz"
    run(["bcftools", "view", "-R", str(regions), str(normalized_path(sample)),
         "-Oz", "-o", str(out)])
    run(["tabix", "-f", "-p", "vcf", str(out)])
    return out


def reannotate_clinvar(sample: str, prev: str, curr: str) -> int:
    """Variant-level incremental: copy everything forward, then re-annotate (via VEP, so the
    rows are byte-identical to a full run) ONLY the carried variants whose ClinVar entry moved."""
    n = copy_forward(sample, prev, curr)
    coords = clinvar_delta_coords(prev, curr)
    keys = _carried_changed(sample, prev, coords) if coords else []
    if not keys:
        log.info("clinvar-delta %s: no carried variant touched by the ClinVar update — "
                 "pure copy-forward.", sample)
        return n
    sites = _extract_sites(sample, curr, keys)
    tsv = annotate_sites(sample, curr, sites)
    load_subset(sample, curr, tsv)
    log.info("clinvar-delta %s: %d carried variant(s) at ClinVar-changed coords re-annotated "
             "(rest copied forward).", sample, len(keys))
    return n


def copy_forward(sample: str, prev: str, curr: str) -> int:
    """Re-stamp the sample's prev rows into curr (idempotent). Byte-identical to re-annotation
    because all three inputs were verified unchanged by `incremental_kind()`."""
    # DELETE-then-INSERT-SELECT is only a re-stamp while the two snapshots are distinct. If
    # prev == curr the DELETE removes the very rows the INSERT would have read, emptying the
    # sample (see `snapshot.previous_snapshot` — the caller must never pass a predecessor
    # that isn't strictly older). Refuse rather than silently return 0 rows.
    if prev >= curr:
        raise SystemExit(
            f"copy_forward({sample}): predecessor {prev!r} is not strictly older than "
            f"{curr!r} — refusing (this would delete the snapshot's rows).")
    con = connect()
    # Atomic, like `load`: a failure between the DELETE and the INSERT would leave the sample
    # emptied for `curr` — worse than not having copied it forward at all, and the scan that
    # was interrupted is exactly when it happens.
    with transaction(con):
        con.execute(f"DELETE FROM {TABLE} WHERE snapshot_id = ? AND sample_id = ?",
                    [curr, sample])
        con.execute(f"INSERT INTO {TABLE} SELECT ?, {_COLS} FROM {TABLE} "
                    f"WHERE snapshot_id = ? AND sample_id = ?", [curr, prev, sample])
    n = con.execute(f"SELECT count(*) FROM {TABLE} WHERE snapshot_id = ? AND sample_id = ?",
                    [curr, sample]).fetchone()[0]
    con.close()
    log.info("copy-forward %s: %d rows re-stamped %s → %s (inputs unchanged).",
             sample, n, prev, curr)
    return n


def _anchor_exists_for_month(con, month: str) -> bool:
    return con.execute(
        "SELECT count(*) FROM scan_inputs WHERE sample_id = ? AND snapshot_id LIKE ?",
        [_ANCHOR_SAMPLE, f"{month}%"]).fetchone()[0] > 0


def mark_anchor(snapshot_id: str) -> None:
    """Record that this snapshot was a full re-annotation — the month's reproducibility anchor."""
    con = connect()
    with transaction(con):   # losing the anchor row silently downgrades the month's guarantee
        con.execute("DELETE FROM scan_inputs WHERE sample_id = ? AND snapshot_id = ?",
                    [_ANCHOR_SAMPLE, snapshot_id])
        con.execute("INSERT INTO scan_inputs (sample_id, snapshot_id, source, digest) "
                    "VALUES (?,?,?,?)", [_ANCHOR_SAMPLE, snapshot_id, "full", snapshot_id])
    con.close()


def scan_mode(prev: str | None, curr: str, force_full: bool) -> tuple[bool, str]:
    """Decide whether this scan is a full re-annotation. Returns (is_full, reason)."""
    if force_full:
        return True, "forced (--full)"
    if prev is None:
        return True, "baseline (no previous snapshot)"
    pm, cm = load_manifest(prev), load_manifest(curr)
    if pm.get("vep_cache_version") != cm.get("vep_cache_version"):
        return True, (f"VEP cache bump {pm.get('vep_cache_version')} → "
                      f"{cm.get('vep_cache_version')}")
    con = connect()
    has_anchor = _anchor_exists_for_month(con, curr[:7])
    con.close()
    if not has_anchor:
        return True, f"monthly full anchor ({curr[:7]})"
    return False, "incremental (copy-forward where inputs unchanged)"
