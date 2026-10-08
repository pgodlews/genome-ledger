"""Phase 2: provenance/confidence metadata + per-snapshot variant history.

Two derived tables, populated from the same bounded "clinically-relevant" slice of `variants`
(anything ClinVar-classified or on a panel gene — thousands of rows, not the ~4M per sample):

  * `annotation_provenance` — for each (variant, source) used in a report/triage: the source's
    pinned version, a compact value summary, a qualitative confidence band, and whether the
    value fed an ACMG evidence code. Makes every number in a report traceable.
  * `variant_history` — one row per variant per snapshot (clinvar sig/review/AF, plus acmg tier
    and model band when available). The cross-snapshot timeline the triage *stability* factor
    reads, and the substrate the validation harness replays over.

Confidence bands are a single deterministic mapping (auditable, adjustable). Both writers are
idempotent: provenance is delete-then-insert per snapshot, history is upserted on its PK.
"""

from __future__ import annotations

from . import (acmg, alphamissense, config, enformer, gnomad_constraint,
               spliceai, triage)
from .load import TABLE, connect, transaction
from .panel import sql_in_list
from .snapshot import load_manifest
from .util import log, review_stars, today

# Sources whose value participates in an ACMG code, by the codes they can drive.
_ACMG_PVS1 = {"PVS1", "PVS1_moderate"}
_ACMG_FREQ = {"BA1", "BS1", "PM2"}
_ACMG_INSILICO = {"PP3", "BP4"}


def _conf_stars(revstat: str | None) -> str:
    s = review_stars(revstat)
    return "High" if s >= 3 else "Medium" if s == 2 else "Low" if s == 1 else "Very Low"


def _conf_enformer(en, ref: str, alt: str) -> str:
    if not en:
        return "Very Low"
    if enformer.calibration_available():
        return "High" if enformer.band(en[0], ref, alt) == "high" else "Medium"
    return "Low"  # raw Δ fallback band — uncalibrated, so weaker confidence


def _relevant(con, snapshot_id: str, panel_genes: set[str]):
    """Distinct clinically-relevant variants for a snapshot (annotation-level — deduped across
    samples, since the annotation is identical per snapshot for a given key)."""
    gl = sql_in_list(panel_genes)
    # Prefer a PANEL gene for the representative row: any_value(gene) could pick a
    # non-panel overlapping gene, after which downstream `gene in panel_genes` checks
    # silently skip ACMG tier/provenance for a variant that IS on-panel.
    return con.execute(f"""
        SELECT chrom, pos, ref, alt,
               coalesce(max(CASE WHEN gene IN ({gl}) THEN gene END),
                        any_value(gene)) AS gene, any_value(consequence) AS csq,
               any_value(clinvar_sig) AS sig, any_value(clinvar_revstat) AS rev,
               any_value(gnomad_af) AS af, any_value(sift) AS sift,
               any_value(polyphen) AS pp
        FROM {TABLE}
        WHERE snapshot_id = ? AND (clinvar_sig IS NOT NULL OR gene IN ({gl}))
        GROUP BY chrom, pos, ref, alt
    """, [snapshot_id]).fetchall()


def record_history(snapshot_id: str) -> int:
    """Upsert one variant_history row per clinically-relevant variant in the snapshot."""
    con = connect()
    panel_genes, gene_modes = acmg.panel_context()
    have_am, have_sa = alphamissense.available(), spliceai.available()
    have_en, have_gc = enformer.available(), gnomad_constraint.available()
    try:
        sdate = load_manifest(snapshot_id).get("created")
    except Exception:  # noqa: BLE001
        sdate = None

    out = []
    for chrom, pos, ref, alt, gene, csq, sig, rev, af, sift, pp in _relevant(
            con, snapshot_id, panel_genes):
        tier = None
        if gene in panel_genes:
            call = acmg.evaluate(gene, chrom, pos, ref, alt, csq, sig, af, sift, pp, None,
                                 disease_genes=panel_genes, gene_modes=gene_modes,
                                 have_am=have_am, have_sa=have_sa, have_gc=have_gc)
            tier = call["tier"] if call else None
        _am, _sa, en = triage.model_lookups(chrom, pos, ref, alt, csq,
                                             have_am, have_sa, have_en)
        band = enformer.band(en[0], ref, alt) if en else None
        out.append([chrom, pos, ref, alt, snapshot_id, sdate, sig, rev, af, tier, band])

    if out:   # DuckDB's executemany rejects an empty parameter list outright, and a snapshot
        # with no history-relevant variants is a legitimate outcome, not an error.
        con.executemany(
            "INSERT OR REPLACE INTO variant_history (chrom, pos, ref, alt, snapshot_id, "
            "snapshot_date, clinvar_sig, clinvar_revstat, gnomad_af, acmg_tier, model_band) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?)", out)
    con.close()
    log.info("variant_history: %d variants recorded for %s.", len(out), snapshot_id)
    return len(out)


def record_provenance(snapshot_id: str) -> int:
    """Replace this snapshot's annotation_provenance rows: one per (variant, source) used."""
    con = connect()
    panel_genes, gene_modes = acmg.panel_context()
    have_am, have_sa = alphamissense.available(), spliceai.available()
    have_en, have_gc = enformer.available(), gnomad_constraint.available()
    try:
        meta = load_manifest(snapshot_id)
    except Exception:  # noqa: BLE001
        meta = {}
    clinvar_ver = f"ClinVar {meta.get('clinvar_date', '?')}"
    vep_ver = f"VEP r{meta.get('vep_cache_version', '?')}"
    gnomad_ver = meta.get("gnomad_version", "?")
    as_of = today()

    rows = []  # (source, version, value_summary, confidence, used_in_acmg)
    for chrom, pos, ref, alt, gene, csq, sig, rev, af, sift, pp in _relevant(
            con, snapshot_id, panel_genes):
        call = None
        if gene in panel_genes:
            call = acmg.evaluate(gene, chrom, pos, ref, alt, csq, sig, af, sift, pp, None,
                                 disease_genes=panel_genes, gene_modes=gene_modes,
                                 have_am=have_am, have_sa=have_sa, have_gc=have_gc)
        codes = set(call["codes"]) if call and not call.get("deferred") else set()

        def add(src, ver, val, conf, used):
            rows.append([snapshot_id, None, chrom, pos, ref, alt, src, ver, val, conf,
                         used, as_of])

        # ClinVar is deliberately held out of the ACMG criteria → never "used_in_acmg".
        if sig:
            add("clinvar", clinvar_ver, sig, _conf_stars(rev), False)
        if csq:
            add("vep_csq", vep_ver, csq, "High", bool(codes & _ACMG_PVS1))
        if af is not None:
            add("gnomad_af", gnomad_ver, f"{af:.3g}", "High", bool(codes & _ACMG_FREQ))
        # Model + constraint provenance only where they were actually consulted (panel-gene
        # coding variants); that's also the only place used_in_acmg is meaningful.
        if call and not call.get("deferred"):
            am, sa, en = triage.model_lookups(chrom, pos, ref, alt, csq,
                                              have_am, have_sa, have_en)
            if am:
                add("alphamissense", "AlphaMissense hg38", f"{am[0]:.3f} {am[1]}",
                    "Medium", bool(codes & _ACMG_INSILICO))
            if sa:
                add("spliceai", "SpliceAI masked",
                    f"Δ{sa[0]:.2f} {sa[1]}", "Medium", bool(codes & _ACMG_INSILICO))
            if en:
                add("enformer", f"{config.INTERPRET_MODEL}",
                    enformer.label(en, ref, alt), _conf_enformer(en, ref, alt), False)
            if have_gc and gnomad_constraint.lookup(gene):
                add("gnomad_constraint", f"gnomAD {config.GNOMAD_CONSTRAINT_VERSION}",
                    gnomad_constraint.label(gene), "Medium", bool(codes & _ACMG_PVS1))

    # Atomic: a failed INSERT would otherwise leave the snapshot with no provenance rows,
    # which is indistinguishable from "nothing was annotated" when the report reads it back.
    with transaction(con):
        con.execute("DELETE FROM annotation_provenance WHERE snapshot_id = ?", [snapshot_id])
        if rows:   # DuckDB's executemany rejects an empty parameter list outright
            con.executemany(
                "INSERT INTO annotation_provenance (snapshot_id, sample_id, chrom, pos, ref, "
                "alt, source, version, value_summary, confidence, used_in_acmg, as_of) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", rows)
    con.close()
    log.info("annotation_provenance: %d rows for %s.", len(rows), snapshot_id)
    return len(rows)


def record(snapshot_id: str) -> None:
    """Both Phase-2 writers for a snapshot (the scan hook)."""
    record_history(snapshot_id)
    record_provenance(snapshot_id)


def backfill_history() -> None:
    """One-pass population of variant_history over every snapshot already loaded in DuckDB —
    so the stability factor has depth on day one rather than after N future scans."""
    con = connect()
    snaps = [r[0] for r in con.execute(
        f"SELECT DISTINCT snapshot_id FROM {TABLE} ORDER BY snapshot_id").fetchall()]
    con.close()
    log.info("Backfilling variant_history over %d snapshot(s): %s",
             len(snaps), ", ".join(snaps))
    for s in snaps:
        record_history(s)
