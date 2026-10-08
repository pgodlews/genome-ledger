"""Stage 6: load VEP output into DuckDB (one table, partitioned by sample+snapshot).

The table is the query surface for reports, diffs, and pedigree analysis. DuckDB pages
it from disk, so a family's rows need not fit in memory.
"""

from __future__ import annotations

from contextlib import contextmanager

import duckdb

from . import config
from . import schema
from .annotate import annotated_path
from .util import log, read_tsv

TABLE = "variants"

SCHEMA = f"""
CREATE TABLE IF NOT EXISTS {TABLE} (
    snapshot_id           VARCHAR,
    sample_id             VARCHAR,
    chrom                 VARCHAR,
    pos                   BIGINT,
    ref                   VARCHAR,
    alt                   VARCHAR,
    gene                  VARCHAR,
    consequence           VARCHAR,
    clinvar_sig           VARCHAR,
    clinvar_revstat       VARCHAR,
    clinvar_alleleid      VARCHAR,
    clinvar_disease       VARCHAR,
    gnomad_af             DOUBLE,
    sift                  VARCHAR,
    polyphen              VARCHAR,
    zygosity              VARCHAR
);
"""


def connect() -> duckdb.DuckDBPyConnection:
    config.DUCKDB_DIR.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(config.DUCKDB_FILE))
    con.execute(SCHEMA)
    for ddl in schema.COMPANION_TABLES:  # derived, canonical-key companion tables
        con.execute(ddl)
    return con


@contextmanager
def transaction(con):
    """Make a DELETE-then-INSERT replacement atomic.

    DuckDB auto-commits per statement, so a failure between the two leaves the target emptied
    of the very rows the INSERT would have restored. `load` did this by hand; copy_forward and
    the ledger/history/triage writers did not, so an interrupted scan could wipe a sample's
    rows for a snapshot rather than leave them as they were. BaseException, not Exception:
    SystemExit and KeyboardInterrupt must roll back too.
    """
    con.execute("BEGIN TRANSACTION")
    try:
        yield con
    except BaseException:
        con.execute("ROLLBACK")
        raise
    con.execute("COMMIT")


# VEP --tab columns, in file order (controlled by --fields in annotate.py).
_VEP_COLUMNS = {
    "uploaded": "VARCHAR", "location": "VARCHAR", "allele": "VARCHAR",
    "symbol": "VARCHAR", "consequence": "VARCHAR",
    "gnomade": "VARCHAR", "gnomadg": "VARCHAR",
    "sift": "VARCHAR", "polyphen": "VARCHAR", "zyg": "VARCHAR",
    "clnsig": "VARCHAR", "clnrevstat": "VARCHAR", "clndn": "VARCHAR",
    "alleleid": "VARCHAR",
}


def _data_select(tsv) -> str:
    """The 14 variants data columns parsed from a VEP --tab TSV (no snapshot_id/sample_id —
    callers prepend those). Single source of truth for VEP-output parsing, shared by the full
    `load` and the incremental `load_subset`. Uploaded_variation is "chrom_pos_ref/alt" in VEP's
    own representation; gnomAD AF is the max of the exome/genome fields ('-' → NULL)."""
    return f"""
        SELECT chrom, pos, ref, alt, gene, consequence,
               clinvar_sig, clinvar_revstat, clinvar_alleleid, clinvar_disease,
               CASE WHEN afe IS NULL THEN afg
                    WHEN afg IS NULL THEN afe
                    ELSE greatest(afe, afg) END AS gnomad_af,
               sift, polyphen, zygosity
        FROM (
            SELECT
                split_part(uploaded, '_', 1) AS chrom,
                try_cast(split_part(uploaded, '_', 2) AS BIGINT) AS pos,
                split_part(split_part(uploaded, '_', 3), '/', 1) AS ref,
                split_part(split_part(uploaded, '_', 3), '/', 2) AS alt,
                nullif(symbol, '-') AS gene,
                nullif(consequence, '-') AS consequence,
                nullif(clnsig, '-') AS clinvar_sig,
                nullif(clnrevstat, '-') AS clinvar_revstat,
                nullif(alleleid, '-') AS clinvar_alleleid,
                nullif(clndn, '-') AS clinvar_disease,
                try_cast(split_part(nullif(gnomade, '-'), '&', 1) AS DOUBLE) AS afe,
                try_cast(split_part(nullif(gnomadg, '-'), '&', 1) AS DOUBLE) AS afg,
                nullif(sift, '-') AS sift,
                nullif(polyphen, '-') AS polyphen,
                nullif(zyg, '-') AS zygosity
            FROM {_read_csv(tsv)}
        )
        WHERE pos IS NOT NULL"""


def _read_csv(tsv) -> str:
    """The read_csv() call over a VEP --tab TSV, shared by the parse and the parse guard."""
    cols = ", ".join(f"'{k}': '{v}'" for k, v in _VEP_COLUMNS.items())
    tsv_sql = str(tsv).replace("'", "''")  # a quote in the path must not break the SQL
    return (f"read_csv('{tsv_sql}', delim='\\t', header=false, comment='#', "
            f"compression='gzip', auto_detect=false, columns={{{cols}}})")


def _guard_parsed(con, tsv) -> None:
    """Refuse to load a TSV whose Uploaded_variation column is not VEP's coordinate form.

    `_data_select` reads chrom/pos/ref/alt out of "chrom_pos_ref/alt", and its trailing
    `WHERE pos IS NOT NULL` used to discard anything else in silence. VEP only synthesizes
    that form when the input VCF has no ID: give it a record carrying "rs334" and the column
    holds "rs334" instead, so every identified variant vanished between annotation and the
    database with no error, no warning and a plausible-looking row count. `normalize` now
    strips ID so this cannot arise, but a TSV annotated before that change — or produced by
    hand — still can, and losing pathogenic variants quietly is the worst outcome available.
    """
    bad, total = con.execute(
        f"SELECT count(*) FILTER (WHERE try_cast(split_part(uploaded, '_', 2) AS BIGINT) "
        f"IS NULL), count(*) FROM {_read_csv(tsv)}").fetchone()
    if not bad:
        return
    example = con.execute(
        f"SELECT uploaded FROM {_read_csv(tsv)} "
        f"WHERE try_cast(split_part(uploaded, '_', 2) AS BIGINT) IS NULL LIMIT 1").fetchone()
    raise SystemExit(
        f"{tsv}: {bad:,} of {total:,} rows have an Uploaded_variation that is not "
        f"'chrom_pos_ref/alt' (e.g. {example[0]!r}) — VEP was given a VCF with a populated "
        f"ID column, so these rows carry an identifier instead of coordinates and cannot be "
        f"loaded. Re-normalize the sample (normalize now strips ID) and re-annotate.")


def load(sample: str, snapshot_id: str) -> int:
    """Load a VEP --tab output straight into DuckDB via read_csv — all parsing
    happens in SQL, which is dramatically faster than row-by-row inserts for the
    millions of variants in a WGS sample."""
    tsv = annotated_path(snapshot_id, sample)
    if not tsv.exists():
        raise SystemExit(f"No annotated output for {sample}/{snapshot_id}.")

    con = connect()
    try:
        _guard_parsed(con, tsv)          # outside the transaction: it only reads
        with transaction(con):
            con.execute(f"DELETE FROM {TABLE} WHERE sample_id = ? AND snapshot_id = ?",
                        [sample, snapshot_id])
            con.execute(f"INSERT INTO {TABLE} SELECT ?, ?, d.* FROM ({_data_select(tsv)}) d",
                        [snapshot_id, sample])
        # Outside the transaction: it is already committed here, so a failure in this
        # read-only count has nothing to roll back — and issuing ROLLBACK with no open
        # transaction raises a TransactionException that replaces the real error, turning a
        # successful load into a misleading traceback.
        n = con.execute(f"SELECT count(*) FROM {TABLE} WHERE sample_id = ? "
                        f"AND snapshot_id = ?", [sample, snapshot_id]).fetchone()[0]
    finally:
        con.close()
    log.info("Loaded %d variant rows for %s/%s.", n, sample, snapshot_id)
    return n


def load_subset(sample: str, snapshot_id: str, tsv) -> int:
    """Merge a re-annotated *subset* TSV into an existing (copied-forward) snapshot: replace
    just the rows whose keys the subset re-annotated, leaving the rest untouched. Used by the
    variant-level incremental scan (`incremental.reannotate_clinvar`)."""
    con = connect()
    try:
        _guard_parsed(con, tsv)          # outside the transaction: it only reads
        with transaction(con):
            con.execute(f"CREATE TEMP TABLE _stage AS SELECT * FROM ({_data_select(tsv)}) d")
            con.execute(f"DELETE FROM {TABLE} WHERE snapshot_id = ? AND sample_id = ? AND "
                        f"(chrom, pos, ref, alt) IN (SELECT chrom, pos, ref, alt FROM _stage)",
                        [snapshot_id, sample])
            con.execute(f"INSERT INTO {TABLE} SELECT ?, ?, * FROM _stage",
                        [snapshot_id, sample])
            n = con.execute("SELECT count(*) FROM _stage").fetchone()[0]
            con.execute("DROP TABLE _stage")
    finally:
        con.close()
    log.info("load_subset: merged %d re-annotated row(s) into %s/%s.", n, sample, snapshot_id)
    return n


def load_all(snapshot_id: str) -> None:
    from .util import find_row

    for row in read_tsv(config.SAMPLES_TSV):
        if annotated_path(snapshot_id, row["sample_id"]).exists():
            load(row["sample_id"], snapshot_id)
