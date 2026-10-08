"""One-time migration: rename ego-centric sample ids to neutral given names.

Samples were once ided `father` / `mother` — relationships seen from one person, ego-centric
strings baked into CRAM/VCF filenames, per-tool output dirs, report names, the manifests and
the DuckDB `sample_id`. This renames them to perspective-neutral given names so the family
graph (pipeline/pedigree.py) can compute relationships relative to any person.

Dry-run by default — it prints every planned operation. `--apply` executes, after backing up
the manifests, pedigree and DuckDB. Idempotent: once renamed, the old id matches nothing and
the plan is empty. Renames key on whole filename *tokens* (so `father` in `mito_father.md` is
matched but `father` inside `grandfather` would not be).

Internal sample names (the CRAM `@RG SM:` tag and the VCF sample column) are left as the old
id: the pipeline keys on filenames / explicit `sample` args / the DuckDB column, not on the
embedded names, so this is cosmetic. `--reheader-vcfs` opts into rewriting the (cheap) VCF
sample columns; a CRAM (tens of GB) is never rewritten.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

from . import config, pedigree
from .util import log

# Generic default mapping template. Override by editing here or extending the CLI.
DEFAULT_MAP = {"father": "SAMPLE001", "mother": "SAMPLE002"}


def _scan_roots() -> list[Path]:
    """Every root whose file/dir names can carry a sample id.

    The reports tree comes from `config.REPORT_DIRS`; the per-sample *data* dirs are
    enumerated explicitly."""
    c = config
    data_roots = [c.ALIGNED_DIR, c.CALLED_DIR, c.NORMALIZED_DIR, c.RAW_DIR,
                  c.ANNOTATED_DIR, c.PLINK2_DIR]
    return data_roots + [d for d in c.REPORT_DIRS if d not in data_roots]


def _token_re(old: str) -> re.Pattern:
    # A whole-token match: `old` not flanked by alphanumerics. Underscores and dots count as
    # delimiters, so `father` matches in `father.cram`, `mito_father.md`, `full_father.md`,
    # but not inside `grandfather`.
    return re.compile(rf"(?<![A-Za-z0-9]){re.escape(old)}(?![A-Za-z0-9])")


def retoken(name: str, mapping: dict[str, str]) -> str:
    for old, new in mapping.items():
        name = _token_re(old).sub(new, name)
    return name


def plan_renames(mapping: dict[str, str], roots: list[Path] | None = None) -> list[tuple[Path, Path]]:
    """[(src, dst)] for every file/dir whose basename carries an old id token, deepest first
    so inner files rename before their parent directory."""
    pats = [_token_re(old) for old in mapping]
    hits: list[Path] = []
    for root in (roots if roots is not None else _scan_roots()):
        if not root.exists():
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            for name in filenames + dirnames:
                if any(p.search(name) for p in pats):
                    hits.append(Path(dirpath) / name)
    hits.sort(key=lambda p: len(p.parts), reverse=True)
    ops = []
    for src in hits:
        dst = src.parent / retoken(src.name, mapping)
        if dst != src:
            ops.append((src, dst))
    return ops


# Per-manifest columns a rename must follow: the id column (exact-match replacement) plus any
# column holding a *path* into the renamed tree (token substitution, the same rule applied to
# the filenames themselves). Rewriting only the id column leaves the paths pointing at files
# plan_renames has just moved — normalize.py reads samples.tsv's raw_path directly, so a
# migrated sample could no longer be re-normalized.
# Columns of family.tsv that hold *another person's* id and so must be remapped too. They are
# id-valued, not path-valued, so they need exact-match replacement like the id column itself —
# not the token substitution _rewrite_row applies to path columns.
_FAMILY_ID_REFS = ("father", "mother", "partner")

_MANIFEST_COLUMNS: list[tuple[Path, str, tuple[str, ...], tuple[str, ...]]] = [
    (config.SAMPLES_TSV, "sample_id", ("source", "raw_path"), ()),
    (config.INCOMING_MANIFEST, "sample", (), ()),
    (config.RUNS_TSV, "sample_id", ("annotated_path",), ()),
    (config.CONSENT_FILE, "sample_id", (), ()),
    # family.tsv is the family GRAPH, and pedigree.ped is regenerated from it at the end of
    # migrate() — so leaving it out meant a migration could rewrite every other manifest and
    # then rebuild the pedigree keyed to the OLD ids, with no backup of the original graph
    # anywhere. Its id appears in four columns: the person's own id and three references to
    # other people, all of which the rename has to follow.
    (config.FAMILY_FILE, "id", (), _FAMILY_ID_REFS),
]

def _column_indexes(header: list[str], id_col: str, path_cols: tuple[str, ...],
                    id_ref_cols: tuple[str, ...] = ()) -> tuple[int, tuple[int, ...],
                                                                tuple[int, ...]]:
    """(id column index, path column indexes, id-reference column indexes) — columns absent
    from an older manifest's header are simply skipped."""
    return (header.index(id_col),
            tuple(header.index(c) for c in path_cols if c in header),
            tuple(header.index(c) for c in id_ref_cols if c in header))


def _rewrite_row(fields: list[str], cols: tuple[int, tuple[int, ...], tuple[int, ...]],
                 mapping: dict[str, str]) -> bool:
    """Rewrite one row's id, path and id-reference cells in place. True if anything changed."""
    ci, path_is, ref_is = cols
    changed = False
    if ci < len(fields) and fields[ci] in mapping:
        fields[ci] = mapping[fields[ci]]
        changed = True
    for pi in path_is:
        if pi < len(fields) and (new := retoken(fields[pi], mapping)) != fields[pi]:
            fields[pi] = new
            changed = True
    for ri in ref_is:   # exact match: these cells hold another person's id, not a path
        if ri < len(fields) and fields[ri] in mapping:
            fields[ri] = mapping[fields[ri]]
            changed = True
    return changed


def _tsv_edits(mapping: dict[str, str]) -> list[tuple[Path, str, tuple[str, ...],
                                                      tuple[str, ...], int]]:
    """[(path, id_col, path_cols, id_ref_cols, n_rows_affected)] for manifests needing a
    rewrite."""
    edits = []
    for path, id_col, path_cols, ref_cols in _MANIFEST_COLUMNS:
        if not path.exists():
            continue
        lines = path.read_text().splitlines()
        if not lines:
            continue
        header = lines[0].split("\t")
        if id_col not in header:
            continue
        cols = _column_indexes(header, id_col, path_cols, ref_cols)
        n = sum(1 for ln in lines[1:] if _rewrite_row(ln.split("\t"), cols, mapping))
        if n:
            edits.append((path, id_col, path_cols, ref_cols, n))
    return edits


def _apply_tsv(path: Path, id_col: str, path_cols: tuple[str, ...],
               ref_cols: tuple[str, ...], mapping: dict[str, str]) -> None:
    lines = path.read_text().splitlines()
    cols = _column_indexes(lines[0].split("\t"), id_col, path_cols, ref_cols)
    out = [lines[0]]
    for ln in lines[1:]:
        f = ln.split("\t")
        _rewrite_row(f, cols, mapping)
        out.append("\t".join(f))
    path.write_text("\n".join(out) + "\n")


# DuckDB tables carrying a sample_id that must follow the rename — otherwise the
# incremental-annotation ledger (scan_inputs) and the DELTA audit trail (triage_results)
# silently stay keyed to the old id.
_SAMPLE_TABLES = ("variants", "triage_results", "annotation_provenance", "scan_inputs")


def _duckdb_counts(mapping: dict[str, str]) -> list[tuple[str, str, int]]:
    """[(table, old_id, n_rows)] across every table with a sample_id column."""
    if not config.DUCKDB_FILE.exists():
        return []
    import duckdb
    con = duckdb.connect(str(config.DUCKDB_FILE), read_only=True)
    try:
        tables = {r[0] for r in con.execute("SHOW TABLES").fetchall()}
        out = []
        for table in _SAMPLE_TABLES:
            if table not in tables:
                continue
            for old in mapping:
                n = con.execute(
                    f"SELECT count(*) FROM {table} WHERE sample_id = ?", [old]
                ).fetchone()[0]
                if n:
                    out.append((table, old, n))
        return out
    finally:
        con.close()


def _apply_duckdb(mapping: dict[str, str]) -> None:
    if not config.DUCKDB_FILE.exists():
        return
    import duckdb
    con = duckdb.connect(str(config.DUCKDB_FILE))
    try:
        tables = {r[0] for r in con.execute("SHOW TABLES").fetchall()}
        for table in _SAMPLE_TABLES:
            if table not in tables:
                continue
            for old, new in mapping.items():
                con.execute(f"UPDATE {table} SET sample_id = ? WHERE sample_id = ?",
                            [new, old])
    finally:
        con.close()


def _reheader_vcfs(mapping: dict[str, str]) -> None:
    """Rewrite the VCF sample column for the renamed VCFs (cheap — header only)."""
    for new in mapping.values():
        for vcf in [config.CALLED_DIR / f"{new}.{config.BUILD}.vcf.gz",
                    config.NORMALIZED_DIR / f"{new}.norm.vcf.gz",
                    config.FORCECALL_DIR / f"{new}.{config.BUILD}.vcf.gz"]:
            if not vcf.exists():
                continue
            samples = vcf.parent / f".{new}.samples.txt"
            samples.write_text(new + "\n")
            tmp = vcf.with_suffix(vcf.suffix + ".reh.tmp")
            try:
                subprocess.run(["bcftools", "reheader", "-s", str(samples), "-o", str(tmp),
                                str(vcf)], check=True, capture_output=True)
                tmp.rename(vcf)
                subprocess.run(["bcftools", "index", "-ft", str(vcf)], capture_output=True)
            except Exception as e:  # noqa: BLE001 — cosmetic, never fatal
                log.warning("VCF reheader failed for %s (%s)", vcf, e)
                tmp.unlink(missing_ok=True)
            finally:
                samples.unlink(missing_ok=True)


def _backup() -> None:
    # Checkpoint first: a plain copy of the .duckdb file misses any un-checkpointed WAL —
    # the safety copy made right before an irreversible rename must be restorable.
    if config.DUCKDB_FILE.exists():
        import duckdb
        try:
            con = duckdb.connect(str(config.DUCKDB_FILE))
            con.execute("CHECKPOINT")
            con.close()
        except Exception as e:  # noqa: BLE001 — checkpoint is best-effort; copy anyway
            log.warning("DuckDB checkpoint before backup failed (%s)", e)
    # FAMILY_FILE included: pedigree.ped is a *generated* artifact rebuilt from it, so a
    # backup of the .ped alone does not preserve the graph it was generated from.
    for p in [config.SAMPLES_TSV, config.INCOMING_MANIFEST, config.PEDIGREE_FILE,
              config.RUNS_TSV, config.CONSENT_FILE, config.FAMILY_FILE, config.DUCKDB_FILE]:
        if p.exists():
            bak = p.with_suffix(p.suffix + ".pre-migrate.bak")
            if bak.exists():
                # A later --apply must not overwrite the *first* migration's safety copy —
                # that is the only pre-rename state that exists anywhere.
                i = 2
                while (alt := p.with_suffix(p.suffix + f".pre-migrate.{i}.bak")).exists():
                    i += 1
                bak = alt
            shutil.copy2(p, bak)
            log.info("backup: %s", bak)


def migrate(apply: bool = False, reheader_vcfs: bool = False,
            mapping: dict[str, str] | None = None) -> None:
    mapping = mapping or DEFAULT_MAP
    renames = plan_renames(mapping)
    tsvs = _tsv_edits(mapping)
    db = _duckdb_counts(mapping)

    log.info("=== family-migrate %s ===  %s",
             "(APPLY)" if apply else "(dry-run)",
             ", ".join(f"{o}->{n}" for o, n in mapping.items()))
    log.info("File/dir renames: %d", len(renames))
    for src, dst in renames:
        log.info("  mv  %s  ->  %s", src.name, dst.name)
    for path, id_col, path_cols, ref_cols, n in tsvs:
        log.info("Edit %s [%s]: %d row(s)", path,
                 ", ".join((id_col, *path_cols, *ref_cols)), n)
    for table, old, n in db:
        log.info("DuckDB: UPDATE %s sample_id '%s' -> '%s' (%d rows)",
                 table, old, mapping[old], n)
    if not renames and not tsvs and not db:
        log.info("Nothing to migrate — already neutral (idempotent).")
        return
    if not apply:
        log.info("Dry-run only. Re-run with --apply to execute.")
        return

    _backup()
    for src, dst in renames:
        if dst.exists():
            log.warning("skip (target exists): %s -> %s", src, dst)
            continue
        src.rename(dst)
    for path, id_col, path_cols, ref_cols, _n in tsvs:
        _apply_tsv(path, id_col, path_cols, ref_cols, mapping)
    _apply_duckdb(mapping)
    pedigree.generate_pedigree_ped()              # regenerate from family.tsv
    if reheader_vcfs:
        _reheader_vcfs(mapping)
    log.info("Migration applied. pedigree.ped regenerated from %s.", config.FAMILY_FILE)
