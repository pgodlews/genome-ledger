"""One-time migration: rename ego-centric sample ids to neutral given names.

Samples were once ided `father` / `mother` — relationships seen from one person, ego-centric
strings baked into CRAM/VCF filenames, per-tool output dirs, report names, the manifests and
the DuckDB `sample_id`. This renames them to perspective-neutral given names so the family
graph (pipeline/pedigree.py) can compute relationships relative to any person.

Dry-run by default — it prints every planned operation. `--apply` executes, after backing up
the manifests, pedigree and DuckDB. Idempotent: once renamed, the old id matches nothing and
the plan is empty. Renames key on whole filename *tokens* (so `father` in `mito_father.md` is
matched but `father` inside `grandfather` would not be), and a longer known id is one token:
renaming `Adam` leaves `Adam_vendor.cram` to the callset it belongs to.

A mapping that could merge or split people is refused before anything changes: chains and
swaps, two ids onto one, a target that already exists, a person renamed without their other
callsets, or a file rename that would land on an existing file.

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
from typing import Iterable

from . import config, pedigree
from .util import log, read_tsv

# Generic default mapping template. Override by editing here or extending the CLI.
DEFAULT_MAP = {"father": "SAMPLE001", "mother": "SAMPLE002"}
_VALID_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")


def _scan_roots() -> list[Path]:
    """Every root whose file/dir names can carry a sample id.

    The reports tree comes from `config.REPORT_DIRS`; the per-sample *data* dirs are
    enumerated explicitly. SHARE_DIR is left out on purpose: those are copies already handed
    to someone, named by `share --name`, not by the sample id."""
    c = config
    data_roots = [c.ALIGNED_DIR, c.CALLED_DIR, c.NORMALIZED_DIR, c.RAW_DIR,
                  c.ANNOTATED_DIR, c.PLINK2_DIR, c.PGX_DIR, c.ANCESTRY_DIR, c.TRAITS_DIR,
                  c.HLA_DIR, c.HLA_TYPE_DIR, c.PRS_DIR, c.REPEATS_DIR, c.SV_DIR,
                  c.CALLABLE_DIR, c.PHASE_DIR, c.FORCECALL_DIR, c.MITO_DIR, c.DENOVO_DIR,
                  c.SEGREGATION_DIR]
    return data_roots + [d for d in c.REPORT_DIRS if d not in data_roots]


def _token_re(ids: Iterable[str]) -> re.Pattern:
    # A whole-token match: an id not flanked by alphanumerics. Underscores and dots count as
    # delimiters, so `father` matches in `father.cram`, `mito_father.md`, `full_father.md`,
    # but not inside `grandfather`. Longest first, because `_` is also part of a callset id:
    # `Adam_vendor` has to match as itself rather than as `Adam` plus a suffix.
    alts = "|".join(re.escape(i) for i in sorted(ids, key=len, reverse=True))
    return re.compile(rf"(?<![A-Za-z0-9])(?:{alts})(?![A-Za-z0-9])")


def retoken(name: str, mapping: dict[str, str], keep: Iterable[str] = ()) -> str:
    """Replace every mapped id in `name` in one pass, so a replacement is never re-matched.
    `keep` are other known ids: they match as whole tokens and stay as they are."""
    return _token_re({*mapping, *keep}).sub(lambda m: mapping.get(m[0], m[0]), name)


def plan_renames(mapping: dict[str, str], roots: list[Path] | None = None,
                 keep: Iterable[str] = ()) -> list[tuple[Path, Path]]:
    """[(src, dst)] for every file/dir whose basename carries an old id token, deepest first
    so inner files rename before their parent directory."""
    keep = set(keep)
    ops: list[tuple[Path, Path]] = []
    for root in (roots if roots is not None else _scan_roots()):
        if not root.exists():
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            for name in filenames + dirnames:
                new = retoken(name, mapping, keep)
                if new != name:
                    ops.append((Path(dirpath) / name, Path(dirpath) / new))
    ops.sort(key=lambda op: len(op[0].parts), reverse=True)
    return ops


def _known_ids() -> set[str]:
    """Every id the manifests, the family graph or the database knows about."""
    ids: set[str] = set()
    for path, id_col, _paths, ref_cols in _MANIFEST_COLUMNS:
        for row in read_tsv(path):
            ids.update(row.get(c) or "" for c in (id_col, *ref_cols))
    if config.DUCKDB_FILE.exists():
        import duckdb
        con = duckdb.connect(str(config.DUCKDB_FILE), read_only=True)
        try:
            tables = {r[0] for r in con.execute("SHOW TABLES").fetchall()}
            for table in _SAMPLE_TABLES:
                if table in tables:
                    ids.update(r[0] for r in con.execute(
                        f"SELECT DISTINCT sample_id FROM {table}").fetchall() if r[0])
        finally:
            con.close()
    ids.discard("")
    return ids


def _refuse_unsafe_mapping(mapping: dict[str, str], known: set[str]) -> None:
    """Stop before anything changes if this mapping could merge or split people."""
    targets = list(mapping.values())
    present = {old for old in mapping if old in known}
    problems = []
    if bad := sorted(n for n in targets if not _VALID_ID.fullmatch(n)):
        problems.append(f"not a usable id: {', '.join(bad)}")
    if chained := sorted(set(mapping) & set(targets)):
        problems.append(f"both renamed and a rename target (chain or swap): {', '.join(chained)}")
    if shared := sorted({n for n in targets if targets.count(n) > 1}):
        problems.append(f"more than one id renamed to: {', '.join(shared)}")
    if taken := sorted({mapping[old] for old in present if mapping[old] in known}):
        problems.append(f"target id already exists: {', '.join(taken)}")
    # The same ownership rule the family graph uses, with the renamed ids counted as people.
    people = {*pedigree.load().people, *present}
    owners = pedigree.Family({p: pedigree.Person(p, p, "", None, None, None) for p in people},
                             sequenced=set())
    left = sorted(k for k in known - set(mapping)
                  if owners.person_of(k) != k and owners.person_of(k) in present)
    if left:
        problems.append("callsets of a renamed person are not in the mapping: "
                        + ", ".join(left))
    if problems:
        raise SystemExit("family-migrate refused, nothing was changed:\n  "
                         + "\n  ".join(problems))


def _refuse_collisions(renames: list[tuple[Path, Path]]) -> None:
    """A rename onto an existing file would strand one of the two, so stop before the first."""
    seen: set[Path] = set()
    clash = []
    for src, dst in renames:
        # samestat: on a case-insensitive disk `adam` "exists" when renaming `Adam` to it.
        # lstat throughout, so a dangling symlink counts as the file it is.
        if dst in seen or (os.path.lexists(dst)
                           and not os.path.samestat(os.lstat(src), os.lstat(dst))):
            clash.append(f"{src} -> {dst}")
        seen.add(dst)
    if clash:
        raise SystemExit(f"family-migrate refused, nothing was changed: {len(clash)} rename(s) "
                         "would land on an existing file:\n  " + "\n  ".join(clash[:10]))


# Per-manifest columns a rename must follow: the id column (exact-match replacement) plus any
# column holding a *path* into the renamed tree (token substitution, the same rule applied to
# the filenames themselves). Rewriting only the id column leaves the paths pointing at files
# plan_renames has just moved — normalize.py reads samples.tsv's raw_path directly, so a
# migrated sample could no longer be re-normalized.
# Columns of family.tsv that hold *another person's* id and so must be remapped too. They are
# id-valued, not path-valued, so they need exact-match replacement like the id column itself —
# not the token substitution _rewrite_row applies to path columns.
_FAMILY_ID_REFS = ("father", "mother", "partner")

class _ManifestColumns:
    """Read from config on each use: `set_genomes_root()` repoints these paths after import."""

    def __iter__(self):
        c = config
        return iter([
            (c.SAMPLES_TSV, "sample_id", ("source", "raw_path"), ()),
            (c.INCOMING_MANIFEST, "sample", (), ()),
            (c.RUNS_TSV, "sample_id", ("annotated_path",), ()),
            (c.CONSENT_FILE, "sample_id", (), ()),
            # family.tsv is the family GRAPH, and pedigree.ped is regenerated from it at the
            # end of migrate() — so leaving it out meant a migration could rewrite every other
            # manifest and then rebuild the pedigree keyed to the OLD ids, with no backup of
            # the original graph anywhere. Its id appears in four columns: the person's own id
            # and three references to other people, all of which the rename has to follow.
            (c.FAMILY_FILE, "id", (), _FAMILY_ID_REFS),
        ])


_MANIFEST_COLUMNS = _ManifestColumns()

def _column_indexes(header: list[str], id_col: str, path_cols: tuple[str, ...],
                    id_ref_cols: tuple[str, ...] = ()) -> tuple[int, tuple[int, ...],
                                                                tuple[int, ...]]:
    """(id column index, path column indexes, id-reference column indexes) — columns absent
    from an older manifest's header are simply skipped."""
    return (header.index(id_col),
            tuple(header.index(c) for c in path_cols if c in header),
            tuple(header.index(c) for c in id_ref_cols if c in header))


def _retoken_path(value: str, mapping: dict[str, str], keep: Iterable[str] = ()) -> str:
    """Follow the rename in a path cell, touching only what the migration actually moves.
    Above GENOMES_ROOT a matching folder (`/home/Adam/...`) is not ours to rename, and a file
    that exists outside the tree stays exactly where it is."""
    root = f"{config.GENOMES_ROOT}/"
    if value.startswith(root):
        return root + retoken(value[len(root):], mapping, keep)
    if Path(value).exists():
        return value
    return retoken(value, mapping, keep)


def _rewrite_row(fields: list[str], cols: tuple[int, tuple[int, ...], tuple[int, ...]],
                 mapping: dict[str, str], keep: Iterable[str] = ()) -> bool:
    """Rewrite one row's id, path and id-reference cells in place. True if anything changed."""
    ci, path_is, ref_is = cols
    changed = False
    if ci < len(fields) and fields[ci] in mapping:
        fields[ci] = mapping[fields[ci]]
        changed = True
    for pi in path_is:
        if pi < len(fields) and (new := _retoken_path(fields[pi], mapping, keep)) != fields[pi]:
            fields[pi] = new
            changed = True
    for ri in ref_is:   # exact match: these cells hold another person's id, not a path
        if ri < len(fields) and fields[ri] in mapping:
            fields[ri] = mapping[fields[ri]]
            changed = True
    return changed


def _tsv_edits(mapping: dict[str, str], keep: Iterable[str] = ()
               ) -> list[tuple[Path, str, tuple[str, ...], tuple[str, ...], int]]:
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
        n = sum(1 for ln in lines[1:] if _rewrite_row(ln.split("\t"), cols, mapping, keep))
        if n:
            edits.append((path, id_col, path_cols, ref_cols, n))
    return edits


def _apply_tsv(path: Path, id_col: str, path_cols: tuple[str, ...],
               ref_cols: tuple[str, ...], mapping: dict[str, str],
               keep: Iterable[str] = ()) -> None:
    lines = path.read_text().splitlines()
    cols = _column_indexes(lines[0].split("\t"), id_col, path_cols, ref_cols)
    out = [lines[0]]
    for ln in lines[1:]:
        f = ln.split("\t")
        _rewrite_row(f, cols, mapping, keep)
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
    from .load import transaction
    con = duckdb.connect(str(config.DUCKDB_FILE))
    try:
        tables = {r[0] for r in con.execute("SHOW TABLES").fetchall()}
        with transaction(con):   # all tables follow the rename, or none of them do
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
    known = _known_ids()
    _refuse_unsafe_mapping(mapping, known)
    keep = known | set(mapping.values())
    renames = plan_renames(mapping, keep=keep)
    _refuse_collisions(renames)
    tsvs = _tsv_edits(mapping, keep)
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
        src.rename(dst)
    for path, id_col, path_cols, ref_cols, _n in tsvs:
        _apply_tsv(path, id_col, path_cols, ref_cols, mapping, keep)
    _apply_duckdb(mapping)
    pedigree.generate_pedigree_ped()              # regenerate from family.tsv
    if reheader_vcfs:
        _reheader_vcfs(mapping)
    log.info("Migration applied. pedigree.ped regenerated from %s.", config.FAMILY_FILE)
