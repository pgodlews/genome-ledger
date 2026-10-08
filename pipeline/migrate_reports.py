"""One-time move of reports from the old layout into reports/md/ (see reportpaths).

Before: per-snapshot reports sat in `reports/<snapshot>/` beside their HTML and the
dashboard pages, and every per-genome stage kept its report in its own working directory
(`pgx/<sample>/pgx_<sample>.md`, …). This finds those markdown files and files each one
where `reportpaths` now puts it, adding the front matter new reports are written with.

Markdown is MOVED, never regenerated: several of these reports cost hours of compute over
a CRAM, and a stage decides whether to re-run by whether its report exists. Old HTML is
deleted — it is derived, and `run.py render` rebuilds all of it in the new tree.

Dry-run by default; `--apply` performs it.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from . import config, reportpaths
from .util import front_matter, log

# kind -> the stage directory its report used to live in
_STAGE_DIRS = {
    "pgx": "PGX_DIR", "hla": "HLA_DIR", "repeats": "REPEATS_DIR", "sv": "SV_DIR",
    "mito": "MITO_DIR", "callable": "CALLABLE_DIR", "phase": "PHASE_DIR",
    "denovo": "DENOVO_DIR", "segregation": "SEGREGATION_DIR", "forcecall": "FORCECALL_DIR",
    "hla_typing": "HLA_TYPE_DIR", "traits": "TRAITS_DIR", "prs": "PRS_DIR",
    "ancestry": "ANCESTRY_DIR",
}
def _new_roots() -> set[str]:
    """Directories under reports/ that belong to the new layout: the two markdown roots and
    one per registered format. Everything else at that level is an old dated directory."""
    from .render import FORMATS
    return {"md", "md-history", *FORMATS}


def _old_snapshot_dirs() -> list[Path]:
    root = config.REPORTS_DIR
    if not root.is_dir():
        return []
    return sorted(d for d in root.iterdir() if d.is_dir() and d.name not in _new_roots())


def plan() -> tuple[list[tuple[Path, Path]], list[Path]]:
    """(moves, deletions): markdown files to re-file, and derived files to remove."""
    moves: list[tuple[Path, Path]] = []
    deletions: list[Path] = []

    for d in _old_snapshot_dirs():
        snap = d.name
        for md in sorted(d.glob("*.md")):
            stem = md.stem
            kind = next((k for k in ("full", "panel", "acmg") if stem.startswith(k + "_")), None)
            if stem == "panel_family" or stem.startswith(("delta_", "validation_")):
                moves.append((md, reportpaths.family_report(stem, snap)))
            elif kind:
                moves.append((md, reportpaths.snapshot_report(kind, stem[len(kind) + 1:], snap)))
            else:
                log.warning("migrate-reports: leaving unrecognised file %s", md)
        deletions += sorted(d.glob("*.html"))
    if (config.REPORTS_DIR / "index.html").exists():
        deletions.append(config.REPORTS_DIR / "index.html")

    for kind, attr in _STAGE_DIRS.items():
        stage = getattr(config, attr)
        if not stage.is_dir():
            continue
        for md in sorted(stage.glob(f"*/{kind}_*.md")):
            sample = md.parent.name
            if md.stem != f"{kind}_{sample}":
                continue
            moves.append((md, reportpaths.genome_report(kind, sample)))
            if md.with_suffix(".html").exists():
                deletions.append(md.with_suffix(".html"))

    if config.INTERPRET_DIR.is_dir():
        for md in sorted(config.INTERPRET_DIR.glob("interpret_*_*.md")):
            model, _, snap = md.stem[len("interpret_"):].partition("_")
            if model and snap:
                moves.append((md, reportpaths.family_report(f"interpret_{model}", snap)))
                if md.with_suffix(".html").exists():
                    deletions.append(md.with_suffix(".html"))
    return moves, deletions


def migrate(apply: bool = False) -> int:
    moves, deletions = plan()
    if not moves and not deletions:
        log.info("migrate-reports: nothing in the old layout — already migrated.")
        return 0
    for src, dest in moves:
        log.info("  %s  ->  %s", src.relative_to(config.GENOMES_ROOT),
                 dest.relative_to(config.GENOMES_ROOT))
    log.info("migrate-reports: %d markdown report(s) to move, %d old HTML file(s) to "
             "delete.", len(moves), len(deletions))
    clashes = [dest for _src, dest in moves if dest.exists()]
    if clashes:
        raise SystemExit("migrate-reports: refusing to overwrite existing report(s): "
                         + ", ".join(str(c) for c in clashes[:5]))
    if not apply:
        log.info("Dry run. Re-run with --apply to perform it, then `run.py render`.")
        return 0

    for src, dest in moves:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(dest))
        text = dest.read_text()
        if not text.startswith("---\n"):
            first = text.splitlines()[0] if text else ""
            title = first[2:].strip() if first.startswith("# ") else ""
            dest.write_text("\n".join(front_matter(dest, title)) + "\n" + text)
    for f in deletions:
        f.unlink(missing_ok=True)
    for d in _old_snapshot_dirs():
        if not any(d.iterdir()):
            d.rmdir()
    snaps = reportpaths.md_snapshots()
    if snaps:
        reportpaths.archive_older(max(snaps))
    log.info("migrate-reports: done. Run `run.py render` to rebuild the HTML tree.")
    return 0
