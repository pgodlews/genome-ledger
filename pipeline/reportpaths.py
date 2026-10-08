"""Where reports live. One tree per format, the same relative path in each.

    reports/
      md/                         the source of truth, and the folder to hand to a retriever
        <snapshot>/                 the LATEST snapshot only
          _family/                    delta_…md, panel_family.md, inheritance.md, …
          <person>/                   full_<callset>.md, panel_<callset>.md, acmg_<callset>.md
        genome/                     reports that depend on the genome, not on a ClinVar release
          <person>/                   pgx_<callset>.md, prs_<callset>.md, hla_<callset>.md, …
      md-history/                 older snapshots' markdown, same shape as md/<snapshot>/
        <snapshot>/…
      html/                       derived from the markdown, plus the navigation pages
        index.html
        <snapshot>/                 index.html, index_<person>.html, inheritance.html, <person>/…
        genome/<person>/…

`md/` holds only the latest snapshot on purpose. A snapshot is one ClinVar release, so an
older snapshot's report states classifications that have since been superseded; a retriever
pointed at a folder containing every snapshot will surface "uncertain significance" for a
variant that is pathogenic today. History is kept, one directory over, where nothing reads
it by accident. HTML keeps every snapshot, because there a person chooses which one to read.

Stage directories (`pgx/`, `prs/`, …) hold working data only — VCFs, tool output, caches —
and no reports. Every writer asks this module for its path; nothing builds one by hand.
See docs/report-formats.md for the file contract and how to add a format.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from . import config
from .util import log

GENOME = "genome"     # directory for per-genome (undated) reports
FAMILY = "_family"    # directory, inside a snapshot, for reports about the whole family


def md_root() -> Path:
    return config.REPORTS_DIR / "md"


def history_root() -> Path:
    return config.REPORTS_DIR / "md-history"


def html_root() -> Path:
    return config.REPORTS_DIR / "html"


def format_root(fmt: str) -> Path:
    """Root of a derived format's tree: reports/<fmt>/."""
    return config.REPORTS_DIR / fmt


def _person(sample: str) -> str:
    """The person a callset id belongs to — reports are filed per person, named per callset.
    Tolerant like consent._person: a missing or unreadable family graph must not stop a
    report being written, and filing under the id itself is the natural fallback."""
    try:
        from . import pedigree
        return pedigree.load().person_of(sample)
    except (Exception, SystemExit):  # noqa: BLE001
        return sample


def md_snapshots() -> list[str]:
    """Snapshot ids that have a markdown directory, in either root."""
    out = set()
    for root in (md_root(), history_root()):
        if root.is_dir():
            out |= {d.name for d in root.iterdir() if d.is_dir() and d.name != GENOME}
    return sorted(out)


def snapshot_md_dir(snap: str) -> Path:
    """The markdown directory for a snapshot: under md/ if it is the latest, md-history/ if
    not. An existing directory wins, so a report is always written beside its siblings;
    `archive_older` is what moves a snapshot across when a newer one arrives."""
    cur, hist = md_root() / snap, history_root() / snap
    if cur.exists():
        return cur
    if hist.exists():
        return hist
    from .snapshot import all_snapshots
    try:
        known = set(all_snapshots())
    except (Exception, SystemExit):  # noqa: BLE001 — no snapshots dir yet
        known = set()
    newest = max(known | set(md_snapshots()) | {snap})
    return cur if snap == newest else hist


def snapshot_report(kind: str, sample: str, snap: str) -> Path:
    """One callset's report for one snapshot: FULL, panel, ACMG."""
    return snapshot_md_dir(snap) / _person(sample) / f"{kind}_{sample}.md"


def family_report(name: str, snap: str) -> Path:
    """A snapshot-level report about the whole family: DELTA, family panel, inheritance."""
    return snapshot_md_dir(snap) / FAMILY / f"{name}.md"


def genome_report(kind: str, sample: str) -> Path:
    """One callset's report for a stage that depends only on the genome."""
    return md_root() / GENOME / _person(sample) / f"{kind}_{sample}.md"


def relative(md: Path) -> Path | None:
    """A report's path relative to its markdown root (md/ or md-history/), or None if the
    file is not inside the reports tree at all."""
    md = Path(md)
    for root in (md_root(), history_root()):
        try:
            return md.resolve().relative_to(root.resolve())
        except ValueError:
            continue
    return None


def derived(md: Path, fmt: str, ext: str) -> Path:
    """Where format `fmt` puts its rendering of `md`: reports/<fmt>/<same path>.<ext>.
    A markdown file outside the reports tree gets a sibling, which is what tests and ad-hoc
    renders want."""
    rel = relative(md)
    if rel is None:
        return Path(md).with_suffix(f".{ext}")
    return (format_root(fmt) / rel).with_suffix(f".{ext}")


def html_for(md: Path) -> Path:
    return derived(md, "html", "html")


def meta(md: Path) -> dict | None:
    """What a report is, read off where it is filed: kind, callset, person, scope, snapshot.
    None for a file outside the reports tree."""
    rel = relative(md)
    if rel is None or len(rel.parts) < 3:
        return None
    top, owner, stem = rel.parts[0], rel.parts[1], rel.stem
    if owner == FAMILY:
        return {"report": stem.split("_")[0] if stem.startswith(("delta_", "validation_",
                                                                  "interpret_")) else stem,
                "scope": "family", "snapshot": top}
    kind = stem[: -len(owner) - 1] if stem.endswith("_" + owner) else None
    callset = owner
    if kind is None:
        # <kind>_<callset> where the callset is an alternate one of this person
        # (full_Adam_vendor under Adam/): the callset starts at the person's id.
        i = stem.find("_" + owner)
        kind, callset = (stem[:i], stem[i + 1:]) if i > 0 else (stem, owner)
    out = {"report": kind, "person": owner, "callset": callset}
    if top == GENOME:
        out["scope"] = "genome"
    else:
        out["scope"], out["snapshot"] = "snapshot", top
    return out


def all_md() -> list[Path]:
    """Every report markdown file, both roots."""
    out: list[Path] = []
    for root in (md_root(), history_root()):
        if root.is_dir():
            out += sorted(root.rglob("*.md"))
    return out


def archive_older(latest: str | None) -> list[str]:
    """Move every snapshot directory in md/ other than `latest` into md-history/, so md/
    never shows a retriever a superseded classification. Returns the snapshots moved."""
    root = md_root()
    if not latest or not root.is_dir():
        return []
    moved = []
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        if d.name in (GENOME, latest) or d.name > latest:
            continue
        dest = history_root() / d.name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            # A partial earlier move: merge file by file rather than nest or clobber.
            for f in sorted(p for p in d.rglob("*") if p.is_file()):
                target = dest / f.relative_to(d)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(f), str(target))
            shutil.rmtree(d)
        else:
            shutil.move(str(d), str(dest))
        moved.append(d.name)
    if moved:
        log.info("Archived superseded snapshot markdown to %s: %s",
                 history_root().name, ", ".join(moved))
    return moved
