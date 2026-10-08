"""Stage 4: create an immutable, dated database snapshot.

A snapshot pins the exact ClinVar VCF (downloaded today) plus the VEP cache /
gnomAD versions in use. Annotation runs reference a snapshot id, making them
reproducible and diffable.
"""

from __future__ import annotations

import fcntl
import json
import os
import re
import shutil
import urllib.request
from pathlib import Path

from . import config
from .util import log, require_tools, run, today


def snapshot_dir(snapshot_id: str) -> Path:
    return config.SNAPSHOTS_DIR / snapshot_id


def all_snapshots() -> list[str]:
    """Every pinned snapshot id, oldest first (ids sort lexicographically == by date)."""
    if not config.SNAPSHOTS_DIR.exists():
        return []
    return sorted(p.name for p in config.SNAPSHOTS_DIR.iterdir()
                  if p.is_dir() and (p / "manifest.json").exists())


# A snapshot dir is *pinned* the moment create_snapshot writes its manifest, but the scan
# that fills it — annotate, load, report, diff — runs afterwards and can fail. This marker
# records that the whole chain finished, so a half-built snapshot can be told apart from a
# finished one.
_SCAN_COMPLETE = ".scan_complete"


def mark_scan_complete(snapshot_id: str, samples: list[str] | None = None,
                       skipped: list[str] | None = None) -> None:
    """Record that `cmd_scan` finished every stage for this snapshot — and for whom.

    "The chain did not raise" is not "every sample was analysed": a sample with no
    normalized VCF is skipped with a warning, and the marker used to say the same thing
    either way. Naming the samples lets `releases.has_new_release` see that one of them has
    since become scannable and is still waiting."""
    (snapshot_dir(snapshot_id) / _SCAN_COMPLETE).write_text(json.dumps(
        {"completed": today(),
         "samples": None if samples is None else sorted(samples),
         "skipped": sorted(skipped or [])}, indent=2) + "\n")


def scan_completed(snapshot_id: str) -> bool:
    return (snapshot_dir(snapshot_id) / _SCAN_COMPLETE).exists()


def scanned_samples(snapshot_id: str) -> set[str] | None:
    """The samples a completed scan actually processed, or None when the marker does not
    say (no marker, or one written before it recorded them — a bare date)."""
    try:
        data = json.loads((snapshot_dir(snapshot_id) / _SCAN_COMPLETE).read_text())
    except (OSError, ValueError):
        return None
    if not isinstance(data, dict) or data.get("samples") is None:
        return None
    return set(data["samples"])


def latest_snapshot(require_complete: bool = False) -> str | None:
    """The newest pinned snapshot id.

    With `require_complete`, the newest one whose scan actually finished — what the release
    check needs, since a snapshot is pinned before the work that fills it. Snapshots created
    before the marker existed carry none, so a database with no markers at all reports its
    newest snapshot as incomplete and earns exactly one catch-up scan.
    """
    ids = all_snapshots()
    if require_complete:
        ids = [i for i in ids if scan_completed(i)]
    elif ids and not scan_completed(ids[-1]) and any(scan_completed(i) for i in ids[:-1]):
        # Only warn once markers are in use, so the pre-marker snapshots stay quiet.
        log.warning("Snapshot %s has no completion marker — its scan did not finish, so "
                    "reports and dashboards built from it may be partial.", ids[-1])
    return ids[-1] if ids else None


def previous_snapshot(curr: str) -> str | None:
    """The newest snapshot STRICTLY older than `curr`, or None if `curr` is the baseline.

    `latest_snapshot()` is the wrong predecessor for a scan: `create_snapshot` reuses an
    existing same-day id, so a second scan on one day would otherwise get prev == curr and
    hand `incremental.copy_forward` a snapshot to copy from itself (delete-then-copy-from-
    itself → every sample emptied). Always derive the predecessor from `curr`, never from
    "whatever is newest".

    Among the older ones, the newest whose scan COMPLETED. A scan that died partway leaves a
    pinned snapshot holding only the samples it reached; as a baseline it makes everyone it
    did not reach look newly added, and the diff leaves new samples out by design — so the
    retry after a failed scan reported no change for exactly the people whose
    reclassification was still undelivered. Falls back to the newest older snapshot only
    when none carries a marker (a database from before markers existed).
    """
    older = [s for s in all_snapshots() if s < curr]
    done = [s for s in older if scan_completed(s)]
    if done:
        if done[-1] != older[-1]:
            log.warning("Skipping unfinished snapshot(s) %s as a baseline for %s — comparing "
                        "against %s, the last scan that completed.",
                        ", ".join(older[older.index(done[-1]) + 1:]), curr, done[-1])
        return done[-1]
    return older[-1] if older else None


def create_snapshot(snapshot_id: str | None = None) -> str:
    require_tools("tabix", "bgzip")
    snapshot_id = snapshot_id or today()
    # Snapshot ordering is lexicographic (latest_snapshot, diff predecessor selection,
    # validate's truth<=replay guard all rely on name order == time order), so custom
    # ids must stay date-prefixed.
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}([._-][A-Za-z0-9._-]*)?", snapshot_id):
        raise SystemExit(
            f"Snapshot id '{snapshot_id}' must be YYYY-MM-DD optionally followed by "
            "-suffix — downstream diff/validate logic orders snapshots lexicographically.")
    d = snapshot_dir(snapshot_id)
    manifest = d / "manifest.json"
    if manifest.exists():
        log.info("Snapshot %s already exists — reusing.", snapshot_id)
        return snapshot_id

    d.mkdir(parents=True, exist_ok=True)
    # Serialise creators of the SAME snapshot id. The per-process .part download already
    # stopped two runs interleaving bytes, but they still raced on everything after it: both
    # renamed onto the final clinvar.vcf.gz, both ran tabix over it, and both wrote the
    # manifest. tabix's .tbi write is not atomic, so one run could leave an index built
    # against the other's file. A scheduled scan overlapping a manual one is the realistic
    # way in. The manifest is re-checked under the lock, so the loser reuses rather than
    # redownloading ~200 MB.
    lock_path = d / ".create.lock"
    with open(lock_path, "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if manifest.exists():
            log.info("Snapshot %s was created by a concurrent run — reusing.", snapshot_id)
            return snapshot_id
        return _create_snapshot_locked(snapshot_id, d, manifest)


def _create_snapshot_locked(snapshot_id: str, d, manifest):
    """The body of create_snapshot, run while holding the per-snapshot creation lock."""
    clinvar = d / "clinvar.vcf.gz"
    log.info("Downloading ClinVar %s → %s", config.BUILD, clinvar)
    # Atomic download: to a .part file + rename, with a timeout — a killed/stalled run
    # must not leave a truncated VCF pinned under the final name forever.
    # Per-process: two create_snapshot runs racing on the same id would otherwise interleave
    # their bytes in one shared .part file and rename the mixture under the immutable name.
    tmp = d / f"{clinvar.name}.{os.getpid()}.part"
    try:
        with urllib.request.urlopen(config.CLINVAR_VCF_URL, timeout=120) as r, \
                open(tmp, "wb") as fh:
            shutil.copyfileobj(r, fh)
        os.replace(tmp, clinvar)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise
    # ClinVar ships bgzipped; (re)build the tabix index locally.
    run(["tabix", "-f", "-p", "vcf", str(clinvar)])

    meta = {
        "snapshot_id": snapshot_id,
        "created": today(),
        "build": config.BUILD,
        "clinvar_url": config.CLINVAR_VCF_URL,
        "clinvar_date": _clinvar_filedate(clinvar),
        "vep_cache_version": config.VEP_CACHE_VERSION,
        "gnomad_version": config.GNOMAD_VERSION,
    }
    manifest.write_text(json.dumps(meta, indent=2))
    log.info("Snapshot %s pinned: ClinVar=%s, VEP cache=%s",
             snapshot_id, meta["clinvar_date"], config.VEP_CACHE_VERSION)
    return snapshot_id


def _clinvar_filedate(clinvar: Path) -> str:
    """Pull the source release date from the ClinVar VCF header (##fileDate)."""
    import gzip

    with gzip.open(clinvar, "rt") as fh:
        for line in fh:
            if line.startswith("##fileDate"):
                return line.strip().split("=", 1)[1]
            if not line.startswith("#"):
                break
    return "unknown"


def load_manifest(snapshot_id: str) -> dict:
    return json.loads((snapshot_dir(snapshot_id) / "manifest.json").read_text())
