"""Stage 10: detect when upstream databases publish a new release.

ClinVar is the fast mover (weekly); gnomAD and the VEP cache change rarely and
are bumped by editing config.py. We compare the newest upstream ClinVar weekly
file against the ClinVar date pinned in our latest snapshot.
"""

from __future__ import annotations

import re
import urllib.request

from . import config
from .normalize import normalized_path
from .snapshot import latest_snapshot, load_manifest, scan_completed, scanned_samples
from .util import log, read_tsv

# Filenames look like: clinvar_20260815.vcf.gz
_WEEKLY_RE = re.compile(r"clinvar_(\d{8})\.vcf\.gz")


def newest_clinvar_release() -> str | None:
    try:
        with urllib.request.urlopen(config.CLINVAR_WEEKLY_INDEX, timeout=20) as resp:
            html = resp.read().decode("utf-8", "replace")
    except Exception as e:  # noqa: BLE001
        log.warning("Could not reach ClinVar weekly index: %s", e)
        return None
    dates = sorted(_WEEKLY_RE.findall(html))
    return dates[-1] if dates else None


def unscanned_samples(snap: str) -> list[str]:
    """Ingested samples that can be scanned now but that completed scan `snap` left out.

    A scan skips a sample with no normalized VCF and still completes for everyone else.
    Without this, normalizing that sample afterwards changed nothing: the release was
    already "processed", so it sat unanalysed until ClinVar next published. Empty when the
    marker predates the sample list — there is nothing to compare against."""
    done = scanned_samples(snap)
    if done is None:
        return []
    return sorted(r["sample_id"] for r in read_tsv(config.SAMPLES_TSV)
                  if r.get("sample_id") and r["sample_id"] not in done
                  and normalized_path(r["sample_id"]).exists())


def has_new_release() -> bool:
    """True if upstream ClinVar is newer than our latest *completed* snapshot.

    Completed, not merely pinned. `create_snapshot` writes the manifest before annotate,
    load, report and diff run, so a scan that died partway still left a snapshot whose
    clinvar_date matched upstream — and this check then reported "up to date" and stood the
    weekly job down. The failed scan suppressed its own retry, while the half-filled
    snapshot stayed newest for reports and dashboards.
    """
    newest_pinned = latest_snapshot()
    if newest_pinned is not None and not scan_completed(newest_pinned):
        log.warning("Snapshot %s was pinned but its scan never completed — re-running "
                    "rather than treating that ClinVar release as processed.", newest_pinned)
        return True
    snap = latest_snapshot(require_complete=True)
    if snap is None:
        log.info("No completed snapshot yet — a first scan is due.")
        return True
    waiting = unscanned_samples(snap)
    if waiting:
        log.warning("%d sample(s) are normalized but were not part of scan %s (%s) — "
                    "running a scan for them rather than waiting for the next release.",
                    len(waiting), snap, ", ".join(waiting))
        return True
    newest = newest_clinvar_release()
    if newest is None:
        log.error("ClinVar release check failed (network?) — skipping this scan; if this "
                  "persists, new releases will be missed silently.")
        return False
    current = load_manifest(snap).get("clinvar_date", "")
    current_digits = re.sub(r"\D", "", current)  # e.g. "2026-08-01" -> "20260801"
    is_new = newest > current_digits
    log.info("ClinVar upstream=%s, pinned=%s → %s",
             newest, current_digits or "?", "NEW" if is_new else "up to date")
    return is_new
