"""Enformer lookup: deep-learning non-coding / regulatory variant effect (DeepMind), GRCh38.

The one in-silico signal in this pipeline that is *computed on demand on the GPU box* rather
than joined from a public precomputed table — because non-coding effect isn't precomputed
genome-wide anywhere usable. Enformer predicts ~5,313 regulatory tracks (CAGE expression,
DNase/ATAC accessibility, histone & TF ChIP) from 196,608 bp of sequence; the variant effect
is the change in those tracks between the reference and alternate allele. This fills VEP +
ClinVar's biggest blind spot: promoter/enhancer/UTR variants a missense or splice predictor
structurally can't speak to. See pipeline/interpret.py (the stage that produces the table)
and docs/hardware-and-scope.md.

The `interpret` stage writes a small per-snapshot table (one row per scored variant) that this
module reads — the same lookup()/label() contract as alphamissense.py/spliceai.py, so reports
enrich themselves when present and degrade gracefully when absent.

Table (TSV, config.ENFORMER_SCORES_FILE): chrom  pos  ref  alt  delta_max  l2_center  top_track  top_desc
"""

from __future__ import annotations

import bisect
import csv
import json

from . import config

# Loaded once, lazily: {(chrom, pos, ref, alt): (delta_max, l2_center, top_track, top_desc)}.
_table: dict | None = None
# Calibration null: {"snv": [sorted Δmax], "indel": [sorted Δmax], ...} or {} if uncalibrated.
_calibration: dict | None = None


def is_snv(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT"


def available() -> bool:
    return config.ENFORMER_SCORES_FILE.exists()


def _load() -> dict:
    global _table
    if _table is not None:
        return _table
    _table = {}
    if not available():
        return _table
    with open(config.ENFORMER_SCORES_FILE) as fh:
        header = fh.readline()  # skip
        # csv.reader — the cache is written with csv.writer (interpret._merge_scores);
        # split("\t") would misparse a quoted field.
        for f in csv.reader(fh, delimiter="\t"):
            if len(f) < 7:
                continue
            try:
                key = (f[0], int(f[1]), f[2], f[3])
                _table[key] = (float(f[4]), float(f[5]), int(f[6]),
                               f[7] if len(f) > 7 else "")
            except ValueError:
                continue
    return _table


def lookup(chrom: str, pos, ref: str, alt: str):
    """(delta_max, l2_center, top_track, top_desc) for a scored variant, or None.
    delta_max is the largest |alt-ref| change across the central regulatory bins/tracks."""
    return _load().get((str(chrom), int(pos), ref, alt))


def calibration_available() -> bool:
    return bool(_load_calibration())


def _load_calibration() -> dict:
    global _calibration
    if _calibration is not None:
        return _calibration
    _calibration = {}
    if config.ENFORMER_CALIBRATION_FILE.exists():
        try:
            data = json.loads(config.ENFORMER_CALIBRATION_FILE.read_text())
            # Keep just the sorted null arrays we need for the percentile lookup.
            _calibration = {k: data[k] for k in ("snv", "indel") if data.get(k)}
        except (ValueError, OSError):
            _calibration = {}
    return _calibration


def percentile(delta_max: float, ref: str, alt: str):
    """Δmax as a percentile (0–100) against the benign null for this variant's stratum, or
    None if uncalibrated / the stratum is empty. The fraction of common (≈benign) variants of
    the same kind whose effect is ≤ this one — higher means a more unusual perturbation."""
    null = _load_calibration().get("snv" if is_snv(ref, alt) else "indel")
    if not null:
        return None
    return 100.0 * bisect.bisect_right(null, delta_max) / len(null)


def band(delta_max: float, ref: str, alt: str) -> str:
    """'high' / 'moderate' / 'low'. By percentile vs the benign null when calibrated;
    otherwise by the raw Δ cut-offs (the pre-calibration fallback)."""
    pct = percentile(delta_max, ref, alt)
    if pct is not None:
        if pct >= config.ENFORMER_PCTL_HIGH:
            return "high"
        return "moderate" if pct >= config.ENFORMER_PCTL_MODERATE else "low"
    if delta_max >= config.ENFORMER_DELTA_HIGH:
        return "high"
    return "moderate" if delta_max >= config.ENFORMER_DELTA_MODERATE else "low"


def is_regulatory_hit(en, ref: str, alt: str) -> bool:
    """A predicted regulatory perturbation at/above the moderate band."""
    return bool(en) and band(en[0], ref, alt) in ("moderate", "high")


def label(en, ref: str = "", alt: str = "") -> str:
    """Human-readable cell, e.g. 'Δ6.98 p99.7 high (CAGE)' or '—'. Includes the percentile
    when calibrated; the top track's assay description names *what* regulatory signal moved."""
    if not en:
        return "—"
    delta, _l2, track, desc = en
    assay = desc.split(":")[0] if desc else f"track {track}"
    pct = percentile(delta, ref, alt) if (ref or alt) else None
    pct_s = f" p{pct:.1f}" if pct is not None else ""
    return f"Δ{delta:.2f}{pct_s} {band(delta, ref, alt)} ({assay})"
