"""AlphaMissense lookup: deep-learning missense pathogenicity (DeepMind), GRCh38.

A tabix point-query over the precomputed table — a computational pathogenicity signal
independent of ClinVar, useful for corroborating findings and (in future) triaging
variants of uncertain significance. Missense SNVs only; returns None for anything not
covered (indels, non-coding, absent). Best-effort: callers degrade gracefully when the
large, optional table isn't installed.
"""

from __future__ import annotations

import subprocess

from . import config

_CLASS = {"likely_pathogenic": "likely path.", "ambiguous": "ambiguous",
          "likely_benign": "likely benign"}


def available() -> bool:
    return config.ALPHAMISSENSE_FILE.exists()


def lookup(chrom: str, pos, ref: str, alt: str):
    """(am_pathogenicity: float, am_class: str) for a missense SNV, or None.
    AlphaMissense uses UCSC (chr1) contigs; our variants are Ensembl-named (1)."""
    if not available():
        return None
    try:
        out = subprocess.run(
            ["tabix", str(config.ALPHAMISSENSE_FILE), f"chr{chrom}:{pos}-{pos}"],
            capture_output=True, text=True, timeout=30).stdout
    except Exception:  # noqa: BLE001 — best-effort enrichment, never fatal
        return None
    best = None
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 10 or f[2] != ref or f[3] != alt:
            continue
        try:
            score = float(f[8])
        except ValueError:
            continue
        if best is None or score > best[0]:  # most pathogenic across transcripts
            best = (score, f[9])
    return best


def label(am) -> str:
    """Human-readable cell, e.g. 'likely path. (0.97)' or '—'."""
    if not am:
        return "—"
    return f"{_CLASS.get(am[1], am[1])} ({am[0]:.2f})"
