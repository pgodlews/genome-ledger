"""gnomAD gene-constraint lookup: loss-of-function intolerance (LOEUF / pLI), GRCh38.

A small gene-keyed table (gene → pLI, LOEUF) reduced by `setup` from gnomAD v4.1's
per-transcript constraint metrics. It answers one question the ACMG engine needs: is a
gene *loss-of-function-intolerant* — i.e. is haploinsufficiency a plausible disease
mechanism? That gates how much weight a null variant's PVS1 criterion carries.

  - **pLI** — probability a gene is intolerant of a single LOF allele (haploinsufficient).
    pLI ≥ 0.9 is the classic "LOF-constrained" call.
  - **LOEUF** — upper bound of the observed/expected LOF ratio; *lower* = more constrained.
    LOEUF < ~0.6 flags genes pLI alone can miss.

Keyed by gene symbol, so it's genome-build-agnostic. Best-effort: callers degrade
gracefully (treat as "no data") when the optional table isn't installed or a gene is
absent. Loaded once and memoised.
"""

from __future__ import annotations

from . import config

_TABLE: dict[str, tuple[float | None, float | None]] | None = None


def available() -> bool:
    return config.GNOMAD_CONSTRAINT_FILE.exists()


def _load() -> dict[str, tuple[float | None, float | None]]:
    """gene → (pLI, LOEUF) from the reduced table; memoised. Missing metrics are None."""
    global _TABLE
    if _TABLE is not None:
        return _TABLE
    table: dict[str, tuple[float | None, float | None]] = {}
    if available():
        for line in config.GNOMAD_CONSTRAINT_FILE.read_text().splitlines():
            if not line or line.startswith("gene\t") or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 3:
                continue
            table[f[0]] = (_num(f[1]), _num(f[2]))
    _TABLE = table
    return table


def _num(s: str) -> float | None:
    try:
        return float(s)
    except (ValueError, TypeError):
        return None


def lookup(gene: str) -> tuple[float | None, float | None] | None:
    """(pLI, LOEUF) for a gene, or None when the gene isn't in the table."""
    return _load().get(gene)


def is_constrained(gene: str) -> bool | None:
    """True/False if `gene` is LOF-constrained; None when there's no usable metric.

    Constrained = pLI ≥ GNOMAD_PLI_CONSTRAINED or LOEUF < GNOMAD_LOEUF_CONSTRAINED. None
    (no data) is deliberately distinct from False so the ACMG engine can decline to penalize
    a gene we simply have no constraint metric for."""
    rec = lookup(gene)
    if rec is None:
        return None
    pli, loeuf = rec
    if pli is None and loeuf is None:
        return None
    return ((pli is not None and pli >= config.GNOMAD_PLI_CONSTRAINED)
            or (loeuf is not None and loeuf < config.GNOMAD_LOEUF_CONSTRAINED))


def label(gene: str) -> str:
    """Report cell, e.g. 'pLI 0.99 / LOEUF 0.21' or '—' when unavailable."""
    rec = lookup(gene)
    if rec is None:
        return "—"
    pli, loeuf = rec
    parts = []
    if pli is not None:
        parts.append(f"pLI {pli:.2f}")
    if loeuf is not None:
        parts.append(f"LOEUF {loeuf:.2f}")
    return " / ".join(parts) or "—"


def _reset_cache() -> None:
    """Drop the memoised table (tests/setup that rewrite the file)."""
    global _TABLE
    _TABLE = None
