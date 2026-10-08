"""SpliceAI lookup: deep-learning splice-disruption prediction (Illumina), GRCh38.

A tabix point-query over the precomputed delta-score VCFs — a splice-impact signal
orthogonal to AlphaMissense's missense pathogenicity. It scores the variants the missense
predictors structurally can't speak to: splice-region, synonymous-but-splicing, and deep
intronic changes that create or destroy a splice site. Returns None for anything not in the
table (absent positions, contigs outside gene bodies). Best-effort: callers degrade
gracefully when the large, optional VCFs aren't installed.

SpliceAI INFO field: SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL.
The four DS_* are delta scores in [0,1] for acceptor/donor gain/loss; the variant's overall
delta is their max, with the recommended interpretation cut-offs in config.
"""

from __future__ import annotations

import subprocess

from . import config

# The four delta-score components in INFO order, with a human label for the maximal one.
_EFFECTS = ("acceptor gain", "acceptor loss", "donor gain", "donor loss")

# Which precomputed VCF naming the installed files use (Ensembl "7" vs UCSC "chr7"). The
# GRCh38 SpliceAI VCFs are Ensembl-named, but we probe once and cache so a re-released
# chr-prefixed build still works without per-call double queries. Cached PER FILE: the SNV
# and indel tables are separate downloads that a re-release could name differently, and one
# shared cache would pin the other file's convention onto both.
_prefix: dict[str, str] = {}


def available() -> bool:
    return config.SPLICEAI_SNV_FILE.exists() or config.SPLICEAI_INDEL_FILE.exists()


def _file_for(ref: str, alt: str):
    """Route SNVs to the SNV table and everything else to the indel table.

    Decided on the NORMALIZED alleles, the same ones `lookup` builds its key from. On the raw
    alleles a padded record like 100 AT>GT looks like an indel and was sent to the indel
    table, while _norm trims it to the SNV 100 A>G that only ever lives in the SNV table —
    so the lookup queried the wrong file and silently found nothing."""
    _, nref, nalt = _norm(0, ref, alt)
    snv = len(nref) == 1 and len(nalt) == 1 and nref != "-" and nalt != "-"
    f = config.SPLICEAI_SNV_FILE if snv else config.SPLICEAI_INDEL_FILE
    return f if f.exists() else None


def _norm(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Minimal representation with VEP-style '-' for empty alleles, so a padded VCF
    record (100 A/AT) compares equal to a VEP allele (101 -/T) — without this, indel
    lookups never match (VEP's Uploaded_variation uses '-', the precomputed VCFs pad)."""
    ref, alt = ("" if a == "-" else a for a in (ref, alt))
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return pos, ref or "-", alt or "-"


def _query(path, chrom, lo, hi) -> list[str]:
    """tabix the region, resolving Ensembl-vs-UCSC contig naming once and caching it.
    An empty result for a known contig is legitimate (no SpliceAI record there); only an
    'unknown chromosome' error triggers a naming retry."""
    known = _prefix.get(str(path))
    prefixes = [known] if known is not None else ["", "chr"]
    for pref in prefixes:
        try:
            r = subprocess.run(["tabix", str(path), f"{pref}{chrom}:{lo}-{hi}"],
                               capture_output=True, text=True, timeout=30)
        except Exception:  # noqa: BLE001 — best-effort enrichment, never fatal
            return []
        if "unknown" in r.stderr.lower() or r.returncode != 0:
            continue  # contig naming mismatch — try the alternate prefix
        lines = r.stdout.splitlines()
        # Only lock the naming in on a query that actually RETURNED something. Some tabix
        # builds exit 0 with empty output for an unknown contig instead of erroring, and
        # locking on that empty result would cache the wrong prefix for the whole process
        # — every later lookup silently missing. An empty-but-successful query is also
        # perfectly normal (no SpliceAI record there), so it simply proves nothing.
        if lines and known is None:
            _prefix[str(path)] = pref
            known = pref
        if lines or known is not None:
            return lines
    return []


def lookup(chrom: str, pos, ref: str, alt: str):
    """(delta: float, effect: str, symbol: str) for the maximal splice effect, or None."""
    path = _file_for(ref, alt) if available() else None
    if path is None:
        return None
    pos = int(pos)
    want = _norm(pos, ref, alt)
    # Indel coordinates shift by one between the padded-VCF and VEP conventions, so
    # query a small window for them; SNVs stay an exact point query.
    lo, hi = (pos - 1, pos + 1) if "-" in (ref, alt) else (pos, pos)
    best = None
    for line in _query(path, chrom, lo, hi):
        f = line.split("\t")
        if len(f) < 8:
            continue
        try:
            if _norm(int(f[1]), f[3], f[4]) != want:
                continue
        except ValueError:
            continue
        for ann in _spliceai_anns(f[7]):
            parts = ann.split("|")
            if len(parts) < 6:
                continue
            symbol = parts[1]
            try:
                ds = [float(x) for x in parts[2:6]]
            except ValueError:
                continue
            top = max(range(4), key=lambda i: ds[i])
            if best is None or ds[top] > best[0]:  # strongest effect across genes/transcripts
                best = (ds[top], _EFFECTS[top], symbol)
    return best


def _spliceai_anns(info: str):
    for field in info.split(";"):
        if field.startswith("SpliceAI="):
            yield from field[len("SpliceAI="):].split(",")


def _band(delta: float) -> str:
    if delta >= config.SPLICEAI_DS_HIGH:
        return "high"
    if delta >= config.SPLICEAI_DS_RECOMMENDED:
        return "likely"
    if delta >= config.SPLICEAI_DS_LOW:
        return "possible"
    return "no effect"


def is_pathogenic(sa) -> bool:
    """Splice-altering at the recommended cut-off — feeds ACMG PP3."""
    return bool(sa) and sa[0] >= config.SPLICEAI_DS_RECOMMENDED


def is_benign(sa) -> bool:
    """Confidently no predicted splice effect — feeds ACMG BP4."""
    return bool(sa) and sa[0] < config.SPLICEAI_DS_LOW


def label(sa) -> str:
    """Human-readable cell, e.g. 'donor loss 0.93 (high)' or '—'."""
    if not sa:
        return "—"
    return f"{sa[1]} {sa[0]:.2f} ({_band(sa[0])})"
