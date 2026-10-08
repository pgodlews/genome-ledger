"""Allele balance: the alt-read fraction behind a reported heterozygous call.

VEP's ZYG column (all the `variants` table keeps) says HET or HOM and nothing about how
many reads said so. A germline het sits near 50% alt reads; a call at 2 of 19 is what
homopolymer slippage and mismapped paralog reads look like, and GATK can still emit and
PASS it. This reads AD back from the
sample's normalized VCF for just the rows a report shows, so the reports can flag such
calls for review. Flag, never drop — see config.LOW_VAF_THRESHOLD.

The `variants` table stores VEP's trimmed allele form ("6:107893551 -/T"), while the VCF
keeps the anchor base ("6:107893550 G/GT"); `vep_key` maps one onto the other. Best-effort:
a missing VCF, a record without AD, or a zero depth yields no entry, and callers show the
call unflagged, exactly as before.
"""

from __future__ import annotations

import os
import subprocess
import tempfile

from . import config
from .normalize import normalized_path


def vep_key(chrom: str, pos: int, ref: str, alt: str) -> tuple[str, int, str, str]:
    """VCF (chrom, pos, ref, alt) → the form VEP writes to Uploaded_variation.

    VEP strips the shared leading base of an indel/complex allele pair and shifts POS by
    one, writing an emptied allele as '-'. SNVs and MNVs of equal length are untouched."""
    if len(ref) != len(alt) and ref[:1] == alt[:1]:
        return chrom, pos + 1, ref[1:] or "-", alt[1:] or "-"
    return chrom, pos, ref, alt


def fractions(sample: str, keys) -> dict[tuple, tuple[int, int]]:
    """{(chrom, pos, ref, alt) in VEP form: (alt_reads, total_reads)} for the given keys.

    One batched bcftools query over the sample's normalized VCF (split multiallelics, so AD
    is exactly ref,alt). Keys not found, or without usable AD, are simply absent."""
    keys = {(str(c), int(p), r, a) for c, p, r, a in keys}
    vcf = normalized_path(sample)
    if not keys or not vcf.exists():
        return {}
    with tempfile.NamedTemporaryFile("w", suffix=".regions", delete=False) as fh:
        # pos-1 covers the anchor base of an indel whose VEP position is shifted by one.
        for c, p, *_ in sorted(keys):
            fh.write(f"{c}\t{max(1, p - 1)}\t{p}\n")
        regions = fh.name
    try:
        out = subprocess.run(
            ["bcftools", "query", "-R", regions, "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n",
             str(vcf)], capture_output=True, text=True, timeout=300).stdout
    except Exception:  # noqa: BLE001 — best-effort enrichment, never fatal
        return {}
    finally:
        os.unlink(regions)
    found: dict[tuple, tuple[int, int]] = {}
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        key = vep_key(f[0], int(f[1]), f[2], f[3])
        if key not in keys:
            continue
        try:
            ad = [int(x) for x in f[4].split(",")]
        except ValueError:
            continue
        if len(ad) >= 2 and sum(ad) > 0:
            found[key] = (ad[1], sum(ad))
    return found


def is_low(zyg: str | None, ab: tuple[int, int] | None) -> bool:
    """True for a HET call whose alt fraction is below config.LOW_VAF_THRESHOLD."""
    if not ab or (zyg or "").upper() != "HET":
        return False
    alt, total = ab
    return alt / total < config.LOW_VAF_THRESHOLD


def label(zyg: str | None, ab: tuple[int, int] | None) -> str:
    """The Zyg cell: 'HET', or 'HET ⚠️ low allele fraction 3/18 (17%) — review'."""
    z = zyg or "-"
    if not is_low(zyg, ab):
        return z
    alt, total = ab
    return f"{z} ⚠️ low allele fraction {alt}/{total} ({alt / total:.0%}) — review"
