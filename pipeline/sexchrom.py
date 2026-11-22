"""How a position is inherited, and how many copies of it a person has.

Autosomes and the pseudoautosomal regions (PAR) are diploid with one copy from each parent.
Outside the PAR a male has one X, from his mother, and one Y, from his father, and everyone's
mitochondria come from their mother. Callers still write these sites as ordinary diploid
genotypes (`config.py` leaves X and Y diploid), so a male's single X copy can be stored as HOM.
Anything that says "two copies" or names a parent has to read the locus through this module.
"""

from __future__ import annotations

from . import config

AUTOSOMAL, X_LINKED, Y_LINKED, MITOCHONDRIAL = "autosomal", "x", "y", "mt"

# 1-based inclusive PAR1 and PAR2 coordinates.
_PAR = {
    "GRCh38": {"X": ((10001, 2781479), (155701383, 156030895)),
               "Y": ((10001, 2781479), (56887903, 57217415))},
    "GRCh37": {"X": ((60001, 2699520), (154931044, 155260560)),
               "Y": ((10001, 2649520), (59034050, 59363566))},
}


def contig(chrom) -> str:
    """'chrX' -> 'X', 'chrM' / 'M' -> 'MT'."""
    c = str(chrom)
    c = c[3:] if c.lower().startswith("chr") else c
    return "MT" if c.upper() in ("M", "MT") else c.upper()


def locus_kind(chrom, pos: int) -> str:
    """AUTOSOMAL (which includes the PAR), X_LINKED, Y_LINKED or MITOCHONDRIAL."""
    c = contig(chrom)
    if c == "MT":
        return MITOCHONDRIAL
    if c not in ("X", "Y"):
        return AUTOSOMAL
    if any(lo <= int(pos) <= hi for lo, hi in _PAR.get(config.BUILD, {}).get(c, ())):
        return AUTOSOMAL
    return X_LINKED if c == "X" else Y_LINKED


def is_male(sex: str | None) -> bool:
    return (sex or "").strip().lower() in ("male", "m", "1")


def is_female(sex: str | None) -> bool:
    return (sex or "").strip().lower() in ("female", "f", "2")


def is_alternate_mt(sample: str, chrom, fam) -> bool:
    """A mitochondrial row from an alternate callset (`Adam_vendor`, `Adam_t2t`).

    Build-to-build chrM liftover is unreliable: a small offset can land a record on a ClinVar
    pathogenic allele. The `mito` stage calls chrM from the reads and is the authority there."""
    return contig(chrom) == "MT" and fam.person_of(sample) != sample
