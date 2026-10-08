"""S2: penetrance / expressivity context for disease-risk findings (pipeline/panel.py).

A pathogenic genotype is elevated risk, not destiny — penetrance is rarely 100% and
severity/age-of-onset (expressivity) varies even among carriers of the same variant.
Same reduced-external-metric idiom as gnomad_constraint.py: a small curated gene-keyed
TSV, joined at report time, degrading gracefully when a gene isn't curated.

Curation is deliberately conservative — a rough, well-established number where one
exists (e.g. BRCA1 lifetime risk ranges), otherwise qualitative language. Genes absent
from this table are not assumed complete-penetrance; panel.py falls back to a generic
caveat (_GENERIC_PENETRANCE_CAVEAT) rather than implying certainty by omission.

**Where a note carries a number, it says where the number came from.** The published
lifetime-risk ranges for BRCA1/2 and the Lynch genes come overwhelmingly from cohorts
recruited *through clinically affected families*, which selects for whatever else makes
those families high-risk; estimates from unselected populations run materially lower. A
person sequenced by choice, rather than referred through a clinic because of a family
history, belongs to the unselected population: for them the ascertained figure is the wrong
reference class and reads high. Such notes are tagged
`ASCERTAINMENT_TAG` and `panel.py` prints `ASCERTAINMENT_NOTE` beside the table whenever a
tagged gene appears — the number is kept (a labelled range beats no range for a reader who
can use it) but it is never presented as *this person's* risk.
"""

from __future__ import annotations

from . import config
from .util import read_tsv

# Marks a note whose number is a family-ascertained estimate (see the module docstring).
# Lives here rather than in panel.py because the table is what carries it.
ASCERTAINMENT_TAG = "[family-ascertained]"

ASCERTAINMENT_NOTE = (
    "_Ranges marked **[family-ascertained]** come from cohorts recruited through clinically "
    "affected families. That selection inflates them relative to an unselected population: "
    "estimates from population cohorts are materially lower. Unless this person was referred "
    "through a clinic because of their family history, treat a tagged range as the upper end "
    "of the published literature, not as this person's risk._")

_TABLE: dict[str, str] | None = None


def available() -> bool:
    return config.PENETRANCE_FILE.exists()


def _load() -> dict[str, str]:
    global _TABLE
    if _TABLE is not None:
        return _TABLE
    _TABLE = {r["gene"].strip(): r["note"].strip()
              for r in read_tsv(config.PENETRANCE_FILE) if r.get("gene")}
    return _TABLE


def note_for(genes: list[str]) -> str | None:
    """Curated note for the first matching gene among `genes`, or None if uncurated.

    A disease can list multiple genes (e.g. Lynch's four MMR genes); genotype variants
    across genes with different notes are rare enough in this panel that "first match"
    is an acceptable simplification over concatenating every gene's note."""
    table = _load()
    for gene in genes:
        if gene in table:
            return table[gene]
    return None


def is_ascertained(note: str | None) -> bool:
    """Does this note quote a family-ascertained figure?"""
    return bool(note) and ASCERTAINMENT_TAG in note


def _reset_cache() -> None:
    global _TABLE
    _TABLE = None
