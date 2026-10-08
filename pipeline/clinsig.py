"""Reading a ClinVar CLNSIG string.

The value is not one of a handful of words. ClinVar aggregates submissions, so a record's
classification is a small expression:

    Pathogenic/Likely_pathogenic                      several germline classifications
    Pathogenic,_low_penetrance                        a classification with a qualifier
    Pathogenic/Pathogenic,_low_penetrance|risk_factor ...plus an independent assertion

`/` joins classifications, `|` appends assertions of another kind (risk_factor,
drug_response, association, other, ...), and `,_low_penetrance` qualifies the term it
follows. Matching the whole string against a short list of exact values therefore drops
every compound record: about 370 in the 2026-09-28 release, HFE C282Y among them — few
records, but some of the commonest clinically relevant variants there are.

A record is *actionable* when any of its terms is (Likely) pathogenic, with or without the
low-penetrance qualifier. The qualifier is carried through to the reader (`flags`,
`cell`) rather than used to hide the record: low penetrance is a reason to read the finding
differently, not a reason to withhold it.

`Conflicting_classifications_of_pathogenicity` is one term, not a compound of
"pathogenic", so it never matches — the split is on the separators, never a substring test.
"""

from __future__ import annotations

import re

from . import config

LOW_PENETRANCE = ",_low_penetrance"

_SEP = re.compile(r"[/|]")
# One expression, used by Python and (as RE2) by DuckDB, so the SQL filters and the triage
# code cannot drift apart about what counts.
ACTIONABLE_RE = ("(^|[/|])(" + "|".join(sorted(config.ACTIONABLE_TERMS)) + ")("
                 + LOW_PENETRANCE + ")?($|[/|])")
_ACTIONABLE = re.compile(ACTIONABLE_RE)


def terms(sig: str | None) -> list[str]:
    """The record's individual terms, qualifiers still attached."""
    return [t for t in _SEP.split(sig or "") if t]


def is_actionable(sig: str | None) -> bool:
    return bool(sig) and bool(_ACTIONABLE.search(sig))


def sql_actionable(col: str = "clinvar_sig") -> str:
    """SQL predicate equivalent to `is_actionable(col)` (false for NULL)."""
    return f"regexp_matches({col}, '{ACTIONABLE_RE}')"


def is_uncertain(sig: str | None) -> bool:
    """No (likely) pathogenic term, and either nothing at all or an explicitly uncertain one.

    `Uncertain_significance|risk_factor` is as uncertain as the bare value; a transition out
    of it into an actionable classification is the same highest-priority event."""
    if is_actionable(sig):
        return False
    ts = terms(sig)
    return not ts or any(t in config.UNCERTAIN_SIG for t in ts)


def flags(sig: str | None) -> list[str]:
    """What qualifies the classification, in plain words, for the reader.

    'low penetrance' when any actionable term carries the qualifier — marked 'some
    submitters' when an unqualified (Likely) pathogenic term sits beside it — followed by
    every term that is not itself a pathogenic classification (risk_factor, ...)."""
    out: list[str] = []
    ts = terms(sig)
    low = [t for t in ts if t.endswith(LOW_PENETRANCE)]
    if low:
        plain = [t for t in ts if t in config.ACTIONABLE_TERMS]
        out.append("low penetrance (some submitters)" if plain else "low penetrance")
    for t in ts:
        base = t.removesuffix(LOW_PENETRANCE)
        if base not in config.ACTIONABLE_TERMS:
            word = t.replace("_", " ").replace(" ,", ",")
            if word not in out:
                out.append(word)
    return out


def rank(sig: str | None) -> int:
    """Sort order within the actionable set: unqualified Pathogenic first, then anything
    containing a Pathogenic term, then the rest."""
    ts = terms(sig)
    if ts and all(t == "Pathogenic" for t in ts):
        return 0
    if any(t.removesuffix(LOW_PENETRANCE) == "Pathogenic" for t in ts):
        return 1
    return 2


def cell(sig: str | None) -> str:
    """The classification as a markdown table cell.

    `|` is the column separator in every report table and in the parsers that read them
    back (render, dashboard), so a raw compound value would shift every later column.
    Assertions are joined with `; ` instead, and any qualifier is spelled out after it."""
    text = (sig or "").replace("|", "; ")
    notes = [f for f in flags(sig) if f.startswith("low penetrance")]
    return f"{text} ⚠️ {notes[0]}" if notes else text
