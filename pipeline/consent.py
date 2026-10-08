"""S1: family opt-in filter for incidental findings.

"Incidental" here means adult-onset, currently-untreatable predictions where learning the
result changes nothing a clinician can act on (e.g. Huntington disease repeat sizing) — as
opposed to the panel's actionable findings (BRCA1/2, Lynch, etc.), which are never gated:
knowing enables surveillance/intervention, so the disclosure trade-off is different.

Consent is per (PERSON, category), recorded in `family/consent.tsv` (config.CONSENT_FILE,
same gitignored-data pattern as family.tsv — one row per person, hand-edited or via
`run.py consent`). Default is **not consented**: an incidental finding is withheld from the
rendered report — replaced by a note that it exists and how to opt in — until a row grants
`consent=yes` for that person+category. This is a disclosure gate, not a data gate:
everything is still computed and stored; only report rendering asks first.

**A grant belongs to the person, not to a callset.** It used to key on the exact
`sample_id`, which meant a second callset of someone who had opted in started out withheld
again: `Jan` was disclosed while `Jan_t2t` — the same genome, called against T2T —
was not, so one person's two reports disagreed about what he was allowed to be told. Ids
now resolve through `pedigree.person_of()` exactly as relationships and dashboards do, and
a revocation reaches the gated outputs of *every* callset of that person.
"""

from __future__ import annotations

import csv
import os

from . import config
from .util import log, read_tsv

# The one incidental category wired up so far — repeat-expansion neurodegenerative
# disorders (pipeline/repeats.py) with no disease-modifying treatment.
CATEGORY_ADULT_ONSET_UNTREATABLE = "adult_onset_untreatable"

_TRUE = {"yes", "true", "1"}


def gated_outputs(sample: str, category: str) -> list:
    """Files whose CONTENT depends on this (sample, category) consent decision.

    Revocation has to reach these. Updating the ledger alone left the previously generated
    report sitting on disk with the findings still in it — readable by the dashboard, by
    `render`, and by anything else that reads the reports tree — so the disclosure continued after consent for it was withdrawn.
    Withdrawing consent is not a preference for next time; it is a statement about what may
    be shown now.
    """
    from . import config as _cfg
    if category != CATEGORY_ADULT_ONSET_UNTREATABLE:
        return []
    from . import reportpaths
    from .render import FORMATS

    def every_format(md) -> list:
        # The report in EVERY format tree: a rendering left behind keeps disclosing.
        return [md, *[reportpaths.derived(md, name, ext)
                      for name, (ext, _r) in FORMATS.items()]]

    out = [*every_format(reportpaths.genome_report("repeats", sample)),
           _cfg.REPEATS_DIR / sample / f"repeats_{sample}.consent-stamp"]
    # The PRS report is gated by the same decision (prs._panel_filters), so it is the same
    # kind of output: a report scored while consent stood lists the severe, non-actionable
    # traits, and it kept listing them — to the dashboard as well — after consent was
    # withdrawn. Only while the S1 gate is on: with PRS_S1_GATE off nothing in that report
    # depends on consent, and removing it would discard a scoring run for no change.
    if _cfg.PRS_S1_GATE:
        out += [*every_format(reportpaths.genome_report("prs", sample)),
                _cfg.PRS_DIR / sample / f"prs_{sample}.consent-stamp",
                _cfg.PRS_DIR / sample / f"prs_{sample}.panel-stamp"]
    return out


def _person(sample: str) -> str:
    """The person a sample id belongs to, or the id itself if the graph cannot say.

    Deliberately tolerant: consent must keep working when family.tsv is missing or
    unreadable, and the fallback (treat the id as its own person) is the pre-existing
    behaviour, so a broken graph narrows nothing and widens nothing."""
    try:
        from . import pedigree
        return pedigree.load().person_of(sample)
    except Exception as e:  # noqa: BLE001 — a disclosure gate must not fail open OR closed
        log.warning("consent: could not resolve '%s' to a person (%s); using the id as-is",
                    sample, e)
        return sample


def _load() -> dict[tuple[str, str], bool]:
    """{(person, category): granted}, resolving every row's id to its owner.

    A row written against a callset id (`Adam_vendor`) is read as that person's grant,
    so ledgers predating the person-keyed rule keep working. An explicit row for the person
    always wins over one written against a callset of theirs — otherwise a legacy row could
    silently override the decision actually recorded for the person."""
    by_person: dict[tuple[str, str], bool] = {}
    from_callset: dict[tuple[str, str], bool] = {}
    for row in read_tsv(config.CONSENT_FILE):
        sample, category = row.get("sample_id", ""), row.get("category", "")
        if not sample or not category:
            continue
        granted = row.get("consent", "").strip().lower() in _TRUE
        person = _person(sample)
        target = by_person if person == sample else from_callset
        target[(person, category)] = granted
    for key, granted in from_callset.items():
        if key in by_person and by_person[key] != granted:
            log.warning("consent: %s/%s is recorded for the person and differently for one "
                        "of their callsets; the person's own row wins.", *key)
        by_person.setdefault(key, granted)
    return by_person


def has_consented(sample: str, category: str) -> bool:
    """Default False — absence of a row means the family member has not opted in.

    `sample` may be any callset of theirs; the grant is looked up for the person."""
    return _load().get((_person(sample), category), False)


_FIELDS = ["sample_id", "category", "consent"]


def callsets_of(sample: str) -> list[str]:
    """Every ingested callset of this person, or just the id if the graph cannot say."""
    try:
        from . import pedigree
        fam = pedigree.load(sequenced=pedigree.sequenced_ids())
        return fam.callsets(sample) or [sample]
    except Exception:  # noqa: BLE001 — see _person
        return [sample]


def set_consent(sample: str, category: str, consent: bool) -> None:
    """Upsert the (PERSON, category) grant/revocation into consent.tsv.

    The row is written under the person `sample` belongs to, and any rows recorded against
    a callset of theirs for this category are folded into it — one decision, one row, so
    the ledger cannot hold two answers for one person. A revocation removes the gated
    outputs of *every* callset: the report that discloses a finding is per callset, and
    leaving the others on disk would keep disclosing what was just withdrawn.

    Written to a temp file and renamed over the original. Opening the real file "w" truncated
    it *before* the rows were validated, so a hand-edited file with an extra column made
    DictWriter raise mid-write and left the disclosure gate's ledger empty — every consent
    grant lost, and the S1 gate silently back to withholding everything. `extrasaction`
    ignores columns this writer doesn't know rather than raising on them, and the rename is
    atomic, so a crash anywhere leaves the previous file intact. Columns the file already
    carries that this writer doesn't manage (a hand-added `note`, say) are preserved rather
    than dropped — losing them quietly would be a smaller version of the same bug.
    """
    person = _person(sample)
    owned = set(callsets_of(person)) | {person, sample}
    rows = read_tsv(config.CONSENT_FILE)
    rows = [r for r in rows
            if not (r.get("sample_id") in owned and r.get("category") == category)]
    rows.append({"sample_id": person, "category": category,
                 "consent": "yes" if consent else "no"})
    rows.sort(key=lambda r: (r.get("sample_id") or "", r.get("category") or ""))

    config.CONSENT_FILE.parent.mkdir(parents=True, exist_ok=True)
    extra = [k for k in dict.fromkeys(k for r in rows for k in r) if k not in _FIELDS]
    tmp = config.CONSENT_FILE.with_name(config.CONSENT_FILE.name + ".part")
    try:
        with open(tmp, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=_FIELDS + extra, delimiter="\t",
                               extrasaction="ignore", restval="")
            w.writeheader()
            w.writerows(rows)
        os.replace(tmp, config.CONSENT_FILE)
        if not consent:
            removed = []
            for callset in sorted(owned):
                for out in gated_outputs(callset, category):
                    if out.exists():
                        out.unlink()
                        removed.append(out.name)
            if removed:
                log.warning("Revoked %s for %s — removed the output(s) that disclosed it: "
                            "%s. Re-run the stage to regenerate them under the new policy.",
                            category, person, ", ".join(sorted(removed)))
            else:
                log.info("Revoked %s for %s (no generated output was disclosing it).",
                         category, person)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    log.info("consent: %s/%s -> %s", person, category, "yes" if consent else "no")
