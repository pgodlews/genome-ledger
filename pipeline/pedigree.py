"""Family relationship graph + perspective engine.

`family.tsv` (config.FAMILY_FILE) is the human-editable source of truth: one row per person —
`id, display_name, sex, father, mother, partner` (the parent/partner cells reference other
rows' ids). Sample ids are perspective-neutral given names, so relationship terms like
"father"/"mother" are **not** identities but are computed *relative to a chosen ego*: the same
person is "father" from one view and "grandfather" from another. `relationship_to(ego, other)`
does that computation; "switching perspective" is just changing the ego.

One person can be sequenced more than once. This repo's id convention is that a bare id
(`Adam`) is our own calls and a suffixed one (`Adam_vendor`, `Jan_t2t`) is an
ALTERNATE CALLSET OF THE SAME PERSON — never a second person. `person_of()` is the graph's
copy of that rule (`denovo.parent_callsets()` and `segregation._parent_calls()` already read
ids this way), so `family.tsv` holds one row per *person* and every relationship, dashboard
and PED row resolves a callset id back to its owner first. Without it the vendor rows became
a duplicate family tree: a person became their own sibling and their children their nieces.

`pedigree.ped` (the standard PED consumed by Mendelian/trio tooling) is *generated* from this
graph — `trio()` + `generate_pedigree_ped()` are the hooks the M2 segregation queries will use
once a child is sequenced.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from . import config
from .util import log, read_tsv

_SEX_CODE = {"male": "1", "m": "1", "female": "2", "f": "2"}


@dataclass
class Person:
    id: str
    display_name: str
    sex: str                  # normalized "male" / "female" / "" (unknown)
    father: str | None
    mother: str | None
    partner: str | None


class Family:
    """An in-memory relationship graph. Construct via `load()` or directly (for tests)."""

    def __init__(self, people: dict[str, Person], sequenced: set[str] | None = None):
        self.people = people
        self._sequenced = sequenced
        # Reverse index: parent id -> [child ids], in insertion order.
        self._children: dict[str, list[str]] = {pid: [] for pid in people}
        for p in people.values():
            for parent in (p.father, p.mother):
                if parent in self._children:
                    self._children[parent].append(p.id)
        # Partners are symmetric: honour an edge declared on either row.
        self._partner: dict[str, str] = {}
        for p in people.values():
            if p.partner and p.partner in people:
                self._partner[p.id] = p.partner
                self._partner.setdefault(p.partner, p.id)

    # --- callset <-> person ---------------------------------------------------
    def person_of(self, sample: str) -> str:
        """The person a callset id belongs to: `Adam_vendor` -> `Adam`.

        A row in `family.tsv` is always its own person, so an id that IS a row is returned
        untouched (a real person whose id happens to contain an underscore is never eaten).
        Otherwise the longest prefix ending at an underscore that names a person wins, which
        also handles a stacked suffix (`Adam_vendor_rerun`). Unknown ids pass through, so
        callers can treat this as a no-op for ids the graph has never heard of."""
        if sample in self.people:
            return sample
        parts = sample.split("_")
        for i in range(len(parts) - 1, 0, -1):
            head = "_".join(parts[:i])
            if head in self.people:
                return head
        return sample

    def _seq(self) -> set[str]:
        if self._sequenced is None:
            self._sequenced = sequenced_ids()
        return self._sequenced

    def is_ingested(self, sample: str) -> bool:
        """True for an exact CALLSET id present in samples.tsv (not person-level)."""
        return sample in self._seq()

    def callsets(self, pid: str) -> list[str]:
        """Every ingested callset belonging to `pid`'s person, bare id first.

        The bare id leads because it is this repo's own GRCh38 calls — the primary view —
        and the suffixed ones are alternates shown beside it (`add, never replace`)."""
        person = self.person_of(pid)
        own = [s for s in self._seq() if self.person_of(s) == person]
        return sorted(own, key=lambda s: (s != person, s))

    # --- graph queries -------------------------------------------------------
    def parents(self, pid: str) -> list[str]:
        p = self.people.get(self.person_of(pid))
        if not p:
            return []
        return [x for x in (p.father, p.mother) if x in self.people]

    def children(self, pid: str) -> list[str]:
        return list(self._children.get(self.person_of(pid), []))

    def partner(self, pid: str) -> str | None:
        return self._partner.get(self.person_of(pid))

    def siblings(self, pid: str) -> list[str]:
        """Anyone sharing at least one parent with pid (excluding pid itself)."""
        pid = self.person_of(pid)
        mine = set(self.parents(pid))
        if not mine:
            return []
        return [q for q in self.people
                if q != pid and mine.intersection(self.parents(q))]

    def is_founder(self, pid: str) -> bool:
        return not self.parents(pid)

    def trio(self, child: str) -> tuple[str, str, str] | None:
        """(child, father, mother) when both parents are known — the M2/segregation hook.

        `child` comes back exactly as passed in, so a CALLSET id keeps its identity
        (`trio("Lena_vendor")` still resolves) while the parents resolve to people;
        callers expand a parent into its callsets themselves (`denovo.parent_callsets`)."""
        p = self.people.get(self.person_of(child))
        if p and p.father in self.people and p.mother in self.people:
            return (child, p.father, p.mother)
        return None

    def is_sequenced(self, pid: str) -> bool:
        """True when the PERSON has at least one ingested callset, under any suffix."""
        person = self.person_of(pid)
        return any(self.person_of(s) == person for s in self._seq())

    # --- perspective ---------------------------------------------------------
    def relationship_to(self, ego: str, other: str) -> str:
        """Kinship term for `other` as seen from `ego` (e.g. 'father', 'mother-in-law').

        Both ids resolve to their person first, so an alternate callset is never a relative
        of the person it belongs to — `relationship_to("Adam", "Adam_vendor")` is
        "self", not "brother"."""
        ego, other = self.person_of(ego), self.person_of(other)
        if ego == other:
            return "self"
        e, o = self.people.get(ego), self.people.get(other)
        if not e or not o:
            return "relative"

        def by_sex(male: str, female: str, neutral: str) -> str:
            return {"male": male, "female": female}.get(o.sex, neutral)

        if other in self.parents(ego):
            return by_sex("father", "mother", "parent")
        if other in self.children(ego):
            return by_sex("son", "daughter", "child")
        if other == self.partner(ego):
            return "partner"
        if other in self.siblings(ego):
            return by_sex("brother", "sister", "sibling")

        grandparents = {gp for p in self.parents(ego) for gp in self.parents(p)}
        if other in grandparents:
            return by_sex("grandfather", "grandmother", "grandparent")
        grandchildren = {gc for c in self.children(ego) for gc in self.children(c)}
        if other in grandchildren:
            return by_sex("grandson", "granddaughter", "grandchild")

        partner = self.partner(ego)
        if partner:
            if other in self.parents(partner):
                return by_sex("father-in-law", "mother-in-law", "parent-in-law")
            if other in self.siblings(partner):
                return by_sex("brother-in-law", "sister-in-law", "sibling-in-law")
        # child's partner
        if any(other == self.partner(c) for c in self.children(ego)):
            return by_sex("son-in-law", "daughter-in-law", "child-in-law")
        # sibling's partner
        if any(other == self.partner(s) for s in self.siblings(ego)):
            return by_sex("brother-in-law", "sister-in-law", "sibling-in-law")

        parent_sibs = {s for p in self.parents(ego) for s in self.siblings(p)}
        if other in parent_sibs:
            return by_sex("uncle", "aunt", "relative")
        niblings = {n for s in self.siblings(ego) for n in self.children(s)}
        if other in niblings:
            return by_sex("nephew", "niece", "relative")
        cousins = {c for pib in parent_sibs for c in self.children(pib)}
        if other in cousins:
            return "cousin"
        return "relative"

    def others(self, ego: str) -> list[tuple[str, str]]:
        """[(id, relationship_to_ego)] for every other PERSON, in file order. Alternate
        callsets are not listed here — they belong to a person, not beside them."""
        ego = self.person_of(ego)
        return [(pid, self.relationship_to(ego, pid))
                for pid in self.people if pid != ego]

    # --- validation ----------------------------------------------------------
    def validate(self) -> list[str]:
        """Non-fatal consistency warnings (dangling refs, sex/parent mismatch, cycles)."""
        issues = []
        # A row whose id resolves to ANOTHER row is an alternate callset that was written
        # into the graph as a second person — the bug that makes someone their own sibling and
        # their children their nieces. Callsets belong in samples.tsv only.
        for pid in self.people:
            owner = Family({k: v for k, v in self.people.items() if k != pid}).person_of(pid)
            if owner != pid:
                issues.append(f"{pid}: looks like a callset of '{owner}' — remove this row "
                              f"from family.tsv (callsets live in samples.tsv)")
        for p in self.people.values():
            for role in ("father", "mother", "partner"):
                ref = getattr(p, role)
                if ref and ref not in self.people:
                    issues.append(f"{p.id}: {role} '{ref}' is not a known person")
            if p.father in self.people and self.people[p.father].sex == "female":
                issues.append(f"{p.id}: father '{p.father}' is recorded female")
            if p.mother in self.people and self.people[p.mother].sex == "male":
                issues.append(f"{p.id}: mother '{p.mother}' is recorded male")
        for pid in self.people:                      # ancestor cycle check
            seen, stack = set(), [pid]
            while stack:
                cur = stack.pop()
                for par in self.parents(cur):
                    if par == pid:
                        issues.append(f"{pid}: is its own ancestor (cycle)")
                        break
                    if par not in seen:
                        seen.add(par)
                        stack.append(par)
        return issues


# --- module helpers ----------------------------------------------------------
def _norm(v: str | None) -> str | None:
    v = (v or "").strip()
    return v or None


def load(path: Path | None = None, sequenced: set[str] | None = None) -> Family:
    people: dict[str, Person] = {}
    for r in read_tsv(path or config.FAMILY_FILE):
        pid = (r.get("id") or "").strip()
        if not pid:
            continue
        people[pid] = Person(
            id=pid,
            display_name=(r.get("display_name") or "").strip() or pid,
            sex=(r.get("sex") or "").strip().lower(),
            father=_norm(r.get("father")),
            mother=_norm(r.get("mother")),
            partner=_norm(r.get("partner")),
        )
    return Family(people, sequenced=sequenced)


def sequenced_ids() -> set[str]:
    return {r["sample_id"] for r in read_tsv(config.SAMPLES_TSV) if r.get("sample_id")}


def generate_pedigree_ped(family: Family | None = None, path: Path | None = None) -> Path:
    """Write a standard PED from the family graph (FamilyID Individual Paternal Maternal Sex
    Pheno), one row per ingested callset. Parent columns reference ids or '0'; replaces
    ingest's stub writer and gives the M2 trio/segregation queries real parent links."""
    family = family or load()
    path = path or config.PEDIGREE_FILE
    for issue in family.validate():
        log.warning("pedigree: %s", issue)

    def _row(sample: str, p: Person) -> str:
        # The suffix that makes this callset ("" for the bare one). A callset's parent
        # column points at the parent callset carrying the SAME suffix when one was
        # ingested (a T2T child against T2T parents), and falls back to the parent's bare
        # id otherwise — PED tooling keys on the sample id as it appears in the VCF.
        suffix = sample[len(p.id):]

        def parent(par: str | None) -> str:
            if not par:
                return "0"
            return par + suffix if family.is_ingested(par + suffix) else par

        return "\t".join(["FAM1", sample, parent(p.father), parent(p.mother),
                          _SEX_CODE.get(p.sex, "0"), "0"])

    # One row per CALLSET, not per person: `Jan` and `Jan_t2t` are one genome called
    # twice, but each carries its own sample id in its own VCF and needs its own PED row.
    # A person with nothing ingested yet still gets their single placeholder row.
    lines = [_row(s, p) for p in family.people.values()
             for s in (family.callsets(p.id) or [p.id])]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path
