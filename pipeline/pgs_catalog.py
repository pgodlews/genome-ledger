"""An offline copy of the PGS Catalog's metadata, and search over it.

The Catalog publishes its whole metadata as five CSV files (scores, performance metrics,
evaluation sample sets, traits, publications: ~12 MB). `sync` downloads them, plus the trait
categories the Catalog's own browser groups traits by, into a dated directory under
config.PGS_CATALOG_DIR. Everything else here reads only that directory, so a machine with no
network can still browse every score: the copy is carried across like a score bundle.

The performance summary is deliberately conservative, after two mistakes a panel review caught:

- **Only covariate-free figures count as the score's own.** A model with age and sex reports an
  AUROC that age earns as much as the score does (age alone separates most late-onset diseases).
  Principal components are ancestry correction, not clinical information, so they are allowed.
  A covariate field recorded as `0`, `-` or `NA` says nothing reliable and is kept apart as
  "unclear", not taken as "none".
- **The median over cohorts, not the best cohort.** Picking the highest of several evaluations
  reports the luckiest one. Cohorts with an impossible case count (as many cases as people) or
  fewer than 100 cases are left out of the summary but still listed by `evaluations`.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import statistics
import urllib.request
from pathlib import Path

from . import config
from .util import log, today

FILES = {
    "scores": "pgs_all_metadata_scores.csv",
    "performance": "pgs_all_metadata_performance_metrics.csv",
    "sample_sets": "pgs_all_metadata_evaluation_sample_sets.csv",
    "traits": "pgs_all_metadata_efo_traits.csv",
    "publications": "pgs_all_metadata_publications.csv",
}
CATEGORIES_FILE = "trait_categories.json"
MANIFEST = "manifest.json"
FORMAT = 1

ANCESTRY = "European"     # the reference the percentiles are computed against (KG_POP=EUR)
MIN_EVAL_N = 1000
MIN_CASES = 100

_LICENCE_TAGS = (         # the panel builder's tags (panels/build_complex_traits.py)
    ("should be cited appropriately", "pgs-catalog-default"),
    ("CC0", "CC0"),
    ("NonCommercial-NoDerivatives", "CC-BY-NC-ND"),
    ("Creative Commons Attribution 4.0", "CC-BY-4.0"),
    ("academic community", "academic-noncommercial"),
    ("academic research purposes only", "academic-noncommercial"),
    ("research purposes only", "research-only"),
)


def licence_tag(text: str) -> str:
    for needle, tag in _LICENCE_TAGS:
        if needle.lower() in (text or "").lower():
            return tag
    return "other"


def _num(x: str | None) -> float | None:
    m = re.match(r"\s*(-?\d+\.?\d*)", x or "")
    return float(m.group(1)) if m else None


# --- covariates --------------------------------------------------------------------------
_UNCLEAR = frozenset({"0", "-", "na", "n/a", "nr", "not reported", "none reported"})
_PC = re.compile(
    r"^(the\s+)?(first\s+|top\s+)?(\d+\s+)?(genetic\s+|genotype\s+|global\s+|ancestry\s+|"
    r"genomic\s+)?(pcs?|principal\s+components?)"
    r"(\s*\(?\s*(pcs?\s*)?\d+\s*([-–]|to)\s*(pcs?\s*)?\d+\s*\)?|\s*\d+)?\.?$", re.I)


def covariate_class(text: str | None) -> str:
    """`none`, `pcs` (ancestry components only), `unclear`, or `adjusted`.

    `none` and `pcs` are the covariate-free evaluations: what the score does on its own."""
    s = (text or "").strip()
    if not s:
        return "none"
    if s.lower() in _UNCLEAR:
        return "unclear"
    return "pcs" if _PC.match(s) else "adjusted"


def genetic_only(text: str | None) -> bool:
    return covariate_class(text) in ("none", "pcs")


# --- the snapshot on disk -----------------------------------------------------------------
def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _get(url: str, timeout: int = 300) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "genome-ledger"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _download(url: str, dest: Path) -> None:
    tmp = dest.with_name(f"{dest.name}.{os.getpid()}.part")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "genome-ledger"})
        with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as fh:
            shutil.copyfileobj(r, fh)
        os.replace(tmp, dest)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def _categories() -> list[dict]:
    """[{label, efotraits: [{id, label}]}], following the REST pagination."""
    out, url = [], config.PGS_REST_URL + "trait_category/all"
    while url:
        page = json.loads(_get(url))
        for c in page.get("results", []):
            out.append({"label": c.get("label", ""),
                        "efotraits": [{"id": t.get("id", ""), "label": t.get("label", "")}
                                      for t in c.get("efotraits", [])]})
        url = page.get("next")
    return out


def _write_manifest(d: Path, release: str, n_scores) -> None:
    files = {name: _sha256(d / name) for name in [*FILES.values(), CATEGORIES_FILE]}
    (d / MANIFEST).write_text(json.dumps(
        {"format": FORMAT, "release": release, "fetched": today(), "scores": n_scores,
         "files": files}, indent=1) + "\n")


def verify(d: Path) -> list[str]:
    """Problems with a snapshot directory; empty when it is complete and intact."""
    try:
        m = json.loads((d / MANIFEST).read_text())
    except (OSError, ValueError):
        return [f"{d}: no readable {MANIFEST}"]
    problems = []
    for name, digest in (m.get("files") or {}).items():
        p = d / name
        if not p.is_file():
            problems.append(f"{name} missing")
        elif _sha256(p) != digest:
            problems.append(f"{name} does not match its checksum")
    if set(m.get("files") or {}) != set(FILES.values()) | {CATEGORIES_FILE}:
        problems.append("the manifest does not list the expected files")
    return problems


def sync(dest_root: Path | None = None, force: bool = False) -> Path:
    """Download the current Catalog metadata into <root>/<release date>/ and return it.

    Needs the network; nothing about the family is sent (these are the Catalog's own public
    files). With no `dest_root` the copy becomes this machine's current snapshot."""
    root = Path(dest_root) if dest_root else config.PGS_CATALOG_DIR
    info = json.loads(_get(config.PGS_REST_URL + "info"))
    latest = info.get("latest_release") or {}
    release = latest.get("date") or today()
    d = root / release
    if d.is_dir() and not verify(d) and not force:
        log.info("PGS Catalog release %s is already here (%s).", release, d)
    else:
        d.mkdir(parents=True, exist_ok=True)
        for name in FILES.values():
            log.info("Downloading %s …", name)
            _download(config.PGS_META_URL + name, d / name)
        log.info("Downloading the trait categories …")
        (d / CATEGORIES_FILE).write_text(json.dumps(_categories()) + "\n")
        _write_manifest(d, release, latest.get("scores"))
        log.info("PGS Catalog release %s: %s scores, saved to %s.", release,
                 latest.get("scores", "?"), d)
    if dest_root is None:
        _set_current(release)
    return d


def _set_current(release: str) -> None:
    config.PGS_CATALOG_DIR.mkdir(parents=True, exist_ok=True)
    (config.PGS_CATALOG_DIR / "current").write_text(release + "\n")


def install(src: Path) -> Path:
    """Copy a carried snapshot directory in, after checking every file against its manifest.
    It becomes current unless this machine already has a newer release."""
    src = Path(src)
    problems = verify(src)
    if problems:
        raise SystemExit(f"Catalog snapshot {src} is not intact: " + "; ".join(problems))
    release = json.loads((src / MANIFEST).read_text())["release"]
    dest = config.PGS_CATALOG_DIR / release
    if dest.resolve() != src.resolve():
        tmp = config.PGS_CATALOG_DIR / f".{release}.{os.getpid()}.part"
        shutil.rmtree(tmp, ignore_errors=True)
        shutil.copytree(src, tmp)
        shutil.rmtree(dest, ignore_errors=True)
        os.replace(tmp, dest)
    cur = current_dir()
    if cur is None or cur.name <= release:
        _set_current(release)
    log.info("PGS Catalog release %s installed.", release)
    return dest


def current_dir() -> Path | None:
    root = config.PGS_CATALOG_DIR
    try:
        name = (root / "current").read_text().strip()
        if name and (root / name / MANIFEST).is_file():
            return root / name
    except OSError:
        pass
    found = sorted(d for d in root.glob("*") if (d / MANIFEST).is_file()) if root.is_dir() else []
    return found[-1] if found else None


# --- reading it ---------------------------------------------------------------------------
def _rows(path: Path) -> list[dict]:
    csv.field_size_limit(10 ** 8)
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _year(date: str | None) -> str:
    return (date or "")[:4]


class Catalog:
    """One snapshot, loaded into memory (a second or two; a few tens of MB)."""

    def __init__(self, d: Path):
        self.dir = Path(d)
        m = json.loads((self.dir / MANIFEST).read_text())
        self.release, self.fetched = m.get("release", ""), m.get("fetched", "")
        self.scores = {r["Polygenic Score (PGS) ID"]: r
                       for r in _rows(self.dir / FILES["scores"])}
        self.trait_label = {r["Ontology Trait ID"]: r["Ontology Trait Label"]
                            for r in _rows(self.dir / FILES["traits"])}
        self.pubs = {r["PGS Publication/Study (PGP) ID"]: r
                     for r in _rows(self.dir / FILES["publications"])}
        self.sets = {r["PGS Sample Set (PSS)"]: r for r in _rows(self.dir / FILES["sample_sets"])}
        self.perf: dict[str, list[dict]] = {}
        for r in _rows(self.dir / FILES["performance"]):
            self.perf.setdefault(r["Evaluated Score"], []).append(r)
        self.by_trait: dict[str, list[str]] = {}
        for pid, r in self.scores.items():
            for efo in self.efos(pid):
                self.by_trait.setdefault(efo, []).append(pid)
                # The metadata labels traits the traits file may lack; keep the score's own.
        for pid, r in self.scores.items():
            for efo, label in zip(self.efos(pid), self._split(r["Mapped Trait(s) (EFO label)"])):
                self.trait_label.setdefault(efo, label)
        try:
            cats = json.loads((self.dir / CATEGORIES_FILE).read_text())
        except (OSError, ValueError):
            cats = []
        self.category_traits = {c["label"]: [t["id"] for t in c["efotraits"]] for c in cats}
        self.trait_categories: dict[str, list[str]] = {}
        for label, efos in self.category_traits.items():
            for efo in efos:
                self.trait_categories.setdefault(efo, []).append(label)

    @staticmethod
    def _split(text: str | None) -> list[str]:
        return [x.strip() for x in (text or "").split("|") if x.strip()]

    def efos(self, pid: str) -> list[str]:
        return self._split(self.scores.get(pid, {}).get("Mapped Trait(s) (EFO ID)"))

    # Traits and categories, for the picker.
    def categories(self) -> list[dict]:
        out = []
        for label, efos in sorted(self.category_traits.items()):
            have = [e for e in efos if e in self.by_trait]
            out.append({"label": label, "traits": len(have),
                        "scores": len({p for e in have for p in self.by_trait[e]})})
        return out

    def traits(self, query: str = "", category: str = "") -> list[dict]:
        q = query.strip().lower()
        pool = (self.category_traits.get(category, []) if category else list(self.by_trait))
        out = []
        for efo in dict.fromkeys(pool):
            if efo not in self.by_trait:
                continue
            label = self.trait_label.get(efo, efo)
            if q and q not in label.lower() and q != efo.lower():
                continue
            out.append({"efo": efo, "label": label, "scores": len(self.by_trait[efo]),
                        "categories": self.trait_categories.get(efo, [])})
        return sorted(out, key=lambda t: (not label_starts(t["label"], q), t["label"].lower()))

    def scores_for(self, efo: str) -> list[dict]:
        return [self.summary(p) for p in sorted(self.by_trait.get(efo, []))]

    def evaluations(self, pid: str) -> list[dict]:
        out = []
        for r in self.perf.get(pid, []):
            s = self.sets.get(r["PGS Sample Set (PSS)"], {})
            pub = self.pubs.get(r["PGS Publication (PGP) ID"], {})
            n, cases = _num(s.get("Number of Individuals")), _num(s.get("Number of Cases"))
            out.append({
                "ancestry": s.get("Broad Ancestry Category", ""),
                "n": int(n) if n else None, "cases": int(cases) if cases else None,
                "cohort": s.get("Cohort(s)", ""),
                "covariates": r["Covariates Included in the Model"],
                "covariate_class": covariate_class(r["Covariates Included in the Model"]),
                "auroc": _num(r["Area Under the Receiver-Operating Characteristic Curve (AUROC)"]),
                "c_index": _num(r["Concordance Statistic (C-index)"]),
                "or": _num(r["Odds Ratio (OR)"]), "hr": _num(r["Hazard Ratio (HR)"]),
                "other": r["Other Metric(s)"],
                "publication": f"{pub.get('First Author', '')} {_year(pub.get('Publication Date'))}".strip(),
            })
        return out

    @staticmethod
    def _plausible(e: dict) -> bool:
        if e["ancestry"] != ANCESTRY or (e["n"] or 0) < MIN_EVAL_N:
            return False
        cases = e["cases"]
        return cases is None or MIN_CASES <= cases < e["n"]

    def performance(self, pid: str) -> dict:
        """What the published evaluations say the score does, conservatively summarised.

        {kind, median, low, high, cohorts, basis}. basis is `covariate-free` when any such
        evaluation exists; otherwise `covariates not stated` (a `0`/`NA` field), and last
        `with covariates`, which bounds the score from above rather than measuring it.
        Empty when no usable European evaluation exists."""
        evals = [e for e in self.evaluations(pid) if self._plausible(e)]
        for basis, keep in (("covariate-free", lambda e: e["covariate_class"] in ("none", "pcs")),
                            ("covariates not stated", lambda e: e["covariate_class"] == "unclear"),
                            ("with covariates", lambda e: e["covariate_class"] == "adjusted")):
            pool = [e for e in evals if keep(e)]
            for kind, key in (("AUROC", "auroc"), ("C-index", "c_index")):
                vals = [e[key] for e in pool if e[key] and 0.5 < e[key] < 1]
                if vals:
                    return {"kind": kind, "median": round(statistics.median(vals), 3),
                            "low": min(vals), "high": max(vals), "cohorts": len(vals),
                            "basis": basis}
        return {}

    def summary(self, pid: str) -> dict:
        r = self.scores.get(pid)
        if r is None:
            return {}
        pub = self.pubs.get(r["PGS Publication (PGP) ID"], {})
        nvar = int(_num(r["Number of Variants"]) or 0)
        efos = self.efos(pid)
        return {
            "pgs_id": pid, "name": r["PGS Name"], "reported_trait": r["Reported Trait"],
            "efo_ids": efos, "traits": [self.trait_label.get(e, e) for e in efos],
            "n_variants": nvar,
            "engine": "plink2" if nvar > config.PRS_PLINK_MIN else "exact",
            "method": r["PGS Development Method"], "weight_type": r["Type of Variant Weight"],
            "interaction_terms": int(_num(r["Number of Interaction Terms"]) or 0),
            "licence": licence_tag(r["License/Terms of Use"]),
            "released": r["Release Date"],
            "publication": {"first_author": pub.get("First Author", ""),
                            "year": _year(pub.get("Publication Date")),
                            "journal": pub.get("Journal Name", ""),
                            "doi": pub.get("digital object identifier (doi)", "")},
            "gwas_ancestry": r["Ancestry Distribution (%) - Source of Variant Associations (GWAS)"],
            "training_ancestry": r["Ancestry Distribution (%) - Score Development/Training"],
            "performance": self.performance(pid),
            "unsupported": unsupported(r["Type of Variant Weight"],
                                       int(_num(r["Number of Interaction Terms"]) or 0)),
        }

    def newer(self, pid: str) -> list[str]:
        """Scores for any of this score's traits released after it, newest first."""
        mine = self.scores.get(pid, {}).get("Release Date", "")
        out = {p for efo in self.efos(pid) for p in self.by_trait.get(efo, [])
               if p != pid and self.scores[p]["Release Date"] > mine}
        return sorted(out, key=lambda p: self.scores[p]["Release Date"], reverse=True)


_LINEAR_WEIGHT = re.compile(r"odds ratio|hazard ratio|^\s*(or|hr)\s*$", re.I)


def unsupported(weight_type: str | None, interaction_terms: int = 0) -> str | None:
    """Why the additive engine cannot score this, from catalog metadata, or None.

    The engine sums dosage × weight. That is the model only for weights on an additive scale
    (beta, log odds, log hazard): a weight published as a plain odds or hazard ratio would
    have to be logged first, and interaction terms are not modelled at all."""
    wt = (weight_type or "").strip()
    if wt and _LINEAR_WEIGHT.search(wt) and "log" not in wt.lower() and "ln(" not in wt.lower():
        return f"weights are published as '{wt}', not on an additive (log) scale"
    if interaction_terms:
        return f"{interaction_terms} interaction term(s), which this engine does not model"
    return None


DISEASE_CATEGORIES = ("Cancer", "Cardiovascular disease", "Digestive system disorder",
                      "Immune system disorder", "Metabolic disorder", "Neurological disorder",
                      "Other disease")
EVIDENCE = ("none", "european", "covariate-free")


def select(cat: "Catalog", categories: list[str] | None = None, evidence: str = "none",
           engine: str = "all", include_unsupported: bool = False) -> list[str]:
    """Score ids by trait category, published evidence and size, for `pgs fetch --from-file`.

    evidence: `european` keeps scores with a usable European evaluation (the performance
    summary is not empty); `covariate-free` keeps those whose summary is covariate-free.
    Scores the additive engine cannot model are left out unless asked for."""
    if evidence not in EVIDENCE:
        raise SystemExit(f"--evidence is one of {', '.join(EVIDENCE)}")
    if categories:
        unknown = [c for c in categories if c not in cat.category_traits]
        if unknown:
            raise SystemExit(f"unknown categor{'y' if len(unknown) == 1 else 'ies'}: "
                             f"{', '.join(unknown)} (known: {', '.join(sorted(cat.category_traits))})")
        pids = {p for c in categories for e in cat.category_traits[c]
                for p in cat.by_trait.get(e, [])}
    else:
        pids = set(cat.scores)
    out = []
    for pid in sorted(pids):
        s = cat.summary(pid)
        if s["unsupported"] and not include_unsupported:
            continue
        if engine != "all" and s["engine"] != {"small": "exact", "genome-wide": "plink2"}[engine]:
            continue
        perf = s["performance"]
        if evidence == "european" and not perf:
            continue
        if evidence == "covariate-free" and perf.get("basis") != "covariate-free":
            continue
        out.append(pid)
    return out


def label_starts(label: str, q: str) -> bool:
    return bool(q) and label.lower().startswith(q)


_loaded: dict[str, Catalog] = {}


def load(d: Path | None = None) -> Catalog:
    """The current snapshot (or `d`), cached per process. SystemExit if there is none."""
    d = Path(d) if d else current_dir()
    if d is None:
        raise SystemExit("No PGS Catalog snapshot here. Run `run.py pgs sync` on a connected "
                         "machine (or `pgs import` a carried one).")
    key = str(d.resolve())
    if key not in _loaded:
        _loaded.clear()
        _loaded[key] = Catalog(d)
    return _loaded[key]


def describe_performance(p: dict) -> str:
    """One line for a report or a table cell."""
    if not p:
        return "no usable European evaluation"
    span = (f"{p['median']:.2f}" if p["cohorts"] == 1
            else f"{p['median']:.2f} median of {p['cohorts']} cohorts ({p['low']:.2f}–{p['high']:.2f})")
    note = {"covariate-free": "", "covariates not stated": ", covariates not stated",
            "with covariates": ", only with age/sex or other covariates"}[p["basis"]]
    return f"{p['kind']} {span}{note}"
