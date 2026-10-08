#!/usr/bin/env python3
"""Generate panels/complex_traits.tsv from PGS Catalog metadata, by objective criteria.

Takes no hand-written trait list as input — the panel falls out of the criteria in
docs/polygenic-panel-design.md §3 applied to PGS Catalog metadata, so it can be regenerated
by anyone from the Catalog alone. Keep it that way.

    python3 panels/build_complex_traits.py [--out panels/complex_traits.tsv]

Criteria: PGS-alone performance metric (no covariates) + ancestry-matched evaluation cohort
>= MIN_EVAL_N. Curation lists below remove ontology artefacts and redundant subtypes the
criteria alone let through, and assign actionability / S1 gating. Score licences are recorded
in a column, not filtered on — see _LICENCE_TAGS.
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import sys
import urllib.request

FTP = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/metadata/"
MIN_EVAL_N = 1000
ANCESTRY = "European"

# Licences are *recorded, not filtered on*. Personal/research use — this project's stated
# posture — satisfies every term in the Catalog: NonCommercial is met by definition, and
# NoDerivatives restricts redistribution, which never happens (prs.py fetches at runtime).
# Filtering here would silently drop good scores for a use their licence explicitly allows.
# The label exists so anyone redeploying commercially can see what needs clearing first.
_LICENCE_TAGS = (
    ("should be cited appropriately", "pgs-catalog-default"),
    ("CC0", "CC0"),
    ("NonCommercial-NoDerivatives", "CC-BY-NC-ND"),
    ("Creative Commons Attribution 4.0", "CC-BY-4.0"),
    ("academic community", "academic-noncommercial"),
    ("academic research purposes only", "academic-noncommercial"),
    ("research purposes only", "research-only"),
)


def _licence_tag(text: str) -> str:
    for needle, tag in _LICENCE_TAGS:
        if needle.lower() in (text or "").lower():
            return tag
    return "other"

# Ontology artefacts, proxy phenotypes and redundant subtypes — see design doc §3.
DROP = {
    "female", "male", "comparative body size at age 10, self-reported", "Back pain",
    "Cleft palate", "cleft lip", "luminal A breast carcinoma", "luminal B breast carcinoma",
    "triple-negative breast carcinoma", "Her2-receptor negative breast cancer",
    "HER2 positive breast carcinoma", "estrogen-receptor positive breast cancer",
    "ovarian serous carcinoma", "childhood onset asthma", "Asthma", "colorectal carcinoma",
    "rheumatoid factor-negative juvenile idiopathic arthritis",
    "oligoarticular juvenile idiopathic arthritis",
    "enthesitis-related juvenile idiopathic arthritis", "melanoma",
    "squamous cell carcinoma", "hemorrhoid", "refractive error measurement",
    "alcoholic liver cirrhosis", "REM sleep behavior disorder", "Brugada syndrome",
}

# screen = surveillance pathway exists; treat = preventive therapy; life = lifestyle; none.
# First-pass judgement — needs review against actual guidelines before it drives a report.
ACTIONABILITY = {
    "type 1 diabetes mellitus": "screen", "celiac disease": "screen",
    "inflammatory bowel disease": "screen",
    "systemic lupus erythematosus": "screen", "coronary artery disorder": "treat",
    "prostate carcinoma": "screen", "brain aneurysm": "screen", "keratoconus": "screen",
    "type 2 diabetes mellitus": "life", "body mass index": "life",
    "multiple sclerosis": "none", "gout": "treat", "hypertensive disorder": "treat",
    "glaucoma": "screen", "atrial fibrillation": "screen",
    "abdominal aortic aneurysm": "screen", "hypothyroidism": "screen",
    "systolic blood pressure": "treat", "triglyceride measurement": "treat",
    "breast carcinoma": "screen", "urate measurement": "treat",
    "venous thromboembolism": "treat",
    "low density lipoprotein cholesterol measurement": "treat",
    "high density lipoprotein cholesterol measurement": "life", "asthma": "treat",
    "juvenile idiopathic arthritis": "screen", "cutaneous melanoma": "screen",
    "angina pectoris": "treat", "vitamin D level": "treat", "ovarian carcinoma": "screen",
    "basal cell carcinoma": "screen", "heart failure": "treat",
    "peripheral arterial disease": "treat", "lung carcinoma": "screen",
    "Ischemic stroke": "treat", "endometriosis": "screen", "stroke disorder": "treat",
    "chronic kidney disease": "screen", "colorectal cancer": "screen",
    "ulcerative colitis": "screen", "renal carcinoma": "none",
    "endometrial carcinoma": "screen", "urinary bladder carcinoma": "none",
    "exocrine pancreatic carcinoma": "none",
}

# Severe + non-actionable: kept in the panel, deferred to the S1 consent filter.
S1_GATE = {"multiple sclerosis", "renal carcinoma", "exocrine pancreatic carcinoma",
           "urinary bladder carcinoma", "lung carcinoma"}


def _fetch(name: str) -> list[dict]:
    with urllib.request.urlopen(FTP + name, timeout=180) as fh:
        text = fh.read().decode("utf-8", "replace")
    return list(csv.DictReader(io.StringIO(text)))


def _num(x: str | None) -> float | None:
    m = re.match(r"\s*(-?\d+\.?\d*)", x or "")
    return float(m.group(1)) if m else None


def build() -> list[dict]:
    csv.field_size_limit(10 ** 7)
    scores = {r["Polygenic Score (PGS) ID"]: r for r in _fetch("pgs_all_metadata_scores.csv")}
    perf = _fetch("pgs_all_metadata_performance_metrics.csv")
    sets = {r["PGS Sample Set (PSS)"]: r
            for r in _fetch("pgs_all_metadata_evaluation_sample_sets.csv")}
    labels = {t["Ontology Trait ID"]: t["Ontology Trait Label"]
              for t in _fetch("pgs_all_metadata_efo_traits.csv")}

    best: dict[str, tuple] = {}
    for p in perf:
        score = scores.get(p["Evaluated Score"])
        if not score:
            continue
        if (p["Covariates Included in the Model"] or "").strip():
            continue                                   # PGS-alone metrics only
        sample = sets.get(p["PGS Sample Set (PSS)"])
        if not sample or sample["Broad Ancestry Category"] != ANCESTRY:
            continue
        n = _num(sample["Number of Individuals"]) or 0
        if n < MIN_EVAL_N:
            continue
        auroc = _num(p["Area Under the Receiver-Operating Characteristic Curve (AUROC)"])
        metric = auroc or _num(p["Concordance Statistic (C-index)"])
        if not metric or not 0.5 < metric < 1:
            continue
        rec = (metric, score["Polygenic Score (PGS) ID"],
               int(_num(score["Number of Variants"]) or 0), int(n),
               "AUROC" if auroc else "C-index",
               _licence_tag(score["License/Terms of Use"]))
        for efo in (e.strip() for e in re.split(r"[|,]", score["Mapped Trait(s) (EFO ID)"] or "")):
            if efo and (efo not in best or metric > best[efo][0]):
                best[efo] = rec

    rows = []
    for efo, (metric, pgs, nvar, n, kind, licence) in best.items():
        trait = labels.get(efo, "")
        if not trait or trait in DROP:
            continue
        rows.append({"efo_id": efo, "trait": trait, "pgs_id": pgs, "n_variants": nvar,
                     "metric": round(metric, 3), "metric_kind": kind, "eval_n": n,
                     "licence": licence,
                     "actionability": ACTIONABILITY.get(trait, "review"),
                     "s1_gate": "yes" if trait in S1_GATE else "no"})
    rows.sort(key=lambda r: -r["metric"])
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="panels/complex_traits.tsv")
    args = ap.parse_args()
    rows = build()
    if not rows:
        print("no traits met the criteria — refusing to write an empty panel", file=sys.stderr)
        return 1
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    unreviewed = sum(1 for r in rows if r["actionability"] == "review")
    print(f"wrote {args.out}: {len(rows)} traits"
          + (f" ({unreviewed} need actionability review)" if unreviewed else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
