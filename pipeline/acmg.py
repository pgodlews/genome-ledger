"""Stage: ACMG/AMP-style variant classification.

Moves from "ClinVar says X" to a *reasoned* call with evidence codes. For each variant in
the monogenic-panel genes it evaluates a transparent subset of the ACMG/AMP criteria from
annotations the pipeline already has in DuckDB — VEP consequence, gnomAD allele frequency,
SIFT/PolyPhen, AlphaMissense (missense), and SpliceAI (splice) — then combines them into a
5-tier classification.

ClinVar is deliberately **held out** of the criteria so the call is independent of it; the
report then compares the two, and its highest-value output is **divergences**: variants this
classifier rates (Likely) Pathogenic that ClinVar leaves VUS/absent (candidate
reclassifications to review), or vice versa.

This is a research-grade *subset*, not a clinical InterVar run. Implemented criteria:
  PVS1  null variant (stop/frameshift/canonical-splice/start-loss) in a panel disease gene,
        graded by gnomAD LOF constraint: full PVS1 in a haploinsufficiency-plausible gene
        (recessive genes, or dominant/XL genes that are LOF-constrained); downgraded to
        PVS1_moderate in a dominant/XL gene that gnomAD shows is LOF-tolerant
  PM2   absent or ultra-rare in gnomAD (AF < 1e-4)            [treated as supporting]
  PP3   in-silico pathogenic (SpliceAI splice-altering, AlphaMissense, or SIFT+PolyPhen
        agree damaging) — also grades near-splice (splice_region) variants via SpliceAI
  BA1   gnomAD AF > 5%  (stand-alone benign)
  BS1   gnomAD AF > 1%  (strong benign)
  BP4   in-silico benign (SpliceAI predicts no splice effect, AlphaMissense, or SIFT+PolyPhen
        agree tolerated)
Not implemented (need data we don't have): PS1/PS3/PM1/PM5/PP1/segregation/functional.
"""

from __future__ import annotations

from pathlib import Path

from . import reportpaths
from . import allele_balance, alphamissense, config, gnomad_constraint, limitations, spliceai
from .load import TABLE, connect
from .panel import load_panel, sql_in_list
from .snapshot import load_manifest
from .util import log, write_report

# Allele-frequency thresholds (generic; BS1 is really disorder-specific).
_BA1_AF = 0.05
_BS1_AF = 0.01
_PM2_AF = 1e-4

# Indels at/above this size are SV-scale: their SNV-style consequence annotation is
# unreliable (often repeat-region artifacts), and real large del/dup is the `sv` stage's
# job (C1). Deferred here rather than confidently called PVS1.
_MAX_INDEL_BP = 50

# VEP consequence terms, most-severe first (consequence is comma/&-joined per variant).
_NULL = ("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
         "stop_gained", "frameshift_variant", "start_lost")
_SEVERITY = _NULL + ("stop_lost", "inframe_deletion", "inframe_insertion",
                     "missense_variant", "protein_altering_variant",
                     "splice_region_variant")
_MISSENSE = ("missense_variant", "protein_altering_variant")
# Consequences worth classifying (skip synonymous/intronic/UTR/intergenic noise).
_CODING = set(_SEVERITY)
# Consequences whose pathogenicity an in-silico predictor can speak to: missense (→ AlphaMissense)
# plus near-splice changes (→ SpliceAI). The canonical splice donor/acceptor nulls already earn
# PVS1, so they're not re-graded here.
_INSILICO = set(_MISSENSE) | {"splice_region_variant"}

_TIER_RANK = {"Pathogenic": 0, "Likely pathogenic": 1, "VUS": 2,
              "Likely benign": 3, "Benign": 4}


def _primary_csq(consequence: str) -> str:
    terms = (consequence or "").replace("&", ",").split(",")
    for sev in _SEVERITY:
        if sev in terms:
            return sev
    return terms[0] if terms else ""


def _path_insilico(sift: str, polyphen: str, am, sa) -> bool:
    if spliceai.is_pathogenic(sa):       # predicted splice-altering (≥ recommended cut-off)
        return True
    if am and am[1] == "likely_pathogenic":
        return True
    return (sift or "").startswith("deleterious") and \
           (polyphen or "").startswith("probably_damaging")


def _benign_insilico(sift: str, polyphen: str, am, sa) -> bool:
    if spliceai.is_pathogenic(sa):       # a real splice signal vetoes a benign call
        return False
    if am and am[1] == "likely_benign":
        return True
    # A near-splice variant with no missense data: SpliceAI predicting no effect is the
    # only in-silico evidence available, so let it stand alone for BP4.
    if spliceai.is_benign(sa) and am is None and not sift and not polyphen:
        return True
    return (sift or "").startswith("tolerated") and (polyphen or "").startswith("benign")


def _pvs1_code(dominant: bool, constrained: bool | None) -> str:
    """Grade PVS1 by loss-of-function-mechanism plausibility (gnomAD constraint).

    A recessive-only gene's LOF is pathogenic without needing haploinsufficiency, so it keeps
    full PVS1. In a dominant/XL gene, full PVS1 holds when the gene is LOF-constrained or we
    have no constraint data (don't penalize the unknown); it drops to PVS1_moderate only when
    gnomAD shows the gene tolerates LOF (haploinsufficiency unlikely → a null is weaker)."""
    if dominant and constrained is False:
        return "PVS1_moderate"
    return "PVS1"


def _codes(csq: str, af, sift, polyphen, am, sa, gene_is_disease: bool,
           dominant: bool = False, constrained: bool | None = None) -> list[str]:
    codes: list[str] = []
    if af is not None and af > _BA1_AF:
        codes.append("BA1")
    elif af is not None and af > _BS1_AF:
        codes.append("BS1")
    if af is None or af < _PM2_AF:
        codes.append("PM2")
    if csq in _NULL and gene_is_disease:
        codes.append(_pvs1_code(dominant, constrained))
    if csq in _INSILICO:
        if _path_insilico(sift, polyphen, am, sa):
            codes.append("PP3")
        elif _benign_insilico(sift, polyphen, am, sa):
            codes.append("BP4")
    return codes


def _combine(codes: list[str]) -> str:
    """Simplified ACMG combine over the implemented codes. PM2 is treated as supporting;
    PVS1_moderate is a constraint-downgraded null (moderate, not very-strong).

    NOTE — "Pathogenic" is unreachable with the codes implemented today, and that is a
    property of the code set, not of the family's variants. The tier needs PVS1 plus two
    supporting codes (Richards 2015 Table 5, Pathogenic (i)(d)), and the only two supporting
    codes here are PM2 and PP3: PVS1 requires a consequence in _NULL, PP3 requires one in
    _INSILICO, and those sets are disjoint, so path_sup can never exceed 1 alongside PVS1.
    The ceiling is therefore Likely pathogenic. The branch is kept (it becomes reachable the
    moment a second supporting code such as PP1/PM5 is implemented) and the report's caveats
    state the ceiling, so an empty Pathogenic count is never read as a finding.
    """
    s = set(codes)
    if "BA1" in s:
        return "Benign"
    path_strong = "PVS1" in s
    path_mod = "PVS1_moderate" in s
    path_sup = sum(c in s for c in ("PM2", "PP3"))
    ben = ("BS1" in s) + ("BP4" in s)
    if (path_strong or path_mod) and ben:         # null variant but a benign signal too
        return "VUS"
    if path_strong:
        return "Pathogenic" if path_sup >= 2 else "Likely pathogenic"
    if path_mod:                                  # moderate null: needs a supporting code
        return "Likely pathogenic" if path_sup >= 1 else "VUS"
    if path_sup and ben:
        return "VUS"
    if ben >= 2:
        return "Likely benign"
    return "VUS"                                   # lone supporting evidence stays uncertain


def _clinvar_class(sig: str | None) -> str:
    s = (sig or "").lower()
    # Checked BEFORE the pathogenic/benign tests: "Conflicting_classifications_of_pathogenicity"
    # contains "pathogenic" (and the older "Conflicting_interpretations_of_pathogenicity" does
    # too), so it used to fall past every branch into "other" and be invisible to _divergence
    # — silently excluding the single case reclassification review exists for: a submitter
    # disagreement that our own reasoning can break the tie on.
    if "conflict" in s:
        return "conflicting"
    if "pathogenic" in s:
        return "path"
    if "benign" in s:
        return "benign"
    if "uncertain" in s or "vus" in s:
        return "vus"
    return "absent" if not s or s == "-" else "other"


def _divergence(tier: str, sig: str | None) -> str:
    """Flag where the reasoned call materially disagrees with ClinVar."""
    cv = _clinvar_class(sig)
    path_tier = tier in ("Pathogenic", "Likely pathogenic")
    ben_tier = tier in ("Benign", "Likely benign")
    # A reasoned P/LP against a ClinVar *benign* call is the sharpest disagreement there is,
    # and was the one direction never flagged — only benign-vs-pathogenic was.
    if path_tier and cv == "benign":
        return "⚠️ contradicts ClinVar (benign) — review"
    if ben_tier and cv == "path":
        return "⚠️ weaker than ClinVar — review"
    if path_tier and cv == "conflicting":
        return "⚠️ resolves a conflicting ClinVar entry toward pathogenic — review"
    if ben_tier and cv == "conflicting":
        return "⚠️ resolves a conflicting ClinVar entry toward benign — review"
    if path_tier and cv in ("vus", "absent"):
        return "⚠️ stronger than ClinVar — review"
    return ""


def panel_context() -> tuple[set[str], dict[str, set[str]]]:
    """Panel disease genes + per-gene inheritance modes — the gene context both the ACMG
    report and the triage layer need. A null in a gene with any dominant/XL association is
    where haploinsufficiency (hence gnomAD constraint) is decisive for PVS1 strength."""
    panel = load_panel()
    disease_genes = {g for d in panel for g in d.genes}
    gene_modes: dict[str, set[str]] = {}
    for d in panel:
        for g in d.genes:
            gene_modes.setdefault(g, set()).add(d.inh)
    return disease_genes, gene_modes


def evaluate(gene, chrom, pos, ref, alt, csq_raw, sig, af, sift, pp, zyg, *,
             disease_genes: set[str], gene_modes: dict[str, set[str]],
             have_am: bool, have_sa: bool, have_gc: bool) -> dict | None:
    """Classify one variant with the implemented ACMG/AMP subset. Returns the call dict, or
    None for consequences not worth classifying (synonymous/intronic/UTR/intergenic). This is
    the single source of truth shared by `acmg_report` and the triage layer (`triage.py`).

    The returned dict carries `deferred=True` for SV-scale indels (consequence annotation
    unreliable — the `sv` stage's job); those have `tier=None`."""
    csq = _primary_csq(csq_raw)
    if csq not in _CODING:
        return None
    base = {"gene": gene, "label": f"{chrom}:{pos} {ref}>{alt}", "csq": csq,
            "af": af, "sig": sig, "zyg": (zyg or "").title(),
            "constraint": gnomad_constraint.label(gene) if have_gc else "—"}
    ref_len = 0 if ref == "-" else len(ref)
    alt_len = 0 if alt == "-" else len(alt)
    if abs(alt_len - ref_len) > _MAX_INDEL_BP:   # SV-scale → defer to the `sv` stage
        return {**base, "tier": None, "codes": [], "div": "", "deferred": True}
    am = (alphamissense.lookup(chrom, pos, ref, alt)
          if have_am and csq in _MISSENSE else None)
    sa = (spliceai.lookup(chrom, pos, ref, alt)
          if have_sa and csq in _INSILICO else None)
    dominant = bool(gene_modes.get(gene, set()) & {"AD", "XL"})
    constrained = gnomad_constraint.is_constrained(gene) if have_gc else None
    codes = _codes(csq, af, sift, pp, am, sa, gene in disease_genes, dominant, constrained)
    tier = _combine(codes)
    return {**base, "tier": tier, "codes": codes, "div": _divergence(tier, sig),
            "deferred": False}


def acmg_report(snapshot_id: str) -> list[Path]:
    disease_genes, gene_modes = panel_context()
    gene_list = sql_in_list(disease_genes)
    con = connect()
    meta = load_manifest(snapshot_id)
    have_am = alphamissense.available()
    have_sa = spliceai.available()
    have_gc = gnomad_constraint.available()
    written = []

    samples = [r[0] for r in con.execute(
        f"SELECT DISTINCT sample_id FROM {TABLE} WHERE snapshot_id = ? ORDER BY 1",
        [snapshot_id]).fetchall()]

    for sample in samples:
        rows = con.execute(f"""
            SELECT gene, chrom, pos, ref, alt, consequence, clinvar_sig,
                   gnomad_af, sift, polyphen, zygosity
            FROM {TABLE}
            WHERE snapshot_id = ? AND sample_id = ? AND gene IN ({gene_list})
        """, [snapshot_id, sample]).fetchall()

        calls = []
        large_indel = 0
        ab = allele_balance.fractions(sample, [r[1:5] for r in rows])
        for gene, chrom, pos, ref, alt, csq_raw, sig, af, sift, pp, zyg in rows:
            call = evaluate(gene, chrom, pos, ref, alt, csq_raw, sig, af, sift, pp, zyg,
                            disease_genes=disease_genes, gene_modes=gene_modes,
                            have_am=have_am, have_sa=have_sa, have_gc=have_gc)
            if call is None:
                continue
            # Display only: the tier stays evidence-from-annotation; a low-fraction het is
            # flagged for review in the zygosity cell, never reclassified here.
            call["zyg"] = allele_balance.label(call["zyg"], ab.get((chrom, pos, ref, alt)))
            if call["deferred"]:   # SV-scale indel → counted, deferred to the `sv` stage
                large_indel += 1
                continue
            calls.append(call)

        calls.sort(key=lambda c: (_TIER_RANK.get(c["tier"], 9), c["gene"]))
        diverging = [c for c in calls if c["div"]]
        by_tier = {t: sum(c["tier"] == t for c in calls) for t in _TIER_RANK}

        def _af(x):
            return "absent" if x is None else (f"{x:.2g}")

        def _row(c):
            return (f"| {c['gene']} | {c['label']} ({c['zyg']}) | {c['csq']} | "
                    f"**{c['tier']}** | {', '.join(c['codes']) or '—'} | {_af(c['af'])} | "
                    f"{c['constraint']} | {c['sig'] or '—'} | {c['div']} |")

        lines = [
            f"# ACMG-style classification — {sample}",
            "",
            f"- **Snapshot:** `{snapshot_id}` (ClinVar {meta['clinvar_date']})",
            f"- **Scope:** coding/splice variants in the {len(disease_genes)} panel genes "
            f"(`{config.PANEL_FILE.name}`); AlphaMissense {'on' if have_am else 'absent'}, "
            f"SpliceAI {'on' if have_sa else 'absent'}, "
            f"gnomAD {config.GNOMAD_CONSTRAINT_VERSION} constraint "
            f"{'on' if have_gc else 'absent'}",
            f"- **Classified:** {len(calls)} — "
            + ", ".join(f"{by_tier[t]} {t}" for t in _TIER_RANK if by_tier[t]),
            f"- **Diverge from ClinVar (review):** {len(diverging)}"
            + (f"  ·  {large_indel} SV-scale indel(s) (>{_MAX_INDEL_BP} bp) deferred to the "
               "`sv` stage" if large_indel else ""),
            "",
            "## Candidate reclassifications (reasoned call ≠ ClinVar)",
            "",
            "| Gene | Variant | Consequence | ACMG | Codes | gnomAD | Constraint | ClinVar | Note |",
            "|------|---------|-------------|------|-------|--------|-----------|---------|------|",
            *([_row(c) for c in diverging] or ["| _none_ | | | | | | | | |"]),
            "",
            "## All classified variants",
            "",
            "| Gene | Variant | Consequence | ACMG | Codes | gnomAD | Constraint | ClinVar | Note |",
            "|------|---------|-------------|------|-------|--------|-----------|---------|------|",
            *([_row(c) for c in calls] or ["| _none_ | | | | | | | | |"]),
            "",
            "## Method & caveats",
            "- **Research-grade subset of ACMG/AMP** (Richards 2015), not a clinical "
            "InterVar run. Implemented codes: PVS1 (constraint-graded), PM2 (supporting), "
            "PP3, BA1, BS1, BP4.",
            "- **The top tier is out of reach here: `Likely pathogenic` is the ceiling.** "
            "`Pathogenic` needs PVS1 plus two supporting codes, and of the two implemented "
            "(PM2, PP3) only PM2 can ever accompany PVS1 — PP3 applies to missense/near-splice "
            "consequences, PVS1 to nulls. So an absence of Pathogenic calls says nothing about "
            "the variants; read `Likely pathogenic` as this report's strongest call.",
            "- **ClinVar is held out** of the criteria, so the call is independent of it and "
            "the comparison is meaningful. A ⚠️ flag means review, not a conclusion.",
            f"- **PVS1 is graded by gnomAD {config.GNOMAD_CONSTRAINT_VERSION} constraint** "
            f"(pLI ≥ {config.GNOMAD_PLI_CONSTRAINED} or LOEUF < {config.GNOMAD_LOEUF_CONSTRAINED} "
            "= LOF-constrained): a null keeps full PVS1 in recessive genes and in "
            "LOF-constrained dominant/XL genes, but drops to **PVS1_moderate** in a dominant/XL "
            "gene gnomAD shows tolerates LOF (so a lone null there lands at VUS, not Likely "
            "pathogenic). Still a constraint proxy, not the full ClinGen PVS1 decision tree.",
            "- **Constraint under-flags two-hit tumor-suppressors** (e.g. BRCA1/2 score as "
            "LOF-tolerant because heterozygous carriers are viable), so a truncating variant "
            "there may grade PVS1_moderate; it still reaches Likely pathogenic when rare (PM2). "
            "Read the constraint column alongside the call.",
            "- Not implemented (data unavailable): PS1/PS3/PM1/PM5/PP1, segregation, "
            "functional, phasing. Most novel variants therefore land at VUS — by design.",
            "- Research-grade, not a clinical diagnosis.",
            "",
        ]
        lines += limitations.render(snapshot_id, report_kind="acmg", sample=sample)
        path = reportpaths.snapshot_report("acmg", sample, snapshot_id)
        write_report(path, lines)
        written.append(path)
        log.info("ACMG report: %s (%d classified, %d diverge from ClinVar)",
                 path, len(calls), len(diverging))
    con.close()
    return written
