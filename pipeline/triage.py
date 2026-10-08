"""Delta triage: multi-factor classification of every changed variant.

Replaces the old fixed 0–4 priority ladder (which keyed almost entirely off ClinVar text
transitions) with a weighted, auditable score over five factors — reasoned ACMG evidence,
in-silico model distance-from-threshold, population rarity, cross-snapshot stability, and
gene context — mapped to four tiers: actionable_now / monitor / log_only / ignore.

`score_delta` is a *pure* function: the orchestrator (`diff.py`) and the retrospective
validation harness (`validate.py`) both call it, so a tier is always reproducible from its
stored `factors` JSON. All weights and cut-offs live in `config.TRIAGE_*`.
"""

from __future__ import annotations

from . import alphamissense, clinsig, config, enformer, spliceai

TIERS = ("actionable_now", "monitor", "log_only", "ignore")

# gnomAD AF boundaries for the rarity factor (mirror the ACMG PM2/BS1 intent).
_AF_ULTRA_RARE = 1e-4   # PM2-scale: absent/ultra-rare → maximal interest
_AF_COMMON = 0.01       # BS1-scale: common → push toward ignore


def is_actionable(sig: str | None) -> bool:
    return clinsig.is_actionable(sig)


def is_uncertain(sig: str | None) -> bool:
    return clinsig.is_uncertain(sig)


def _review_stars(revstat: str | None) -> int:
    return config.REVIEW_STATUS_RANK.get((revstat or "").replace(" ", "_"), 0)


def change_kind(prev_sig, curr_sig, prev_rev, curr_rev, prev_af, curr_af):
    """Detect *what changed* between snapshots. Returns (kind, description) or None when the
    variant's annotation is unchanged in any way triage cares about. Generalizes the old
    diff._classify ladder — the ClinVar transition still defines the change, the score (below)
    decides how loudly to surface it."""
    if is_actionable(curr_sig) and not is_actionable(prev_sig):
        if is_uncertain(prev_sig) or prev_sig is None:
            return ("new_pathogenic", f"'{prev_sig or 'absent'}' → '{curr_sig}'")
        return ("upgraded", f"'{prev_sig}' → '{curr_sig}'")
    if is_actionable(prev_sig) and not is_actionable(curr_sig):
        return ("downgraded", f"'{prev_sig}' → '{curr_sig}'")
    if (prev_sig or "") != (curr_sig or ""):
        return ("reclassified", f"'{prev_sig or 'absent'}' → '{curr_sig or 'absent'}'")
    if is_actionable(curr_sig) and _review_stars(curr_rev) > _review_stars(prev_rev):
        return ("reclassified", f"confidence ↑: '{prev_rev}' → '{curr_rev}'")
    if prev_af is not None and curr_af is not None:
        was_rare = prev_af < config.RARE_AF_THRESHOLD
        now_rare = curr_af < config.RARE_AF_THRESHOLD
        if was_rare != now_rare:
            return ("af_shift", f"gnomAD AF {prev_af:.2g} → {curr_af:.2g} "
                                f"({'now rare' if now_rare else 'now common'})")
    return None


def _noncoding(consequence: str | None) -> bool:
    terms = set((consequence or "").replace("&", ",").split(","))
    return bool(terms & set(config.INTERPRET_CONSEQUENCES))


def _am_applicable(terms: str) -> bool:
    """Whether AlphaMissense speaks to this consequence. Imported from acmg so the two stages
    cannot drift apart again (they disagreed on protein_altering_variant)."""
    from .acmg import _MISSENSE

    return bool(set(terms.split(",")) & set(_MISSENSE))


def model_lookups(chrom, pos, ref, alt, consequence, have_am, have_sa, have_en):
    """Resolve the in-silico signals relevant to this variant's consequence class — the
    gating that keeps a coding variant from drawing (false) confidence off Enformer, and a
    promoter variant off AlphaMissense. Returns (am, sa, en), any of which may be None."""
    terms = (consequence or "").replace("&", ",")
    # Same consequence set acmg._MISSENSE gates AlphaMissense on. A substring test for
    # "missense" excluded `protein_altering_variant`, so the ACMG stage looked AM up for
    # such a variant and the triage stage did not — one verdict computed from two different
    # evidence sets, for the same variant, in the same report.
    am = (alphamissense.lookup(chrom, pos, ref, alt)
          if have_am and _am_applicable(terms) else None)
    sa = (spliceai.lookup(chrom, pos, ref, alt)
          if have_sa and "splice" in terms else None)
    en = (enformer.lookup(chrom, pos, ref, alt)
          if have_en and _noncoding(consequence) else None)
    return am, sa, en


def _clinvar_evidence(curr_sig, curr_rev) -> float:
    """ClinVar's *own* current call as evidence, scaled by review confidence (stars). An
    expert-reviewed Pathogenic is strong evidence even where our limited reasoner — which only
    covers panel genes and lands most novel variants at VUS — has nothing to say."""
    if not is_actionable(curr_sig):
        return 0.0
    return min(1.0, 0.5 + 0.125 * _review_stars(curr_rev))  # 0★→0.5 … 4★→1.0


def _f_acmg(acmg_call: dict | None, clinvar_ev: float) -> float:
    """Evidence-strength factor → [0, ~1.3]. The stronger of the reasoned ACMG tier and
    ClinVar's own confidence-scaled call; a divergence between the two (candidate
    reclassification) earns a bonus on top."""
    base = clinvar_ev
    if acmg_call and acmg_call.get("tier"):
        base = max(base, config.TRIAGE_ACMG_SCORE.get(acmg_call["tier"], 0.3))
        if acmg_call.get("div"):
            base += config.TRIAGE_DIVERGENCE_BONUS
    return base


def _f_model(am, sa, en, ref, alt) -> float:
    """Strength of in-silico evidence **for pathogenicity**, in [0, 1]. Takes the strongest of
    the applicable predictors (already gated by consequence in model_lookups).

    One-sided on purpose. This factor carries a positive weight into a score whose whole
    question is "how much should this be acted on", so only evidence pointing at pathogenic
    may raise it. The old `abs(am - 0.5) * 2` was symmetric: a confidently BENIGN
    AlphaMissense call (0.02) scored 0.96, exactly as much as a confidently pathogenic one
    (0.98). A rare on-panel missense AM called benign then summed to 0.602 against the 0.60
    `actionable_now` cut-off — promoted to the top tier by benign evidence. SpliceAI's ramp
    below was already one-sided; this now matches it."""
    cands = []
    if am:  # (pathogenicity 0–1, class) → ramp above the 0.5 boundary; benign side → 0
        cands.append(max(0.0, min(1.0, (am[0] - 0.5) * 2)))
    if sa:  # SpliceAI delta 0–1 → ramp from the recommended cut-off (keep in sync w/ config)
        cut = config.SPLICEAI_DS_RECOMMENDED
        cands.append(max(0.0, min(1.0, (sa[0] - cut) / (1.0 - cut))))
    if en:  # Enformer/Borzoi: calibrated band (percentile vs benign null, or raw fallback)
        cands.append({"high": 1.0, "moderate": 0.5, "low": 0.1}[enformer.band(en[0], ref, alt)])
    return max(cands) if cands else 0.0


def _f_freq(curr_af) -> float:
    """Population rarity. Absent/ultra-rare → 1.0; common → negative (pushes toward ignore)."""
    if curr_af is None or curr_af < _AF_ULTRA_RARE:
        return 1.0
    if curr_af < config.RARE_AF_THRESHOLD:
        return 0.6
    if curr_af >= _AF_COMMON:
        return -1.0
    return 0.0


def path_ordinal(sig: str | None) -> int:
    """Coarse pathogenicity level for a ClinVar significance: benign −1, vus/absent 0,
    pathogenic +1. Enough granularity to detect a monotone trend vs oscillation (stability
    factor) and a confirm-vs-revert trajectory (validation harness)."""
    s = (sig or "").lower()
    if "benign" in s:
        return -1
    if "pathogenic" in s and "conflict" not in s:
        return 1
    return 0


def stability_factor(sigs: list[str | None]) -> float:
    """How steady a variant's classification has been across snapshots (oldest → newest). A
    monotone strengthening toward pathogenic is more trustworthy (+0.5); ≥2 direction reversals
    means the call is flaky (−0.5); otherwise neutral. <2 data points → 0.0 (no trend)."""
    ords = [path_ordinal(s) for s in sigs]
    if len(ords) < 2:
        return 0.0
    diffs = [b - a for a, b in zip(ords, ords[1:]) if b != a]
    if not diffs:
        return 0.0
    reversals = sum(1 for x, y in zip(diffs, diffs[1:]) if (x > 0) != (y > 0))
    if reversals >= 2:
        return -0.5
    if all(d > 0 for d in diffs):   # steadily strengthening toward pathogenic
        return 0.5
    return 0.0


def _f_gene(in_panel: bool, constrained, is_null: bool) -> float:
    """Gene context. Panel membership is itself strong context — these are the curated,
    clinically-actionable monogenic disease genes — so it carries most of the weight; a
    LOF-constrained gene hit by a null variant (haploinsufficiency plausible) tops it up.
    Note: two-hit tumour suppressors (BRCA1/2) read as LOF-tolerant in gnomAD, so they rely on
    panel membership alone here — same constraint caveat the ACMG stage documents."""
    s = 0.7 if in_panel else 0.0
    if is_null and constrained:
        s += 0.3
    return min(1.0, s)


def score_delta(*, prev_sig, curr_sig, prev_rev, curr_rev, prev_af, curr_af,
                consequence=None, ref="", alt="", acmg_call=None,
                am=None, sa=None, en=None, in_panel=False, constrained=None,
                stability=0.0) -> dict | None:
    """Score one changed variant. Returns {tier, score, factors, change_kind, rationale} or
    None if nothing notable changed. Pure: same inputs → same verdict, fully reconstructable
    from the returned `factors`."""
    change = change_kind(prev_sig, curr_sig, prev_rev, curr_rev, prev_af, curr_af)
    if change is None:
        return None
    kind, desc = change

    is_null = bool(acmg_call) and acmg_call.get("csq") in (
        "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
        "stop_gained", "frameshift_variant", "start_lost")
    factors = {
        "acmg": round(_f_acmg(acmg_call, _clinvar_evidence(curr_sig, curr_rev)), 3),
        "model": round(_f_model(am, sa, en, ref, alt), 3),
        "freq": round(_f_freq(curr_af), 3),
        "stability": round(float(stability), 3),
        "gene": round(_f_gene(in_panel, constrained, is_null), 3),
    }
    w = config.TRIAGE_WEIGHTS
    score = sum(w[k] * factors[k] for k in factors)

    # `actionable_now` demands real clinical context, not just a high model/freq score on an
    # unknown gene — that lands at `monitor` instead.
    has_context = in_panel or is_actionable(curr_sig)
    cut = config.TRIAGE_CUTOFFS
    if score >= cut["actionable_now"] and has_context:
        tier = "actionable_now"
    elif score >= cut["monitor"]:
        tier = "monitor"
    elif score >= cut["log_only"]:
        tier = "log_only"
    else:
        tier = "ignore"

    top = max(factors, key=lambda k: abs(factors[k]))
    rationale = f"{kind.replace('_', ' ')}: {desc}; top factor {top}={factors[top]:+.2f}"
    return {"tier": tier, "score": round(score, 4), "factors": factors,
            "change_kind": kind, "rationale": rationale}
