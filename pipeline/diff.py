"""Stage 8: DELTA report — what changed between two snapshots, triaged.

This is the core deliverable. For every variant the family actually carries, it compares the
annotation under `prev` vs `curr` and, for each real change, runs the multi-factor triage
layer (`triage.score_delta`) to assign one of four tiers: actionable_now / monitor / log_only
/ ignore. Verdicts are persisted to `triage_results` (auditable, queryable, replayable) and
rendered grouped by tier. Because both runs used pinned, dated databases over identically
normalized genomes, a detected change is a real database change.
"""

from __future__ import annotations

import json
from pathlib import Path

from . import reportpaths
from . import (acmg, alphamissense, config, enformer, gnomad_constraint,
               limitations, spliceai, triage)
from .load import TABLE, connect, transaction
from .snapshot import load_manifest
from .util import disease_name, log, review_stars, write_report

# Tier → (badge, heading) for the rendered report. `ignore` is persisted but not shown.
_TIER_BADGE = {"actionable_now": "🔴", "monitor": "🟠", "log_only": "⚪", "ignore": "·"}


def _stability_window(con, curr: str, n: int) -> dict:
    """Per-variant ClinVar-significance trail (oldest → newest) over the last `n` snapshots up
    to `curr`, from variant_history. Empty when Phase-2 history hasn't been populated yet — the
    stability factor then degrades to neutral. Keyed by (chrom, pos, ref, alt)."""
    snaps = [r[0] for r in con.execute(
        "SELECT DISTINCT snapshot_id FROM variant_history WHERE snapshot_id <= ? "
        "ORDER BY snapshot_id DESC LIMIT ?", [curr, n]).fetchall()]
    if not snaps:
        return {}
    ph = ",".join("?" * len(snaps))
    window: dict = {}
    for chrom, pos, ref, alt, snap, sig in con.execute(
            f"SELECT chrom, pos, ref, alt, snapshot_id, clinvar_sig FROM variant_history "
            f"WHERE snapshot_id IN ({ph})", snaps).fetchall():
        window.setdefault((chrom, pos, ref, alt), []).append((snap, sig))
    return {k: [sig for _, sig in sorted(v)] for k, v in window.items()}  # snapshot-id sorts as date


def cache_mismatch(prev: str, curr: str) -> tuple[str, str] | None:
    """(prev_version, curr_version) when the two snapshots were annotated under different VEP
    cache versions, else None.

    A DELTA across a cache bump is not a database-change report. Only `curr` is re-annotated
    (`incremental.scan_mode` forces a full run on a bump), so `prev` keeps the OLD cache's
    gnomAD allele frequencies and the OLD `--pick` transcript choice. That produces two
    artefact classes the triage layer cannot tell from real change:

      * every variant whose gnomAD AF moved between cache releases becomes an `af_shift`
        (thousands of rows, none of them a ClinVar event);
      * every variant whose picked gene changed fails the gene equality in the join and
        appears TWICE — once curr-only (reads as a new finding) and once prev-only (reads
        as a downgrade).

    The ClinVar columns are unaffected: they come from the snapshot's own ClinVar VCF via
    VEP `--custom`, not from the cache. So the honest thing is to keep the ClinVar-driven
    verdicts and suppress the cache-driven artefacts, not to skip the report.
    """
    def _ver(sid: str):
        try:
            return load_manifest(sid).get("vep_cache_version")
        except (OSError, ValueError):
            # No/unreadable manifest (an old snapshot, or a synthetic pair in the validation
            # harness): we cannot tell, so don't claim a mismatch and don't take down the diff.
            return None

    pv, cv = _ver(prev), _ver(curr)
    if pv is None or cv is None or pv == cv:
        return None
    return (str(pv), str(cv))


def _persist(con, prev: str, curr: str, findings: list[dict]) -> None:
    """Replace this snapshot-pair's triage rows (idempotent re-runs)."""
    # Atomic: the DELETE alone would leave the pair with no triage rows at all if the INSERT
    # failed, which reads downstream as "this snapshot pair produced no findings".
    with transaction(con):
        con.execute("DELETE FROM triage_results WHERE snapshot_id = ? AND prev_snapshot = ?",
                    [curr, prev])
        if not findings:
            return
        con.executemany(
            "INSERT INTO triage_results (snapshot_id, prev_snapshot, sample_id, chrom, pos, "
            "ref, alt, gene, tier, score, factors, change_kind, rationale) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            [[curr, prev, f["sample"], f["chrom"], f["pos"], f["ref"], f["alt"], f["gene"],
              f["tier"], f["score"], json.dumps(f["factors"]), f["change_kind"],
              f["rationale"]] for f in findings])


def score_pair(con, prev: str, curr: str) -> list[dict]:
    """Triage every changed variant between two snapshots and return the findings — the shared
    core of the DELTA report and the retrospective validation harness (which replays this with
    *today's* code over an old snapshot pair). Pure read: persists nothing, renders nothing."""
    xcache = cache_mismatch(prev, curr)
    if xcache:
        log.warning(
            "DELTA %s → %s spans a VEP cache bump (%s → %s): gnomAD AF and --pick gene come "
            "from the cache, so af_shift rows and gene reassignments here are annotation "
            "artefacts, not database changes. Suppressing af_shift and joining without gene; "
            "ClinVar-driven verdicts are unaffected.", prev, curr, xcache[0], xcache[1])
    # Only diff samples annotated in BOTH snapshots. A sample present in just one (e.g. a
    # newly-sequenced family member) would otherwise flood the DELTA with false "new"
    # findings — its variants belong in the FULL report instead.
    common = [r[0] for r in con.execute(f"""
        SELECT sample_id FROM {TABLE} WHERE snapshot_id = ?
        INTERSECT
        SELECT sample_id FROM {TABLE} WHERE snapshot_id = ?
    """, [curr, prev]).fetchall()]
    new_samples = [r[0] for r in con.execute(f"""
        SELECT DISTINCT sample_id FROM {TABLE} WHERE snapshot_id = ?
          AND sample_id NOT IN (SELECT sample_id FROM {TABLE} WHERE snapshot_id = ?)
    """, [curr, prev]).fetchall()]
    if new_samples:
        log.info("DELTA skips %d sample(s) new since %s: %s (see their FULL report).",
                 len(new_samples), prev, ", ".join(new_samples))
    if not common:
        log.warning("No samples shared between snapshots %s and %s — no deltas.", prev, curr)
        return []

    placeholders = ",".join("?" * len(common))
    # Full outer join the two snapshots on (sample, variant, gene). Pull the curr-state
    # annotations triage scores on (consequence/sift/polyphen feed the ACMG + model factors).
    #
    # The WHERE clause is a cheap pre-gate for `triage.change_kind`, and it is what makes this
    # stage tractable on whole genomes: a snapshot pair shares essentially all of its rows, so
    # without it DuckDB ships every carried variant to Python (millions of rows per sample)
    # and the loop below runs the per-variant model + ACMG work on all of them. On two
    # consecutive snapshots with no annotation change at all, that was tens of thousands of
    # `tabix` subprocesses (one per `alphamissense`/`spliceai` lookup on a missense/splice
    # row) spent to produce an empty report. The same question
    # answered in SQL takes 0.1s.
    #
    # Deliberately a *superset* of change_kind's conditions: any difference in significance,
    # review status, or allele frequency passes, and score_delta stays the authority on
    # whether that difference is notable (it returns None otherwise). So this can only ever
    # discard rows change_kind would have rejected anyway. Rows where both sides are NULL for
    # all three are exactly change_kind's None case, including the FULL-OUTER-JOIN rows for a
    # variant present in only one snapshot with no ClinVar record on either side.
    # Across a VEP cache bump the picked gene can change for the same variant; joining on
    # gene would then split it into a phantom curr-only + prev-only pair. `--pick` gives one
    # row per variant, so dropping the gene equality is safe and keeps the pair together.
    gene_join = "" if xcache else "AND c.gene IS NOT DISTINCT FROM p.gene"
    rows = con.execute(f"""
        SELECT
            COALESCE(c.sample_id, p.sample_id)  AS sample_id,
            COALESCE(c.gene, p.gene)            AS gene,
            COALESCE(c.chrom, p.chrom)          AS chrom,
            COALESCE(c.pos, p.pos)              AS pos,
            COALESCE(c.ref, p.ref)              AS ref,
            COALESCE(c.alt, p.alt)              AS alt,
            p.clinvar_sig, c.clinvar_sig,
            p.clinvar_revstat, c.clinvar_revstat,
            p.gnomad_af, c.gnomad_af,
            c.zygosity, c.clinvar_disease,
            c.consequence, c.sift, c.polyphen
        FROM (SELECT * FROM {TABLE}
              WHERE snapshot_id = ? AND sample_id IN ({placeholders})) c
        FULL OUTER JOIN (SELECT * FROM {TABLE}
              WHERE snapshot_id = ? AND sample_id IN ({placeholders})) p
          ON  c.sample_id = p.sample_id AND c.chrom = p.chrom AND c.pos = p.pos
          AND c.ref = p.ref AND c.alt = p.alt {gene_join}
        WHERE c.clinvar_sig     IS DISTINCT FROM p.clinvar_sig
           OR c.clinvar_revstat IS DISTINCT FROM p.clinvar_revstat
           OR c.gnomad_af       IS DISTINCT FROM p.gnomad_af
    """, [curr, *common, prev, *common]).fetchall()
    log.info("DELTA %s → %s: %d candidate row(s) with an annotation difference to triage.",
             prev, curr, len(rows))

    disease_genes, gene_modes = acmg.panel_context()
    have_am, have_sa = alphamissense.available(), spliceai.available()
    have_en, have_gc = enformer.available(), gnomad_constraint.available()
    stability = _stability_window(con, curr, config.TRIAGE_STABILITY_WINDOW)

    findings: list[dict] = []
    n_af_suppressed = 0
    for (sample, gene, chrom, pos, ref, alt, ps, cs, prv, crv, paf, caf,
         zyg, disease, csq, sift, pp) in rows:
        # Ask what changed BEFORE doing any per-variant work. `change_kind` is pure and reads
        # only the six columns above, while `model_lookups` shells out to tabix and
        # `acmg.evaluate` runs the reasoner — neither is worth paying for on a row that
        # score_delta is about to discard. (score_delta re-derives this internally; it stays
        # the single source of truth for the verdict.)
        if triage.change_kind(ps, cs, prv, crv, paf, caf) is None:
            continue
        am, sa, en = triage.model_lookups(chrom, pos, ref, alt, csq,
                                           have_am, have_sa, have_en)
        acmg_call = acmg.evaluate(gene, chrom, pos, ref, alt, csq, cs, caf, sift, pp, zyg,
                                  disease_genes=disease_genes, gene_modes=gene_modes,
                                  have_am=have_am, have_sa=have_sa, have_gc=have_gc)
        in_panel = gene in disease_genes
        constrained = (gnomad_constraint.is_constrained(gene)
                       if have_gc and gene else None)
        verdict = triage.score_delta(
            prev_sig=ps, curr_sig=cs, prev_rev=prv, curr_rev=crv, prev_af=paf, curr_af=caf,
            consequence=csq, ref=ref, alt=alt, acmg_call=acmg_call,
            am=am, sa=sa, en=en, in_panel=in_panel, constrained=constrained,
            stability=triage.stability_factor(stability.get((chrom, pos, ref, alt), [])))
        if verdict is None:
            continue
        if xcache and verdict["change_kind"] == "af_shift":
            # gnomAD AF is a cache column: across a bump this is the cache moving, not the
            # database. Dropped rather than reported at a lower tier — there is no evidence
            # here to triage, and thousands of them would bury the real ClinVar findings.
            n_af_suppressed += 1
            continue
        findings.append({
            "sample": sample, "gene": gene, "chrom": chrom, "pos": pos, "ref": ref,
            "alt": alt, "zyg": zyg, "disease": disease, "stars": review_stars(crv),
            "prev_sig": ps, "curr_sig": cs,
            **verdict,
        })
    if n_af_suppressed:
        log.warning("DELTA %s → %s: %d af_shift row(s) suppressed as VEP-cache artefacts.",
                    prev, curr, n_af_suppressed)
    return findings


def diff(prev: str, curr: str) -> tuple[Path, int, str]:
    con = connect()
    findings = score_pair(con, prev, curr)
    _persist(con, prev, curr, findings)
    con.close()

    # Sort: tier (actionable first), then score desc, then review stars desc.
    order = {t: i for i, t in enumerate(triage.TIERS)}
    findings.sort(key=lambda f: (order[f["tier"]], -f["score"], -f["stars"],
                                 f["sample"], f["gene"] or ""))
    shown = [f for f in findings if f["tier"] != "ignore"]
    by_tier = {t: [f for f in findings if f["tier"] == t] for t in triage.TIERS}

    pmeta, cmeta = load_manifest(prev), load_manifest(curr)
    path = reportpaths.family_report(f"delta_{prev}_to_{curr}", curr)

    n_act = len(by_tier["actionable_now"])
    head = "| ! | Sample | Gene | Disease | Variant | Change | Score | Review | Zyg |"
    sep = "|---|--------|------|---------|---------|--------|:-----:|:------:|-----|"

    def _row(f) -> str:
        badge = _TIER_BADGE.get(f["tier"], "·")
        var = f"{f['chrom']}:{f['pos']} {f['ref']}>{f['alt']}"
        return (f"| {badge} | {f['sample']} | {f['gene'] or '?'} | "
                f"{disease_name(f['disease'])} | {var} | {f['rationale']} | "
                f"{f['score']:.2f} | {'★' * f['stars'] or '0★'} | {f['zyg'] or '-'} |")

    def _section(title: str, items: list[dict], blurb: str = "") -> list[str]:
        out = [f"## {title}", ""]
        if blurb:
            out += [blurb, ""]
        out += [head, sep]
        out += [_row(f) for f in items] or ["| | | _none_ | | | | | | |"]
        out += [""]
        return out

    xcache = cache_mismatch(prev, curr)
    banner = ([
        f"> ⚠️ **Annotation baseline changed between these snapshots** (VEP cache "
        f"{xcache[0]} → {xcache[1]}). gnomAD allele frequencies and the picked transcript "
        f"come from the cache, so frequency shifts here would be the cache moving rather "
        f"than the database: `af_shift` changes are suppressed and the variant join ignores "
        f"the gene, keeping a reassigned variant as one row instead of a phantom "
        f"new/downgraded pair. ClinVar-driven findings below are unaffected. For a clean "
        f"frequency comparison, re-annotate both snapshots under one cache version.",
        "",
    ] if xcache else [])

    lines = [
        f"# DELTA report — {prev} → {curr}",
        "",
        *banner,
        f"- **From:** ClinVar {pmeta['clinvar_date']}  →  **To:** ClinVar {cmeta['clinvar_date']}",
        f"- **Changes triaged:** {len(findings)}  |  "
        f"🔴 {n_act} actionable now  ·  🟠 {len(by_tier['monitor'])} monitor  ·  "
        f"⚪ {len(by_tier['log_only'])} log-only  ·  · {len(by_tier['ignore'])} ignored",
        "",
    ]
    if by_tier["actionable_now"]:
        lines += _section(
            "🔴 Actionable now",
            by_tier["actionable_now"],
            "_Highest combined evidence (ACMG + model + rarity + gene context) on a carried "
            "variant — review these first._")
    if by_tier["monitor"]:
        lines += _section("🟠 Monitor closely", by_tier["monitor"],
                          "_Material change worth watching; evidence not yet decisive._")
    if by_tier["log_only"]:
        lines += _section("⚪ Log only", by_tier["log_only"],
                          "_Recorded for the trail; low actionability._")
    if not shown:
        lines.append("_No notable changes since the previous snapshot._")
    lines += ["", "_🔴 actionable now · 🟠 monitor · ⚪ log only · score = weighted triage "
              "(ACMG · model · rarity · stability · gene); ★ = current ClinVar review "
              "confidence. Ignored changes are persisted in `triage_results`, not shown._", ""]
    lines += limitations.render(curr, report_kind="delta")
    write_report(path, lines)
    log.info("DELTA report: %s (%d triaged, %d actionable, %d monitor)",
             path, len(findings), n_act, len(by_tier["monitor"]))

    alert_summary = ""
    if by_tier["actionable_now"]:
        a = by_tier["actionable_now"][0]
        more = f" (+{n_act - 1} more)" if n_act > 1 else ""
        alert_summary = f"{n_act} ACTIONABLE: {a['gene'] or '?'} in {a['sample']}{more}"
    return path, len(findings), alert_summary
