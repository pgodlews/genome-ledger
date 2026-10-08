"""Phase 3: retrospective validation harness — were past deltas useful or noise?

Replays an old snapshot's triage with *today's* code (re-scoring the diff into `replay` via
`diff.score_pair`, so it reflects current rules/weights, not whatever was stored at the time),
then asks a later, retained snapshot what actually happened: did each flagged change hold up in
ClinVar, or get reverted? "Truth" is therefore later-ClinVar consensus — automatable from data
already on disk, and honestly *not* biological ground truth (the report says so).

Outcome per finding, from the variant's ClinVar pathogenicity trajectory replay → truth:
  * confirmed     — an upgrade prediction that stayed (likely) pathogenic at truth (useful)
  * reverted      — the change flipped back by truth (noise: premature)
  * stable_noise  — a low-signal change that never moved pathogenicity (noise)
  * unknown       — variant absent from the truth snapshot (survivorship; excluded from precision)

Headline metric: precision of the `actionable_now` tier = confirmed / (confirmed + reverted)
among the variants triage would have escalated. Persisted to validation_runs / _findings.

A run with zero deltas is a legitimate result, not an error: it means the replayed pair moved
no variant's ClinVar significance, so there is nothing to judge. It still records a
validation_runs row (n_deltas = 0, precision NULL) — the absence of a judgeable delta is part
of the audit trail, and a silent gap would be indistinguishable from a scan that never ran.
"""

from __future__ import annotations

from pathlib import Path

from . import reportpaths
from . import config, triage
from .diff import score_pair
from .load import TABLE, connect, transaction
from .util import log, write_report

_MISSING = object()  # variant not present at the truth snapshot (distinct from a NULL sig)


def _predecessor(con, snapshot_id: str) -> str | None:
    row = con.execute(
        f"SELECT DISTINCT snapshot_id FROM {TABLE} WHERE snapshot_id < ? "
        f"ORDER BY snapshot_id DESC LIMIT 1", [snapshot_id]).fetchone()
    return row[0] if row else None


def _outcome(change_kind: str, replay_sig, truth_sig) -> str:
    """Classify one prediction against the truth snapshot's ClinVar call."""
    if truth_sig is _MISSING:
        return "unknown"
    r, t = triage.path_ordinal(replay_sig), triage.path_ordinal(truth_sig)
    if change_kind in ("new_pathogenic", "upgraded"):     # predicted: heading pathogenic
        return "confirmed" if t >= 1 else "reverted"
    if change_kind == "downgraded":                        # predicted: heading benign/away
        return "confirmed" if t <= 0 else "reverted"
    # reclassified / af_shift: lower-signal. Useful only if ClinVar moved further the same way.
    if t > r:
        return "confirmed"
    if t < r:
        return "reverted"
    return "stable_noise"


def _classify(run_id: str, findings: list[dict], truth_sigs: dict) -> tuple[list, dict, dict]:
    """Classify each finding against `truth_sigs` (keyed by (chrom,pos,ref,alt)). Pure and
    testable in isolation from the DB. Returns (rows-for-insert, outcome counts, actionable-tier
    confirmed/reverted counts)."""
    rows, counts = [], {"confirmed": 0, "reverted": 0, "stable_noise": 0, "unknown": 0}
    act = {"confirmed": 0, "reverted": 0}  # for actionable-tier precision
    for f in findings:
        key = (f["chrom"], f["pos"], f["ref"], f["alt"])
        ts = truth_sigs.get(key, _MISSING)
        # Pass ts through unchanged (not None) so _outcome's `is _MISSING` check actually
        # fires — a variant absent from truth must classify as "unknown", not fall through
        # to the ordinal comparison (where None ordinal-ranks as 0, i.e. "benign", so a
        # new_pathogenic/upgraded prediction would be wrongly scored "reverted").
        outcome = _outcome(f["change_kind"], f["curr_sig"], ts)
        counts[outcome] += 1
        if f["tier"] == "actionable_now" and outcome in act:
            act[outcome] += 1
        rows.append([run_id, f["chrom"], f["pos"], f["ref"], f["alt"], f["tier"], outcome,
                     None if ts is _MISSING else ts])  # _MISSING is a sentinel, not storable
    return rows, counts, act


def validate(replay: str, truth: str) -> Path:
    con = connect()
    prev = _predecessor(con, replay)
    if prev is None:
        con.close()
        raise SystemExit(f"Cannot replay {replay}: no earlier snapshot loaded to diff against.")
    if truth <= replay:
        con.close()
        raise SystemExit(f"--truth ({truth}) must be later than --replay ({replay}).")

    findings = score_pair(con, prev, replay)  # today's triage over the historical pair
    truth_sigs = {(c, p, r, a): s for c, p, r, a, s in con.execute(
        f"SELECT chrom, pos, ref, alt, any_value(clinvar_sig) FROM {TABLE} "
        f"WHERE snapshot_id = ? GROUP BY chrom, pos, ref, alt", [truth]).fetchall()}

    run_id = f"{replay}_vs_{truth}"
    rows, counts, act = _classify(run_id, findings, truth_sigs)

    n = len(findings)
    n_useful = counts["confirmed"]
    n_noise = counts["reverted"] + counts["stable_noise"]
    denom = act["confirmed"] + act["reverted"]
    precision = (act["confirmed"] / denom) if denom else None

    # Atomic across both tables: a run row without its findings (or findings without their
    # run) is a corrupt audit trail, and this is the harness that judges the triage model.
    with transaction(con):
        con.execute("DELETE FROM validation_findings WHERE run_id = ?", [run_id])
        con.execute("DELETE FROM validation_runs WHERE run_id = ?", [run_id])
        con.execute(
            "INSERT INTO validation_runs (run_id, replay_snapshot, truth_snapshot, n_deltas, "
            "n_useful, n_noise, precision_actionable) VALUES (?,?,?,?,?,?,?)",
            [run_id, replay, truth, n, n_useful, n_noise, precision])
        if rows:  # DuckDB's executemany rejects an empty parameter list outright
            con.executemany(
                "INSERT INTO validation_findings (run_id, chrom, pos, ref, alt, "
                "predicted_tier, outcome, truth_sig) VALUES (?,?,?,?,?,?,?,?)", rows)
    con.close()

    path = _render(replay, prev, truth, n, counts, act, precision, findings, rows)
    log.info("VALIDATION %s: %d deltas — %d confirmed, %d reverted, %d stable-noise, %d unknown "
             "(actionable precision %s)", run_id, n, counts["confirmed"], counts["reverted"],
             counts["stable_noise"], counts["unknown"],
             f"{precision:.0%}" if precision is not None else "n/a")
    return path


def _render(replay, prev, truth, n, counts, act, precision, findings, rows) -> Path:
    path = reportpaths.family_report(f"validation_{replay}_vs_{truth}", replay)
    by_tier = {}
    for r in rows:
        by_tier.setdefault(r[5], {}).setdefault(r[6], 0)
        by_tier[r[5]][r[6]] += 1

    prec_s = f"**{precision:.0%}**" if precision is not None else "_n/a (no actionable_now alerts)_"
    lines = [
        f"# Validation replay — triage at {replay}, judged against ClinVar {truth}",
        "",
        f"- **Replayed delta:** {prev} → {replay} (re-scored with today's triage rules).",
        f"- **Truth:** ClinVar as of snapshot {truth}.",
        f"- **Actionable-now precision:** {prec_s}  "
        f"({act['confirmed']} confirmed / {act['confirmed'] + act['reverted']} judged).",
        f"- **All deltas ({n}):** {counts['confirmed']} confirmed · {counts['reverted']} "
        f"reverted · {counts['stable_noise']} stable-noise · {counts['unknown']} unknown.",
        "",
        "## Outcome by predicted tier",
        "",
        "| Tier | Confirmed | Reverted | Stable-noise | Unknown |",
        "|------|:---------:|:--------:|:------------:|:-------:|",
    ]
    for tier in triage.TIERS:
        c = by_tier.get(tier, {})
        if not c:
            continue
        lines.append(f"| {tier} | {c.get('confirmed', 0)} | {c.get('reverted', 0)} | "
                     f"{c.get('stable_noise', 0)} | {c.get('unknown', 0)} |")
    lines += [
        "",
        "## Method & caveats",
        "- **Replayed with today's code**, so this measures the *current* rules/weights against "
        "history — re-run after tuning `config.TRIAGE_*` to see the effect.",
        "- **\"Truth\" is later-ClinVar consensus, not biological ground truth.** A reverted call "
        "may simply reflect ClinVar churn; a confirmed one, durable submitter agreement.",
        "- **Survivorship:** variants absent from the truth snapshot are `unknown` and excluded "
        "from precision.",
        "",
    ]
    write_report(path, lines)
    return path
