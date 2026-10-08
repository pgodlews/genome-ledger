"""Auto-generated "Assumptions & Known Limitations" block, appended to every report.

Sourced from the run's pinned manifest + actual inputs (model availability, panel size,
family composition) so a caveat can never silently contradict reality — if AlphaMissense is
absent or the family has no child genome, the block says so. This is the truth-over-comfort
principle made structural: every report carries its own disclosure of what it can and cannot
see. Best-effort throughout — a missing input degrades to a plain statement, never an error.
"""

from __future__ import annotations

from . import alphamissense, config, enformer, spliceai
from .panel import load_panel
from .snapshot import load_manifest
from .util import log, read_tsv


def _on(flag: bool) -> str:
    return "on" if flag else "absent"


def _family_note() -> str:
    """Whether a child/trio genome is present — gates de-novo and phasing inferences."""
    try:
        fam = read_tsv(config.FAMILY_FILE)
        n = len(read_tsv(config.SAMPLES_TSV)) or len(fam)
        has_child = any((r.get("father") or "").strip() or (r.get("mother") or "").strip()
                        for r in fam)
        if has_child:
            return f"{n} sample(s); a child genome is present (de-novo / phasing possible)."
        return (f"{n} sample(s); **no child genome** — de-novo calling and read-backed "
                "phasing are unavailable.")
    except Exception as e:  # noqa: BLE001 — best-effort
        log.debug("limitations: family note unavailable: %s", e)
        return "family composition unavailable."


def render(snapshot_id: str, report_kind: str = "full", sample: str | None = None) -> list[str]:
    """Return the markdown lines for the limitations section. `report_kind` ∈
    {full, delta, acmg, panel} tailors a closing note. Never raises."""
    try:
        meta = load_manifest(snapshot_id)
    except Exception as e:  # noqa: BLE001
        log.debug("limitations: manifest unavailable for %s: %s", snapshot_id, e)
        meta = {}

    have_am, have_sa = alphamissense.available(), spliceai.available()
    have_en = enformer.available()
    calib = " (calibrated)" if (have_en and enformer.calibration_available()) else \
            (" (uncalibrated)" if have_en else "")
    try:
        n_genes = len({g for d in load_panel() for g in d.genes})
    except Exception:  # noqa: BLE001
        n_genes = 0

    lines = [
        "## Assumptions & Known Limitations",
        "",
        "_Auto-generated from this run's pinned versions and inputs — read alongside every "
        "finding._",
        "",
        f"- **Snapshot:** `{snapshot_id}` · build {meta.get('build', '?')} · "
        f"ClinVar {meta.get('clinvar_date', '?')} · VEP cache r{meta.get('vep_cache_version', '?')} "
        f"· gnomAD {meta.get('gnomad_version', '?')}.",
        "- **ClinVar is consensus, not ground truth.** A classification reflects current "
        "submitter agreement; it can move, and the absence of a pathogenic assertion is not "
        "evidence of benignity.",
        f"- **In-silico models:** AlphaMissense {_on(have_am)} (missense only) · "
        f"SpliceAI {_on(have_sa)} (splice only) · {config.INTERPRET_MODEL} non-coding "
        f"regulatory {_on(have_en)}{calib} (~{config.ENFORMER_SEQ_LEN:,} bp context; cannot "
        "speak to coding-sequence consequence). Scores are predictions, not measurements.",
        f"- **ACMG scope:** reasoned classification covers coding/splice variants in the "
        f"{n_genes} monogenic-panel genes only; non-panel and non-coding variants are not "
        "ACMG-classified.",
        f"- **Callability:** key-gene coverage is reported per sample (a base counts as "
        f"callable at ≥{config.CALLABLE_MIN_DEPTH}×; genes below "
        f"{config.CALLABLE_MIN_FRAC:.0%} callable are flagged in `callable_<sample>.md`).",
        f"- **Allele balance:** a heterozygous call where fewer than "
        f"{config.LOW_VAF_THRESHOLD:.0%} of reads carry the variant is marked "
        "\"⚠️ low allele fraction — review\" in the FULL and ACMG reports and \"low AF\" on "
        "the inheritance page. A real germline het sits near 50%; below the line is the "
        "usual signature of sequencer slippage (homopolymers) or mismapped reads. Flagged, "
        "never removed — at ~45× a true het falls below it about 1 time in 2,500. It does "
        "not catch paralog artifacts that sit near 50%; those need region-level review.",
        "- **Variant calling filters:** our own calls (bare sample id) are GATK "
        "HaplotypeCaller output, emitted at QUAL ≥ 30 (≤ 1-in-1,000 model error) and "
        "hard-filtered with GATK's standard thresholds (e.g. QD < 2.0 fails); `_vendor` / "
        "`_t2t` callsets carry their source pipeline's own filters. Either model's error "
        "estimate assumes independent read errors and is overconfident in homopolymers and "
        "duplicated regions.",
        f"- **Family:** {_family_note()}",
    ]
    if report_kind == "delta":
        lines.append(
            "- **Triage is a heuristic.** Tiers come from a weighted, adjustable score "
            "(`config.TRIAGE_*`); they rank attention, they do not diagnose. The full "
            "per-factor breakdown for every change is in the `triage_results` table.")
    lines += ["- **Research-grade, not a clinical diagnosis.**", ""]
    return lines
