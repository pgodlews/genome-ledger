"""Companion (derived) tables, all keyed on the canonical variant tuple.

These hold *derived* state — triage verdicts, provenance/confidence, per-snapshot
history, the incremental-scan digest ledger, and validation runs. Every one is
rebuildable from the immutable snapshots plus this code, so they never weaken the
reproducibility guarantee. They are created idempotently by `load.connect()` and
join the central `variants` table on (chrom, pos, ref, alt) [+ snapshot/sample].

`factors` is stored as a JSON *string* in a VARCHAR (no JSON-extension dependency):
auditable, greppable, and reconstructable without loading an extension.
"""

from __future__ import annotations

# One row per delta finding that triage scored (Feature 1).
TRIAGE_RESULTS = """
CREATE TABLE IF NOT EXISTS triage_results (
    snapshot_id     VARCHAR,   -- the "curr" snapshot of the diff
    prev_snapshot   VARCHAR,
    sample_id       VARCHAR,
    chrom           VARCHAR,
    pos             BIGINT,
    ref             VARCHAR,
    alt             VARCHAR,
    gene            VARCHAR,
    tier            VARCHAR,    -- actionable_now | monitor | log_only | ignore
    score           DOUBLE,     -- composite triage score (see triage.score_delta)
    factors         VARCHAR,    -- JSON: per-factor breakdown {acmg, model, freq, stability, gene}
    change_kind     VARCHAR,    -- new_pathogenic | upgraded | reclassified | downgraded | af_shift | model_shift
    rationale       VARCHAR,    -- human-readable one-liner
    created_at      TIMESTAMP DEFAULT now()
);
"""

# Source / version / confidence per (variant, source) (Feature 2). Scoped to
# variants that actually reach a report or triage — not all ~4M rows per sample.
ANNOTATION_PROVENANCE = """
CREATE TABLE IF NOT EXISTS annotation_provenance (
    snapshot_id     VARCHAR,
    sample_id       VARCHAR,
    chrom           VARCHAR,
    pos             BIGINT,
    ref             VARCHAR,
    alt             VARCHAR,
    source          VARCHAR,    -- clinvar | gnomad_af | alphamissense | spliceai | enformer | vep_csq | gnomad_constraint
    version         VARCHAR,    -- e.g. 'ClinVar 2026-06-28', 'VEP r112'
    value_summary   VARCHAR,    -- compact rendering of the value used
    confidence      VARCHAR,    -- High | Medium | Low | Very Low
    used_in_acmg    BOOLEAN,    -- did this value feed an ACMG evidence code?
    as_of           DATE
);
"""

# Per-snapshot timeline for each variant (Features 1-stability, 3, 4). Upserted
# every scan; the triage stability factor reads the recent window of this.
VARIANT_HISTORY = """
CREATE TABLE IF NOT EXISTS variant_history (
    chrom           VARCHAR,
    pos             BIGINT,
    ref             VARCHAR,
    alt             VARCHAR,
    snapshot_id     VARCHAR,
    snapshot_date   DATE,
    clinvar_sig     VARCHAR,
    clinvar_revstat VARCHAR,
    gnomad_af       DOUBLE,
    acmg_tier       VARCHAR,    -- from the acmg stage when the gene is on-panel
    model_band      VARCHAR,    -- enformer/borzoi band when scored
    PRIMARY KEY (chrom, pos, ref, alt, snapshot_id)
);
"""

# Digest ledger driving incremental scan (Feature 4): what input produced each
# sample's annotation, so an unchanged input can be copied forward instead of
# re-running VEP.
SCAN_INPUTS = """
CREATE TABLE IF NOT EXISTS scan_inputs (
    sample_id       VARCHAR,
    snapshot_id     VARCHAR,
    source          VARCHAR,    -- clinvar | vep_cache | normalized_vcf | alphamissense | spliceai | enformer
    digest          VARCHAR,    -- sha / version string that drove this annotation
    rows_loaded     BIGINT,
    as_of           TIMESTAMP DEFAULT now()
);
"""

# Retrospective validation (Feature 3): replay an old snapshot's triage, score it
# against a later snapshot's ClinVar.
VALIDATION_RUNS = """
CREATE TABLE IF NOT EXISTS validation_runs (
    run_id               VARCHAR,
    replay_snapshot      VARCHAR,
    truth_snapshot       VARCHAR,
    n_deltas             BIGINT,
    n_useful             BIGINT,
    n_noise              BIGINT,
    precision_actionable DOUBLE,
    created_at           TIMESTAMP DEFAULT now()
);
"""

VALIDATION_FINDINGS = """
CREATE TABLE IF NOT EXISTS validation_findings (
    run_id          VARCHAR,
    chrom           VARCHAR,
    pos             BIGINT,
    ref             VARCHAR,
    alt             VARCHAR,
    predicted_tier  VARCHAR,
    outcome         VARCHAR,    -- confirmed | reverted | stable_noise | unknown
    truth_sig       VARCHAR
);
"""

# Applied in order by load.connect() after the variants table.
COMPANION_TABLES = (
    TRIAGE_RESULTS,
    ANNOTATION_PROVENANCE,
    VARIANT_HISTORY,
    SCAN_INPUTS,
    VALIDATION_RUNS,
    VALIDATION_FINDINGS,
)
