#!/usr/bin/env python3
"""Genome Ledger — single CLI entry point.

    python run.py setup
    python run.py ingest <vcf> --sample dad --relation father --sex male
    python run.py normalize dad
    python run.py snapshot
    python run.py scan                      # the periodic entry point
    python run.py check-releases            # run by launchd; scans only if DB changed
    python run.py query queries/de_novo.sql --param child=kid1 --param father=dad ...

See README.md for the full workflow.
"""

from __future__ import annotations

import argparse
import sys

from pipeline import (ancestry as ancestry_mod, annotate as annotate_mod,
                      call as call_mod, denovo as denovo_mod, diff as diff_mod,
                      segregation as segregation_mod,
                      ingest as ingest_mod, load as load_mod, normalize as norm_mod,
                      panel as panel_mod, pgx as pgx_mod, releases,
                      report as report_mod, setup as setup_mod,
                      snapshot as snap_mod, traits as traits_mod,
                      hla as hla_mod, prs as prs_mod, repeats as repeats_mod,
                      sv as sv_mod, callability as callability_mod, acmg as acmg_mod,
                      phasing as phasing_mod, forcecall as forcecall_mod,
                      mito as mito_mod, family_migrate as family_migrate_mod,
                      interpret as interpret_mod, provenance as provenance_mod,
                      validate as validate_mod, incremental as incr_mod,
                      consent as consent_mod, hla_typing as hla_typing_mod,
                      notify as notify_mod, share as share_mod)
from pipeline.util import log, read_tsv, setup_logging
from pipeline import config


def _all_samples() -> list[str]:
    return [r["sample_id"] for r in read_tsv(config.SAMPLES_TSV)]


def _build_dashboards(snap: str | None = None) -> None:
    """Per-PERSON dashboards + the family inheritance tree, into the snapshot's dated dir,
    plus the two index.html landing pages (dated dir and reports root).

    Iterates people, not samples: `Adam` and `Adam_vendor` are one person's two
    callsets and share one page (pipeline/render.py:person_reports).

    `snap` defaults to the latest snapshot. Naming an older one rebuilds that dated
    directory's pages in place — which is also how a machine holding report files synced
    from elsewhere gives them their dashboards without having scanned them itself. The
    inheritance tree is the one part that needs the variant rows, so it is skipped (with a
    warning, never fatally) when this database has nothing for that snapshot."""
    from pipeline import pedigree, reportpaths
    from pipeline.render import (build_index, build_landing, build_snapshot_index,
                                 current_snapshot)

    # md/ holds the latest snapshot only; anything older moves to md-history/ first, so
    # the pages built below link each report where it now is.
    reportpaths.archive_older(current_snapshot())
    fam = pedigree.load(sequenced=pedigree.sequenced_ids())
    # The inheritance tree goes first: the pages built after it link to it only when it is
    # on disk, so building it last made every one of them omit the link on a fresh run.
    try:
        from pipeline import inheritance
        inheritance.build_page(snap)
    # SystemExit too: it is this repo's idiom for a user-facing precondition failure, and
    # it does NOT derive from Exception — an uncaught one here took the snapshot index and
    # the landing page down with it, which is exactly what this guard exists to prevent.
    except (Exception, SystemExit) as e:  # noqa: BLE001 — never fail the scan
        log.warning("inheritance tree skipped: %s", e)
    for pid in fam.people:
        if fam.is_sequenced(pid):
            build_index(pid, snap=snap)
    build_snapshot_index(snap)
    # The pages just written list the older snapshots; the older snapshots' pages were
    # written before this one existed, so bring their lists up to date as well.
    from pipeline.render import refresh_other_snapshots
    refresh_other_snapshots(snap or current_snapshot())
    build_landing()


def _refresh_person_pages(sample: str) -> None:
    """Rebuild a person's dashboard page in every dated directory that has one.

    A revocation deletes the reports that disclosed the finding, but the dashboard had
    already copied a summary out of them (the "highest polygenic score" tile), and that page
    is a file of its own. Rebuilt from what is on disk now, it no longer has anything to
    quote. Every dated directory, because older snapshots' pages show the same undated
    per-sample reports."""
    from pipeline import pedigree, reportpaths
    from pipeline.render import build_index
    try:
        person = pedigree.load().person_of(sample)
        root = reportpaths.html_root()
        dated = [d for d in sorted(root.iterdir())
                 if (d / f"index_{person}.html").exists()] if root.is_dir() else []
        for d in dated:
            build_index(person, snap=d.name)
    except (Exception, SystemExit) as e:  # noqa: BLE001 — the revocation itself has succeeded
        log.warning("Could not rebuild %s's dashboard pages (%s) — they may still show a "
                    "summary of the removed report(s). Run `run.py render` to rebuild them.",
                    sample, e)


def cmd_scan(args) -> int:
    """Periodic chain: snapshot → (incremental) annotate+load all → report → diff(prev).

    Incremental by default: a sample whose inputs (genome, ClinVar, VEP cache) are unchanged
    since the previous snapshot is copied forward instead of re-annotated. A full re-annotation
    is forced by `--full`, a VEP cache bump, or the monthly anchor (see `incremental.scan_mode`)."""
    curr = snap_mod.create_snapshot()
    # Strictly older than `curr` — `create_snapshot` reuses an existing same-day id, so
    # `latest_snapshot()` here would return `curr` itself on a second scan in one day.
    prev = snap_mod.previous_snapshot(curr)
    samples = _all_samples()
    if not samples:
        log.warning("No samples ingested yet — nothing to scan.")
        return 1
    explicit_full = getattr(args, "full", False)
    full, reason = incr_mod.scan_mode(prev, curr, force_full=explicit_full)
    log.info("Scan %s: %s.", curr, "FULL re-annotation — " + reason if full else reason)
    tally = {"copy_forward": 0, "clinvar_delta": 0, "full": 0}
    done, skipped = [], []
    for s in samples:
        if not norm_mod.normalized_path(s).exists():
            log.warning("%s not normalized; skipping. Run `normalize %s`.", s, s)
            skipped.append(s)
            continue
        kind = "full" if full else incr_mod.incremental_kind(s, prev, curr)
        if kind == "copy_forward":
            rows = incr_mod.copy_forward(s, prev, curr)
        elif kind == "clinvar_delta":
            rows = incr_mod.reannotate_clinvar(s, prev, curr)
        else:
            # Force only when a human passed --full. `annotate` now re-annotates on its own
            # whenever an input moved, so the automatic monthly anchor need not discard work
            # it can verify: a stamp-matching TSV is byte-identical to what a re-run would
            # produce, which is the whole point of the anchor. Forcing there would make every
            # retried anchor scan re-run VEP over the samples it had already finished.
            annotate_mod.annotate(s, curr, force=explicit_full)
            rows = load_mod.load(s, curr)
        tally[kind] += 1
        incr_mod.record_inputs(s, curr, incr_mod.source_digests(s, curr), rows)
        done.append(s)
    if not done:
        # Nothing was analysed, so there is nothing to report on and nothing to mark: a
        # completion marker here would tell `check-releases` this ClinVar release is handled.
        log.error("Scan %s: none of the %d ingested sample(s) is normalized — nothing "
                  "scanned. Run `normalize <sample>` first.", curr, len(samples))
        return 1
    if full:
        incr_mod.mark_anchor(curr)  # this snapshot is the month's reproducibility anchor
    log.info("Scan %s: %d copied forward, %d ClinVar-Δ re-annotated, %d full.",
             curr, tally["copy_forward"], tally["clinvar_delta"], tally["full"])
    report_mod.full_report(curr)
    panel_mod.panel_report(curr)
    acmg_mod.acmg_report(curr)
    # Provenance + history *before* the diff, so the triage stability factor sees curr.
    provenance_mod.record(curr)

    if prev and prev != curr:
        path, n, alert = diff_mod.diff(prev, curr)
        msg = f"Genome scan {curr}: {n} change(s) vs {prev}. Report: {path.name}"
        if alert:
            msg = f"⚠️ {alert} — {msg}"
            log.warning(msg)
        else:
            log.info(msg)
    else:
        log.info("No previous snapshot to diff against — FULL report only.")
        msg = f"Genome scan {curr}: baseline established ({len(done)} sample(s))."
        log.info(msg)
    if skipped:
        # Said in the notification, not only the log: the scan is complete for the others,
        # and a reader of "0 changes" must not take that to cover these.
        msg += (f" NOT scanned (no normalized VCF): {', '.join(skipped)} — "
                "normalize them and the next check-releases will pick them up.")

    _build_dashboards()
    # Last statement in the chain: everything above completed, so this ClinVar release can
    # now be treated as processed. Anything that raises earlier leaves the snapshot unmarked
    # and `check-releases` re-runs it instead of standing the weekly job down. The marker
    # names who was scanned, so a skipped sample is not counted as processed either.
    snap_mod.mark_scan_complete(curr, done, skipped)
    notify_mod.notify(msg)
    return 0


def _rescore_stale_prs() -> int:
    """Re-score the callsets whose PRS report predates a change to the trait panel.

    Polygenic scores do not depend on ClinVar, so the scan never touches them; what makes
    one stale is the panel (`panels/complex_traits.tsv`) gaining or changing a score. The
    report carries a stamp of the panel it was scored against (`prs.panel_token`), and this
    is where the scheduled job notices the mismatch — registering a score is an edit to a
    file, so there is no command to hang an immediate rescore on.

    One callset failing does not stop the rest. Returns 1 if any failed."""
    if not config.PRS_AUTO_RESCORE:
        return 0
    stale = prs_mod.stale_reports(_all_samples())
    if not stale:
        return 0
    log.info("PRS: %d callset(s) were scored against an older trait panel — re-scoring: %s",
             len(stale), ", ".join(stale))
    failed = []
    for s in stale:
        try:
            prs_mod.prs_report(s)
        except (Exception, SystemExit) as e:  # noqa: BLE001 — isolate each callset
            log.error("PRS re-score failed for %s: %s", s, e)
            failed.append(s)
    try:
        _build_dashboards()
    except (Exception, SystemExit) as e:  # noqa: BLE001 — the reports themselves are written
        log.warning("Dashboards not rebuilt after the PRS re-score (%s); run `render`.", e)
    msg = (f"Polygenic scores: {len(stale) - len(failed)} callset(s) re-scored against the "
           "updated trait panel.")
    if failed:
        msg += f" FAILED: {', '.join(failed)} — they stay stale and are retried next time."
    notify_mod.notify(msg)
    return 1 if failed else 0


def cmd_check_releases(args) -> int:
    rc = 0
    if releases.has_new_release():
        log.info("New database release detected — running scan.")
        rc = cmd_scan(args)
    else:
        log.info("Databases unchanged — nothing to do.")
    # After the scan, never instead of it: the two share the CPUs and the dashboards.
    return _rescore_stale_prs() or rc


def cmd_query(args) -> int:
    import duckdb
    bad = [p for p in (args.param or []) if "=" not in p]
    if bad:
        # dict() on a 1-element split raises a bare "dictionary update sequence element #0
        # has length 1" ValueError, which tells the user nothing about --param.
        raise SystemExit(f"--param must be name=value; got: {', '.join(repr(b) for b in bad)}")
    params = dict(p.split("=", 1) for p in (args.param or []))
    sql = open(args.sqlfile).read()
    for k, v in params.items():
        sql = sql.replace(f"${{{k}}}", v)
    con = duckdb.connect(str(config.DUCKDB_FILE), read_only=True)
    rows = con.execute(sql).fetchall()
    cols = [d[0] for d in con.description]
    print("\t".join(cols))
    for r in rows:
        print("\t".join("" if v is None else str(v) for v in r))
    con.close()
    return 0


def cmd_render(args) -> int:
    """(Re)generate HTML from existing markdown reports + per-person dashboards,
    without re-running any expensive stage. Only reads *.md (never touches the
    PharmCAT HTML)."""
    from pipeline import reportpaths
    from pipeline.render import render_formats
    n = 0
    for md in reportpaths.all_md():
        n += len(render_formats(md))
    _build_dashboards(getattr(args, "snapshot", None))
    log.info("render: wrote %d rendered report file(s) + dashboards.", n)
    return 0


def cmd_serve(args) -> int:
    from pipeline import serve as serve_mod
    serve_mod.serve(args.host, args.port)
    return 0


def cmd_demo(args) -> int:
    """Run a zero-data end-to-end demonstration using the synthetic family fixtures.

    Every page it writes is marked as demo output (config.DEMO_NOTICE), and the settings it
    changes for that are put back, so a caller in the same process is left as it was."""
    from pipeline import demo as demo_mod
    saved = (config.DEMO_NOTICE, config.PRS_S1_GATE)
    config.DEMO_NOTICE = demo_mod.NOTICE
    try:
        return _run_demo(args)
    finally:
        config.DEMO_NOTICE, config.PRS_S1_GATE = saved


def _run_demo(args) -> int:
    from pipeline import demo as demo_mod
    import json
    from pathlib import Path
    import shutil
    from pipeline.util import today

    demo_root = Path(args.out_dir).resolve()
    if args.clean and demo_root.exists():
        log.info("Cleaning existing demo directory: %s", demo_root)
        shutil.rmtree(demo_root)

    log.info("================================================================================")
    log.info("Starting Genome Ledger Demo (Zero-Data Mode)")
    log.info("Output directory: %s", demo_root)
    log.info("================================================================================")

    config.set_genomes_root(demo_root)
    for d in (config.REPORTS_DIR, config.MANIFEST_DIR, config.FAMILY_FILE.parent,
              config.RAW_DIR, config.NORMALIZED_DIR, config.CALLED_DIR,
              config.DUCKDB_DIR, config.SNAPSHOTS_DIR, config.TRAITS_DIR):
        d.mkdir(parents=True, exist_ok=True)

    fixtures_dir = Path(__file__).resolve().parent / "tests" / "fixtures" / "trio"
    if not fixtures_dir.exists():
        log.error("Trio fixtures directory not found: %s", fixtures_dir)
        return 1

    # 1. Family pedigree setup (3 generations, 6 members)
    log.info("1/6 Setting up family graph and pedigree...")
    shutil.copy2(fixtures_dir / "family.tsv", config.FAMILY_FILE)
    from pipeline import pedigree as pedigree_mod
    pedigree_mod.generate_pedigree_ped()

    # 2. Ingest synthetic family cohort samples
    samples_info = [
        ("father", "father", "male"),
        ("mother", "mother", "female"),
        ("daughter", "daughter", "female"),
        ("partner", "partner", "male"),
        ("grandson", "grandson", "male"),
        ("granddaughter", "granddaughter", "female"),
    ]
    log.info("2/6 Ingesting synthetic family fixtures (6 members, 3 generations)...")
    for s, rel, sex in samples_info:
        src_vcf = fixtures_dir / f"{s}.vcf.gz"
        ingest_mod.ingest(str(src_vcf), s, rel, sex, warn_unmirrored=False)
        # Populate called and normalized dirs so downstream genotype lookups resolve
        shutil.copy2(src_vcf, config.CALLED_DIR / f"{s}.GRCh38.vcf.gz")
        shutil.copy2(fixtures_dir / f"{s}.vcf.gz.tbi", config.CALLED_DIR / f"{s}.GRCh38.vcf.gz.tbi")
        shutil.copy2(src_vcf, config.NORMALIZED_DIR / f"{s}.norm.vcf.gz")
        shutil.copy2(fixtures_dir / f"{s}.vcf.gz.tbi", config.NORMALIZED_DIR / f"{s}.norm.vcf.gz.tbi")

    # Also populate child alias if child.vcf.gz exists for backwards compatibility
    child_vcf = fixtures_dir / "child.vcf.gz"
    if child_vcf.exists():
        shutil.copy2(child_vcf, config.CALLED_DIR / "child.GRCh38.vcf.gz")
        shutil.copy2(fixtures_dir / "child.vcf.gz.tbi", config.CALLED_DIR / "child.GRCh38.vcf.gz.tbi")
        shutil.copy2(child_vcf, config.NORMALIZED_DIR / "child.norm.vcf.gz")
        shutil.copy2(fixtures_dir / "child.vcf.gz.tbi", config.NORMALIZED_DIR / "child.norm.vcf.gz.tbi")

    # The family is synthetic, so the demo shows everything: gate off, and every member
    # recorded as opted in so the consent-aware stages agree if they are run by hand later.
    config.PRS_S1_GATE = False
    for s, _, _ in samples_info:
        consent_mod.set_consent(s, consent_mod.CATEGORY_ADULT_ONSET_UNTREATABLE, True)

    # 3. Single-SNP traits & wellness reports
    log.info("3/6 Generating single-SNP traits reports...")
    for s, _, _ in samples_info:
        traits_mod.traits_report(s, force=True)

    # 4. The per-person engines. HLA tag SNPs are read from the demo genomes for real; the
    # stages that need reads or reference downloads get bundled sample results instead,
    # run through the real report code and labelled as such (pipeline/demo.py).
    log.info("4/6 HLA tag SNPs (real), plus illustrative repeats / PGx / haplogroup / PRS reports...")
    for s, _, _ in samples_info:
        hla_mod.hla_report(s, force=True)
    demo_mod.toy_reports([s for s, _, _ in samples_info], fixtures_dir.parent / "demo")

    # 5. DuckDB & Monogenic disease panel screening
    log.info("5/6 Populating analytical database and rendering disease panel reports...")
    snapshot_id = f"{today()}-demo"
    snap_dir = config.SNAPSHOTS_DIR / snapshot_id
    snap_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = snap_dir / "manifest.json"
    manifest_data = {
        "snapshot_id": snapshot_id,
        "created": today(),
        "build": config.BUILD,
        "clinvar_url": "synthetic-demo-fixtures",
        "clinvar_date": today(),
        "vep_cache_version": config.VEP_CACHE_VERSION,
        "gnomad_version": config.GNOMAD_VERSION,
    }
    manifest_path.write_text(json.dumps(manifest_data, indent=2))

    con = load_mod.connect()
    con.execute(f"DELETE FROM {load_mod.TABLE} WHERE snapshot_id = ?", [snapshot_id])
    records = [
        # G1 Father: APOB het (At-Risk), HFE het (Carrier)
        (snapshot_id, "father", "2", 21002409, "G", "A", "APOB", "missense_variant",
         "Pathogenic", "reviewed_by_expert_panel", "17822",
         "Familial hypercholesterolemia", 0.0005, "deleterious(0.01)", "probably_damaging(0.99)", "HET"),
        (snapshot_id, "father", "6", 26093141, "G", "A", "HFE", "missense_variant",
         "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", "24440",
         "Hemochromatosis type 1", 0.05, "deleterious(0.01)", "probably_damaging(0.98)", "HET"),

        # G1 Mother: HFE het (Carrier), CFTR het (Carrier)
        (snapshot_id, "mother", "6", 26093141, "G", "A", "HFE", "missense_variant",
         "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", "24440",
         "Hemochromatosis type 1", 0.05, "deleterious(0.01)", "probably_damaging(0.98)", "HET"),
        (snapshot_id, "mother", "7", 117548628, "G", "A", "CFTR", "stop_gained",
         "Pathogenic", "reviewed_by_expert_panel", "7134",
         "Cystic fibrosis", 0.001, "deleterious(0.0)", "probably_damaging(1.0)", "HET"),

        # G2 Daughter: APOB het (At-Risk), HFE hom (Affected), CFTR het (Carrier)
        (snapshot_id, "daughter", "2", 21002409, "G", "A", "APOB", "missense_variant",
         "Pathogenic", "reviewed_by_expert_panel", "17822",
         "Familial hypercholesterolemia", 0.0005, "deleterious(0.01)", "probably_damaging(0.99)", "HET"),
        (snapshot_id, "daughter", "6", 26093141, "G", "A", "HFE", "missense_variant",
         "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", "24440",
         "Hemochromatosis type 1", 0.05, "deleterious(0.01)", "probably_damaging(0.98)", "HOM"),
        (snapshot_id, "daughter", "7", 117548628, "G", "A", "CFTR", "stop_gained",
         "Pathogenic", "reviewed_by_expert_panel", "7134",
         "Cystic fibrosis", 0.001, "deleterious(0.0)", "probably_damaging(1.0)", "HET"),

        # G2 Partner (Marry-in Founder): HFE het (Carrier)
        (snapshot_id, "partner", "6", 26093141, "G", "A", "HFE", "missense_variant",
         "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", "24440",
         "Hemochromatosis type 1", 0.05, "deleterious(0.01)", "probably_damaging(0.98)", "HET"),

        # G3 Grandson: APOB het (At-Risk), HFE hom (Affected), CFTR het (Carrier), PEX10 (de novo)
        (snapshot_id, "grandson", "2", 21002409, "G", "A", "APOB", "missense_variant",
         "Pathogenic", "reviewed_by_expert_panel", "17822",
         "Familial hypercholesterolemia", 0.0005, "deleterious(0.01)", "probably_damaging(0.99)", "HET"),
        (snapshot_id, "grandson", "6", 26093141, "G", "A", "HFE", "missense_variant",
         "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", "24440",
         "Hemochromatosis type 1", 0.05, "deleterious(0.01)", "probably_damaging(0.98)", "HOM"),
        (snapshot_id, "grandson", "7", 117548628, "G", "A", "CFTR", "stop_gained",
         "Pathogenic", "reviewed_by_expert_panel", "7134",
         "Cystic fibrosis", 0.001, "deleterious(0.0)", "probably_damaging(1.0)", "HET"),
        (snapshot_id, "grandson", "1", 10000000, "A", "T", "PEX10", "missense_variant",
         "Uncertain_significance", "criteria_provided,_single_submitter", "999999",
         "Peroxisome biogenesis disorder", None, "deleterious(0.02)", "probably_damaging(0.95)", "HET"),

        # G3 Granddaughter: HFE het (Carrier) — normal APOB and normal CFTR
        (snapshot_id, "granddaughter", "6", 26093141, "G", "A", "HFE", "missense_variant",
         "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", "24440",
         "Hemochromatosis type 1", 0.05, "deleterious(0.01)", "probably_damaging(0.98)", "HET"),
    ]
    con.executemany(f"""
        INSERT INTO {load_mod.TABLE} (
            snapshot_id, sample_id, chrom, pos, ref, alt, gene, consequence,
            clinvar_sig, clinvar_revstat, clinvar_alleleid, clinvar_disease,
            gnomad_af, sift, polyphen, zygosity
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, records)

    # An earlier snapshot of the same family, so the DELTA report — the pipeline's core
    # deliverable — has something to compare. Same people and genotypes; only what the
    # variant database said about three of their variants differs, the way ClinVar moves
    # between releases. The comparison itself is the real diff + triage code.
    from datetime import date, timedelta
    prev_id = f"{(date.today() - timedelta(days=90)).isoformat()}-demo"
    earlier = {  # gene -> (clinvar_sig, clinvar_revstat) as of the earlier snapshot
        "APOB": ("Uncertain_significance", "criteria_provided,_single_submitter"),
        "CFTR": ("Pathogenic", "criteria_provided,_single_submitter"),
        "PEX10": ("Likely_benign", "criteria_provided,_single_submitter"),
    }
    prev_records = []
    for r in records:
        r = list(r)
        r[0] = prev_id
        if r[6] in earlier:
            r[8], r[9] = earlier[r[6]]
        prev_records.append(tuple(r))
    prev_dir = config.SNAPSHOTS_DIR / prev_id
    prev_dir.mkdir(parents=True, exist_ok=True)
    (prev_dir / "manifest.json").write_text(json.dumps(
        {**manifest_data, "snapshot_id": prev_id, "created": prev_id[:10],
         "clinvar_date": prev_id[:10]}, indent=2))
    con.execute(f"DELETE FROM {load_mod.TABLE} WHERE snapshot_id = ?", [prev_id])
    con.executemany(f"""
        INSERT INTO {load_mod.TABLE} (
            snapshot_id, sample_id, chrom, pos, ref, alt, gene, consequence,
            clinvar_sig, clinvar_revstat, clinvar_alleleid, clinvar_disease,
            gnomad_af, sift, polyphen, zygosity
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, prev_records)
    con.close()

    # The earlier snapshot gets its own reports and pages, so the navigation tree's
    # "Snapshots" section has something to switch to.
    with demo_mod.notice(demo_mod.EARLIER_NOTICE):
        report_mod.full_report(prev_id)
        panel_mod.panel_report(prev_id)
    report_mod.full_report(snapshot_id)
    panel_mod.panel_report(snapshot_id)
    # History for both snapshots first, so the triage stability factor sees the trail.
    for sid in (prev_id, snapshot_id):
        provenance_mod.record(sid)
    with demo_mod.notice(demo_mod.DELTA_NOTICE):
        diff_mod.diff(prev_id, snapshot_id)
    snap_mod.mark_scan_complete(snapshot_id)

    # 5. HTML dashboards & index rendering
    log.info("6/6 Building the dashboards and navigation pages...")
    _build_dashboards(prev_id)
    _build_dashboards(snapshot_id)

    from pipeline import reportpaths
    html_dir = reportpaths.html_root() / snapshot_id
    father_dash = html_dir / "index_father.html"
    mother_dash = html_dir / "index_mother.html"
    daughter_dash = html_dir / "index_daughter.html"
    partner_dash = html_dir / "index_partner.html"
    grandson_dash = html_dir / "index_grandson.html"
    granddaughter_dash = html_dir / "index_granddaughter.html"

    print("\n" + "=" * 80)
    print("🎉 Genome Ledger Demo Complete! (3 Generations, 6 Family Members)")
    print("=" * 80)
    print(f"Data & reports created in: {demo_root}\n")
    print("Key Multi-Generational Findings:")
    print("  • Autosomal Dominant (AD) Transmission — Familial Hypercholesterolemia (APOB):")
    print("    - Gen 1: Father (At-Risk: APOB heterozygous)")
    print("    - Gen 2: Daughter (At-Risk, inherited from Father)")
    print("    - Gen 3: Grandson (At-Risk, inherited from Daughter) vs Granddaughter (Normal / Spared)")
    print("  • Autosomal Recessive (AR) Transmission & Carrier Matching (HFE Hemochromatosis):")
    print("    - Gen 1: Father & Mother both carriers (25% risk) → Daughter Affected (HOM)")
    print("    - Gen 2: Daughter (Affected) pairs with Partner (Carrier) → 50% child risk in panel_family.md!")
    print("    - Gen 3: Grandson Affected (HOM) and Granddaughter Carrier (HET)")
    print("  • Carrier Screening Peace of Mind (CFTR Cystic Fibrosis):")
    print("    - Daughter is a Carrier, but Partner is Normal (0/0) → 0% disease risk for children")
    print("  • HLA (read from the demo genomes): celiac DQ2.5 haplotype Father → Daughter → Grandson")
    print("  • What changed (DELTA vs an earlier demo snapshot): APOB reclassified VUS → Pathogenic,")
    print("    CFTR review status raised, PEX10 Likely benign → VUS — triaged by the real diff code")
    print("  • Single-SNP Traits & Wellness:")
    print("    - Lactase: Grandparents Persistent/Intolerant → Grandchildren segregation")
    print("    - Eye colour: Grandfather (Hazel) × Grandmother (Blue) → Blue & Hazel transmission\n")
    print("Illustrative sections (bundled sample results — NOT computed from the demo genomes;")
    print("each page says so): repeat expansions, pharmacogenomics, haplogroups, polygenic scores.\n")
    print(f"Markdown reports (the source of truth; point a retriever here): {reportpaths.md_root()}")
    print("Interactive HTML Dashboards:")
    print(f"  • Grandfather:   file://{father_dash}")
    print(f"  • Grandmother:   file://{mother_dash}")
    print(f"  • Daughter:      file://{daughter_dash}")
    print(f"  • Partner:       file://{partner_dash}")
    print(f"  • Grandson:      file://{grandson_dash}")
    print(f"  • Granddaughter: file://{granddaughter_dash}")
    print(f"  • Family inheritance tree: file://{html_dir / 'inheritance.html'}\n")
    print("To launch the local web UI and navigate the family tree:")
    print(f"  GENOMES_ROOT={demo_root} uv run python run.py serve")
    print("=" * 80 + "\n")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="run.py", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-v", "--verbose", action="store_true")
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("setup", help="install tools + download reference data")

    sp = sub.add_parser("align", help="FASTQ → marked-dup CRAM (fastp|bwa-mem|markdup)")
    sp.add_argument("sample")
    sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("call", help="CRAM → VCF (GATK HaplotypeCaller, scattered)")
    sp.add_argument("sample")
    sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("ingest", help="register a raw VCF as immutable source")
    sp.add_argument("vcf")
    sp.add_argument("--sample", required=True)
    sp.add_argument("--relation", required=True,
                    help="father/mother/self/son/daughter/…")
    sp.add_argument("--sex", required=True, help="male/female")

    sp = sub.add_parser("normalize", help="split multiallelics + left-align a sample")
    sp.add_argument("sample")
    sp.add_argument("--force", action="store_true")

    sub.add_parser("snapshot", help="create a dated, pinned DB snapshot")

    sp = sub.add_parser("annotate", help="VEP-annotate a sample against a snapshot")
    sp.add_argument("sample")
    sp.add_argument("snapshot")
    sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("load", help="load annotated output into DuckDB")
    sp.add_argument("sample")
    sp.add_argument("snapshot")

    sp = sub.add_parser("report", help="write FULL reports for a snapshot")
    sp.add_argument("snapshot")

    sp = sub.add_parser("acmg", help="ACMG/AMP-style classification of panel-gene variants")
    sp.add_argument("snapshot")

    sp = sub.add_parser("panel-report",
                        help="write monogenic disease-panel reports for a snapshot")
    sp.add_argument("snapshot")

    sp = sub.add_parser("pgx", help="pharmacogenomics report (PharmCAT) for a sample")
    sp.add_argument("sample")
    sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("ancestry", help="mtDNA + Y haplogroup report for a sample")
    sp.add_argument("sample")
    sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("traits", help="single-SNP traits/wellness report for a sample")
    sp.add_argument("sample")
    sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("hla", help="HLA risk (B27, celiac) via tag SNPs for a sample")
    sp.add_argument("sample")
    sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("hla-type", help="full classical HLA typing (arcasHLA, Docker) "
                        "for a sample")
    sp.add_argument("sample")
    sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("mirror-1kg", help="local mirror of the 1000G phased autosome VCFs "
                        "(~34.6GB; speeds up prs region queries)")
    sp.add_argument("chroms", nargs="*", help="specific chromosomes (default: all 1-22)")

    sp = sub.add_parser("prs", help="polygenic risk scores (PGS Catalog) for a sample")
    sp.add_argument("sample")
    sp.add_argument("--force", action="store_true")
    sp.add_argument("--max-variants", type=int, metavar="N",
                    help="skip scores with more than N variants. Half the trait panel is "
                         "genome-wide, and each such score needs its own cached 1000G "
                         "reference-panel run; --max-variants 50000 scores the cheap half in "
                         "minutes and leaves the rest to a later unrestricted run.")
    sp.add_argument("--s1-gate", action=argparse.BooleanOptionalAction, default=None,
                    help="withhold severe non-actionable traits from anyone who has not "
                         "opted in (S1 incidental-findings consent filter). On by default; "
                         "--no-s1-gate discloses everything to everyone (or PRS_S1_GATE=0).")

    sp = sub.add_parser("repeats", help="repeat-expansion report (ExpansionHunter) for a sample")
    sp.add_argument("sample")
    sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("sv", help="structural-variant / CNV report (Manta, Docker) for a sample")
    sp.add_argument("sample")
    sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("callable", help="callability / coverage report (mosdepth) for a sample")
    sp.add_argument("sample")
    sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("mito",
                        help="mitochondrial heteroplasmy report (Mutect2 mito-mode) for a sample")
    sp.add_argument("sample")
    sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("consent",
                        help="grant/revoke a family member's opt-in for an incidental-findings "
                             "category (S1); omit category/decision to list current grants")
    sp.add_argument("sample", nargs="?")
    sp.add_argument("category", nargs="?",
                    choices=[consent_mod.CATEGORY_ADULT_ONSET_UNTREATABLE])
    sp.add_argument("decision", nargs="?", choices=["yes", "no"])

    sp = sub.add_parser("family-migrate",
                        help="rename relationship-style sample ids (father/mother) to given names")
    sp.add_argument("--apply", action="store_true",
                    help="execute the renames (default is dry-run)")
    sp.add_argument("--reheader-vcfs", action="store_true",
                    help="also rewrite the VCF sample columns (cheap; CRAM @RG left as-is)")
    sp.add_argument("--map", dest="id_map", metavar="OLD=NEW[,OLD=NEW...]",
                    help="rename these ids instead of the built-in father/mother map, e.g. "
                         "--map SAMPLE1=SAMPLE1_vendor,SAMPLE2=SAMPLE2_vendor")

    sp = sub.add_parser("denovo", help="trio de-novo variants (child vs both parents, "
                                      "confirmed off the parents' CRAMs)")
    sp.add_argument("child")
    sp.add_argument("--force", action="store_true")
    sp.add_argument("--no-confirm", action="store_true",
                    help="skip the parental force-call pass (candidates only — "
                         "FP-dominated, for inspecting the filter, not for reading counts)")

    sp = sub.add_parser("segregation",
                        help="trio segregation: resolve compound hets by which parent "
                             "transmitted each variant (confirmed off the parents' CRAMs)")
    sp.add_argument("child")
    sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("phase", help="read-backed phasing cis/trans (WhatsHap) for a sample")
    sp.add_argument("sample")
    sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("force-call",
                        help="all-sites genotyping at curated positions (GATK) for a sample")
    sp.add_argument("sample")
    sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("interpret",
                        help="non-coding/regulatory effect (Enformer on the GPU box) for a snapshot")
    sp.add_argument("snapshot")
    sp.add_argument("--force", action="store_true")
    sp.add_argument("--calibrate", action="store_true",
                    help="build the benign-null calibration (percentile bands) instead of a report")

    sp = sub.add_parser("provenance",
                        help="record provenance/confidence + variant history for a snapshot")
    sp.add_argument("snapshot", nargs="?", help="omit with --backfill")
    sp.add_argument("--backfill", action="store_true",
                    help="(re)build variant_history over every snapshot already in DuckDB")

    sp = sub.add_parser("diff", help="write DELTA report between two snapshots")
    sp.add_argument("prev")
    sp.add_argument("curr")

    sp = sub.add_parser("validate",
                        help="replay an old snapshot's triage vs a later ClinVar snapshot")
    sp.add_argument("--replay", required=True, help="snapshot whose deltas to judge")
    sp.add_argument("--truth", required=True, help="later snapshot providing ground-truth ClinVar")

    sp = sub.add_parser("scan", help="periodic chain (snapshot→…→diff); incremental by default")
    sp.add_argument("--full", action="store_true",
                    help="force full re-annotation of every sample (the reproducibility anchor)")
    sub.add_parser("check-releases", help="scan only if upstream DB changed (launchd)")

    sp = sub.add_parser("query", help="run a parameterized SQL file against DuckDB")
    sp.add_argument("sqlfile")
    sp.add_argument("--param", action="append", help="key=value (repeatable)")

    sp = sub.add_parser("share", help="copy a person's VCF or CRAM for handing to someone "
                        "else: header rebuilt without command lines, paths or internal ids")
    sp.add_argument("sample")
    sp.add_argument("--cram", action="store_true", help="the CRAM instead of the VCF")
    sp.add_argument("--name", help="sample name to write into the file (default: the "
                    "person's id from family.tsv)")

    sp = sub.add_parser("render",
                        help="(re)generate HTML + dashboards from existing markdown")
    sp.add_argument("--snapshot", help="rebuild this dated report dir instead of the latest")

    sp = sub.add_parser("migrate-reports", help="one-time: move reports from the old layout "
                        "into reports/md/ (dry-run unless --apply)")
    sp.add_argument("--apply", action="store_true")

    sp = sub.add_parser("serve", help="local web server for the report dashboards")
    sp.add_argument("--host", default="127.0.0.1")
    sp.add_argument("--port", type=int, default=8765)

    sp = sub.add_parser("demo", help="run zero-data end-to-end demonstration using synthetic trio fixtures")
    sp.add_argument("--out-dir", default="demo_output",
                    help="directory for demo data and reports (default: demo_output)")
    sp.add_argument("--clean", action="store_true",
                    help="clean output directory before running demo")
    return p


def main(argv: list[str]) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(args.verbose)
    match args.cmd:
        case "setup":
            setup_mod.setup(); return 0
        case "align":
            call_mod.align(args.sample, args.force); return 0
        case "call":
            call_mod.call(args.sample, args.force); return 0
        case "ingest":
            ingest_mod.ingest(args.vcf, args.sample, args.relation, args.sex); return 0
        case "normalize":
            norm_mod.normalize(args.sample, args.force); return 0
        case "snapshot":
            print(snap_mod.create_snapshot()); return 0
        case "annotate":
            annotate_mod.annotate(args.sample, args.snapshot, args.force); return 0
        case "load":
            load_mod.load(args.sample, args.snapshot); return 0
        case "report":
            report_mod.full_report(args.snapshot); return 0
        case "panel-report":
            panel_mod.panel_report(args.snapshot); return 0
        case "acmg":
            acmg_mod.acmg_report(args.snapshot); return 0
        case "pgx":
            pgx_mod.pgx_report(args.sample, args.force); return 0
        case "ancestry":
            ancestry_mod.ancestry_report(args.sample, args.force); return 0
        case "traits":
            traits_mod.traits_report(args.sample, args.force); return 0
        case "hla":
            hla_mod.hla_report(args.sample, args.force); return 0
        case "hla-type":
            hla_typing_mod.hla_typing_report(args.sample, args.force); return 0
        case "mirror-1kg":
            setup_mod.mirror_1000g(args.chroms or None); return 0
        case "prs":
            if args.max_variants is not None:
                config.PRS_MAX_VARIANTS = args.max_variants
            if args.s1_gate is not None:
                config.PRS_S1_GATE = args.s1_gate
            prs_mod.prs_report(args.sample, args.force); return 0
        case "repeats":
            repeats_mod.repeats_report(args.sample, args.force); return 0
        case "sv":
            sv_mod.sv_report(args.sample, args.force); return 0
        case "callable":
            callability_mod.callability_report(args.sample, args.force); return 0
        case "mito":
            mito_mod.mito_report(args.sample, args.force); return 0
        case "consent":
            if args.sample and args.category and args.decision:
                consent_mod.set_consent(args.sample, args.category, args.decision == "yes")
                if args.decision != "yes":
                    _refresh_person_pages(args.sample)
            else:
                rows = read_tsv(config.CONSENT_FILE)
                if args.sample:
                    rows = [r for r in rows if r["sample_id"] == args.sample]
                if not rows:
                    print("no consent grants recorded" + (f" for {args.sample}" if args.sample else ""))
                for r in rows:
                    print(f"{r['sample_id']}\t{r['category']}\t{r['consent']}")
            return 0
        case "family-migrate":
            # A rename is only safe while no *other* artifact already carries the target
            # id, so the mapping is explicit rather than inferred.
            mapping = None
            if args.id_map:
                mapping = {}
                for pair in args.id_map.split(","):
                    old_id, _, new_id = pair.partition("=")
                    old_id, new_id = old_id.strip(), new_id.strip()
                    if not old_id or not new_id or old_id == new_id:
                        raise SystemExit(f"--map: bad entry {pair!r}; expected OLD=NEW")
                    mapping[old_id] = new_id
            family_migrate_mod.migrate(apply=args.apply,
                                       reheader_vcfs=args.reheader_vcfs,
                                       mapping=mapping); return 0
        case "denovo":
            denovo_mod.denovo(args.child, args.force,
                              confirm_parents=not args.no_confirm); return 0
        case "segregation":
            segregation_mod.segregation_report(args.child, args.force); return 0
        case "phase":
            phasing_mod.phase_report(args.sample, args.force); return 0
        case "force-call":
            forcecall_mod.force_call(args.sample, args.force); return 0
        case "interpret":
            if args.calibrate:
                interpret_mod.calibrate(args.snapshot, args.force)
            else:
                interpret_mod.interpret_report(args.snapshot, args.force)
            return 0
        case "provenance":
            if args.backfill:
                provenance_mod.backfill_history()
            if args.snapshot:
                provenance_mod.record(args.snapshot)
            elif not args.backfill:
                raise SystemExit("provenance: give a snapshot id or --backfill")
            return 0
        case "diff":
            diff_mod.diff(args.prev, args.curr); return 0
        case "validate":
            validate_mod.validate(args.replay, args.truth); return 0
        case "scan":
            return cmd_scan(args)
        case "check-releases":
            return cmd_check_releases(args)
        case "query":
            return cmd_query(args)
        case "share":
            share_mod.share(args.sample, "cram" if args.cram else "vcf", args.name); return 0
        case "render":
            return cmd_render(args)
        case "migrate-reports":
            from pipeline import migrate_reports
            return migrate_reports.migrate(args.apply)
        case "serve":
            return cmd_serve(args)
        case "demo":
            return cmd_demo(args)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
