"""Stage: non-coding / regulatory variant-effect scoring via Enformer on the GPU box — I3.

The interpretation-side counterpart to the GPU front-end in docs/hardware-and-scope.md:
where AlphaMissense (missense) and SpliceAI (splicing) are joins over public precomputed
tables, *non-coding regulatory* effect isn't precomputed genome-wide anywhere usable — so we
run DeepMind's Enformer ourselves on a GPU box. It predicts ~5,313 regulatory tracks
(expression, accessibility, TF/histone marks) from a 196,608 bp context; the per-variant
effect is the change in those tracks between the reference and alternate allele.

Division of labour (keeps the box stateless — a clean "hand it a job, get scores back" node):
  • Local (here): select the family's carried non-coding variants from DuckDB, extract the
                 196,608 bp ref + alt windows from the reference FASTA (samtools faidx), and
                 ship a compact gzipped job over SSH.
  • Box (worker): the baked `enformer-worker` image runs pipeline/interpret_worker.py on the
                 GPU and returns a per-variant effect table.
Results merge into a genome-level cache (config.ENFORMER_SCORES_FILE) joined by
(chrom,pos,ref,alt) — see pipeline/enformer.py — so re-runs only score *new* variants, and
the FULL report enriches itself from the cache (degrading gracefully when it's absent).
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import subprocess
from pathlib import Path

import duckdb

from . import reportpaths
from . import config
from . import enformer as enformer_mod
from .util import log, require_tools, run, write_report

SEQ_LEN = config.ENFORMER_SEQ_LEN
HALF = SEQ_LEN // 2
_SCORE_COLS = ["chrom", "pos", "ref", "alt", "delta_max", "l2_center", "top_track", "top_desc"]


def _select_variants(snapshot: str, already: set) -> list[dict]:
    """Carried non-coding variants for the snapshot, deduped to genome level (one row per
    distinct variant, with its carriers), rare-first, capped — and not already scored."""
    cons_filter = " OR ".join(["consequence LIKE ?"] * len(config.INTERPRET_CONSEQUENCES))
    params = [snapshot, config.INTERPRET_MAX_AF]
    params += [f"%{c}%" for c in config.INTERPRET_CONSEQUENCES]
    con = duckdb.connect(str(config.DUCKDB_FILE), read_only=True)
    rows = con.execute(f"""
        SELECT chrom, pos, ref, alt,
               any_value(gene)        AS gene,
               any_value(consequence) AS consequence,
               min(gnomad_af)         AS af,
               string_agg(DISTINCT sample_id, ',') AS carriers
        FROM variants
        WHERE snapshot_id = ?
          AND zygosity IN ('HET', 'HOM')
          AND (gnomad_af IS NULL OR gnomad_af <= ?)
          AND ({cons_filter})
        GROUP BY chrom, pos, ref, alt
        ORDER BY (min(gnomad_af) IS NULL) DESC, af ASC, chrom, pos, ref, alt
    """, params).fetchall()
    cols = [d[0] for d in con.description]
    con.close()
    out = []
    for r in rows:
        d = dict(zip(cols, r))
        if (str(d["chrom"]), int(d["pos"]), d["ref"], d["alt"]) in already:
            continue
        out.append(d)
        if len(out) >= config.INTERPRET_MAX_VARIANTS:
            break
    return out


def _faidx(chrom: str, start1: int, end1: int) -> str:
    """Reference bases for chrom:start-end (1-based, inclusive) as an uppercase string, or ""
    if the region can't be fetched.

    `util.run` uses check=True, and nothing here caught CalledProcessError — so a single
    variant on a contig absent from REF_FASTA (the pipeline builds against the *primary
    assembly*, while a vendor callset or VEP cache can carry patch/alt contigs) killed the
    whole interpret/calibrate stage, hours of GPU work in. `_windows` already treats an
    unusable window as "skip this variant", which is the right outcome here too."""
    proc = subprocess.run(
        ["samtools", "faidx", str(config.REF_FASTA), f"{chrom}:{start1}-{end1}"],
        capture_output=True, text=True)
    if proc.returncode != 0:
        log.warning("faidx could not fetch %s:%d-%d (%s) — variant skipped, not scored.",
                    chrom, start1, end1,
                    " | ".join(proc.stderr.strip().splitlines()[-1:]) or "no stderr")
        return ""
    return "".join(l.strip() for l in proc.stdout.splitlines()
                   if not l.startswith(">")).upper()


def _windows(v: dict) -> tuple[str, str] | None:
    """Build the (ref_window, alt_window) pair of exactly SEQ_LEN bases centred on the
    variant. SNVs are first-class; small indels are best-effort (re-anchored to fixed
    length). Returns None if the variant can't be placed (edge of contig, ref mismatch,
    too large an indel). The variant's REF allele must match the reference — otherwise the
    window is wrong and we skip rather than emit a bogus score."""
    chrom, pos = str(v["chrom"]), int(v["pos"])
    # VEP encodes a pure deletion's alt (and a pure insertion's ref) as "-" (empty allele);
    # treat those as the empty string so the splice removes/inserts bases cleanly rather than
    # injecting a literal "-" into the sequence.
    ref = "" if v["ref"] == "-" else v["ref"]
    alt = "" if v["alt"] == "-" else v["alt"]
    # Window so the variant's first base sits at the centre index HALF.
    start = pos - HALF
    end = start + SEQ_LEN - 1
    if start < 1:
        return None                              # near contig start; pad logic deferred
    win = _faidx(chrom, start, end)
    if len(win) != SEQ_LEN:                      # ran off the contig end
        return None
    c = HALF                                     # index of `pos` within the window
    if win[c:c + len(ref)] != ref:               # reference disagreement → don't trust it
        return None
    alt_win = win[:c] + alt + win[c + len(ref):]
    # Re-trim/pad to the fixed length around the centre (indels change length).
    if len(alt_win) > SEQ_LEN:
        alt_win = alt_win[:SEQ_LEN]
    elif len(alt_win) < SEQ_LEN:
        pad_end = _faidx(chrom, end + 1, end + 1 + (SEQ_LEN - len(alt_win)))
        alt_win = (alt_win + pad_end)[:SEQ_LEN]
    if len(alt_win) != SEQ_LEN or len(win) != SEQ_LEN:
        return None
    return win, alt_win


def _vid(v: dict) -> str:
    return f"{v['chrom']}_{v['pos']}_{v['ref']}_{v['alt']}"


def _write_job(variants: list[dict], job_path: Path) -> dict[str, dict]:
    """Build windows and write the gzipped job TSV. Returns the id→variant map of those
    actually placed (some are skipped for edges/ref-mismatch/oversize indels)."""
    placed: dict[str, dict] = {}
    skipped = 0
    with gzip.open(job_path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        for v in variants:
            if len(v["ref"]) > 200 or len(v["alt"]) > 200:   # SV-scale → not Enformer's job
                skipped += 1
                continue
            win = _windows(v)
            if win is None:
                skipped += 1
                continue
            vid = _vid(v)
            w.writerow([vid, win[0], win[1]])
            placed[vid] = v
    log.info("interpret: %d variants placed into job, %d skipped (edge/mismatch/indel).",
             len(placed), skipped)
    return placed


def _run_remote(job_path: Path, out_path: Path) -> None:
    """scp the job to the box, run the baked enformer-worker on the GPU, scp results back."""
    require_tools("ssh", "scp")
    host = config.INTERPRET_SSH_HOST
    rdir = config.INTERPRET_REMOTE_DIR
    cache = config.INTERPRET_CACHE_DIR
    rjob = f"{rdir}/{job_path.name}"
    rout = f"{rdir}/{out_path.name}"
    model = config.INTERPRET_MODEL
    targets = f"{cache}/targets_{model}.txt"
    run(["ssh", host, f"mkdir -p {rdir} {cache}/hf"])
    run(["scp", "-q", str(job_path), f"{host}:{rjob}"])
    # Best-effort: ensure the model's track-description file is in the cache for human labels.
    run(["ssh", host, f"test -f {targets} || curl -fsSL -o {targets} "
         f"{config.ENFORMER_TARGETS_URL} || true"])
    docker = (
        f"docker run --rm --gpus all "
        f"-v {cache}:/cache -v {rdir}:/jobs "
        f"{config.ENFORMER_IMAGE} --model {model} --targets /cache/targets_{model}.txt "
        f"--job /jobs/{job_path.name} --out /jobs/{out_path.name} "
        f"--batch-size {config.ENFORMER_BATCH_SIZE}"
    )
    log.info("interpret: running %s on %s …", model, host)
    with open(config.LOGS_DIR / "interpret_worker.log", "wb") as logf:  # no fd leak
        run(["ssh", host, docker], stdout=logf, stderr=subprocess.STDOUT)
    run(["scp", "-q", f"{host}:{rout}", str(out_path)])


def _merge_scores(out_path: Path, placed: dict[str, dict]) -> list[dict]:
    """Fold the worker's results into the genome-level cache (keyed by chrom,pos,ref,alt);
    return the freshly scored rows (joined back to gene/consequence/carriers for the report)."""
    config.ENFORMER_SCORES_FILE.parent.mkdir(parents=True, exist_ok=True)
    cache: dict[tuple, list] = {}
    if config.ENFORMER_SCORES_FILE.exists():
        with open(config.ENFORMER_SCORES_FILE) as fh:
            next(fh, None)
            # csv.reader matches the csv.writer below — a plain split("\t") misreads a
            # quoted field (tab/quote in a track description) and each cache round-trip
            # would compound the corruption.
            for f in csv.reader(fh, delimiter="\t"):
                if len(f) >= 7:
                    cache[(f[0], f[1], f[2], f[3])] = f

    fresh = []
    with open(out_path) as fh:
        next(fh, None)
        for f in csv.reader(fh, delimiter="\t"):
            if len(f) < 6:
                continue
            vid, delta_max, l2c, track, _top_delta, *rest = f
            desc = rest[0] if rest else ""
            v = placed.get(vid)
            if not v:
                continue
            key = (str(v["chrom"]), str(v["pos"]), v["ref"], v["alt"])
            cache[key] = [key[0], key[1], key[2], key[3], delta_max, l2c, track, desc]
            fresh.append({**v, "delta_max": float(delta_max), "l2_center": float(l2c),
                          "top_track": int(track), "top_desc": desc})

    # Temp + rename. Opening the real file "w" truncated the accumulated scores in place:
    # a crash (or an OOM kill on the GPU box's result merge) mid-write destroyed every score
    # ever computed, and enformer._load() reading concurrently would parse a half-written
    # file into a partial lookup table — silently, since its per-row guards skip short lines.
    dest = config.ENFORMER_SCORES_FILE
    tmp = dest.with_name(f"{dest.name}.{os.getpid()}.part")
    try:
        with open(tmp, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(_SCORE_COLS)
            for row in sorted(cache.values(), key=lambda r: (r[0], int(r[1]))):
                w.writerow(row)
        os.replace(tmp, dest)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    return fresh


def _report(snapshot: str, fresh: list[dict], n_selected: int, summary: Path) -> Path:
    fresh.sort(key=lambda d: d["delta_max"], reverse=True)
    calibrated = enformer_mod.calibration_available()
    n_hit = sum(1 for d in fresh
                if enformer_mod.band(d["delta_max"], d["ref"], d["alt"]) in ("moderate", "high"))
    hit_basis = (f"≥ {config.ENFORMER_PCTL_MODERATE:g}th pctl vs benign null" if calibrated
                 else f"Δ ≥ {config.ENFORMER_DELTA_MODERATE:g}, uncalibrated")
    lines = [
        f"# Non-coding regulatory effect ({config.INTERPRET_MODEL}) — {snapshot}",
        "",
        f"- **Model:** {config.INTERPRET_MODEL} (`{config.INTERPRET_MODEL_HF}`), "
        f"{config.ENFORMER_SEQ_LEN:,} bp context, predicted regulatory tracks; run on the "
        "GPU box.",
        f"- **Scored:** {len(fresh)} of {n_selected} selected non-coding/regulatory variants "
        f"(carried, rare ≤ {config.INTERPRET_MAX_AF:g} AF; the rest skipped at edge/large-indel).",
        f"- **Calibration:** {'percentile vs a benign/common-variant null (SNV & indel strata)' if calibrated else '**none yet** — run `interpret calibrate ' + snapshot + '`; bands below are raw-Δ heuristics'}.",
        f"- **Predicted regulatory hits ({hit_basis}):** {n_hit}.",
        "",
        "## Top predicted regulatory perturbations",
        "",
    ]
    if fresh:
        lines += [
            "| Variant (GRCh38) | Gene | Consequence | Δmax | Pctl | Band | Top track | Carriers |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for d in fresh[:40]:
            assay = d["top_desc"] or f"track {d['top_track']}"
            pct = enformer_mod.percentile(d["delta_max"], d["ref"], d["alt"])
            pct_s = f"{pct:.1f}" if pct is not None else "—"
            lines.append(
                f"| {d['chrom']}:{d['pos']} {d['ref']}>{d['alt']} | {d.get('gene') or '—'} | "
                f"{(d.get('consequence') or '—').replace('_',' ')} | {d['delta_max']:.2f} | "
                f"{pct_s} | {enformer_mod.band(d['delta_max'], d['ref'], d['alt'])} | "
                f"{assay} | {d.get('carriers','')} |")
    else:
        lines.append("_No variants scored this run (all already cached, or none selected)._")
    lines += [
        "",
        "## Caveats",
        "- **Δmax is a model-predicted track change, not a calibrated pathogenicity.** It is *not* on "
        "a [0,1] scale like SpliceAI. When calibrated, **Pctl** is the percentile of Δmax against "
        "common (≈benign) variants of the same kind (SNV/indel) — \"more disruptive than N% of "
        "benign variants of its type\" — and the bands key off it (≥"
        f"{config.ENFORMER_PCTL_HIGH:g}th high, ≥{config.ENFORMER_PCTL_MODERATE:g}th moderate).",
        "- Predicts *regulatory* effect (expression, accessibility, TF/histone marks), not "
        "coding or splice impact — pair with AlphaMissense (missense) and SpliceAI (splice).",
        "- Effect summarised over the central bins of the 114 kb receptive field; very distal "
        "enhancer effects and edge-of-contig variants are out of scope here.",
        "- Research-grade, not a clinical diagnosis.",
        "",
    ]
    write_report(summary, lines)
    log.info("interpret: report %s (%d scored, %d hits, calibrated=%s).",
             summary, len(fresh), n_hit, calibrated)
    return summary


def _select_common(snapshot: str) -> list[dict]:
    """Common (≈benign) non-coding variants for the null, in two strata (SNV / indel), each
    capped at INTERPRET_CALIBRATION_N and spread across the genome (ordered by a stable hash)."""
    cons_filter = " OR ".join(["consequence LIKE ?"] * len(config.INTERPRET_CONSEQUENCES))
    cons_params = [f"%{c}%" for c in config.INTERPRET_CONSEQUENCES]
    snv_pred = ("length(ref)=1 AND length(alt)=1 AND ref IN ('A','C','G','T') "
                "AND alt IN ('A','C','G','T')")
    con = duckdb.connect(str(config.DUCKDB_FILE), read_only=True)

    def q(extra: str) -> list[dict]:
        rows = con.execute(f"""
            SELECT chrom, pos, ref, alt,
                   any_value(gene) AS gene, any_value(consequence) AS consequence
            FROM variants
            WHERE snapshot_id = ? AND gnomad_af >= ? AND ({cons_filter}) AND {extra}
            GROUP BY chrom, pos, ref, alt
            ORDER BY hash(chrom, pos, ref, alt)
            LIMIT ?
        """, [snapshot, config.INTERPRET_COMMON_AF, *cons_params,
              config.INTERPRET_CALIBRATION_N]).fetchall()
        cols = [d[0] for d in con.description]
        return [dict(zip(cols, r)) for r in rows]

    out = q(snv_pred) + q(f"NOT ({snv_pred})")
    con.close()
    return out


def calibrate(snapshot: str, force: bool = False) -> Path:
    """Build the benign null: score a sample of common non-coding variants and save the sorted
    Δmax distribution per stratum (SNV / indel). Subsequent reports express each variant's Δ as
    a percentile against this null. Genome/model-level — run once, reused for the whole family."""
    require_tools("samtools")
    if not config.REF_FASTA.exists():
        raise SystemExit(f"{config.REF_FASTA} missing. Run `setup` first.")
    if config.ENFORMER_CALIBRATION_FILE.exists() and not force:
        log.info("interpret: calibration exists (%s) — skipping (use --force).",
                 config.ENFORMER_CALIBRATION_FILE)
        return config.ENFORMER_CALIBRATION_FILE

    variants = _select_common(snapshot)
    if not variants:
        raise SystemExit(f"No common (AF≥{config.INTERPRET_COMMON_AF:g}) non-coding variants "
                         f"in {snapshot} to build a null from.")
    out_dir = config.INTERPRET_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    job = out_dir / f"caljob_{snapshot}.tsv.gz"
    out = out_dir / f"calscores_{snapshot}.tsv"
    placed = _write_job(variants, job)
    if not placed:
        raise SystemExit("interpret calibrate: no variants could be placed into a job.")
    _run_remote(job, out)
    enformer_mod._table = None
    fresh = _merge_scores(out, placed)            # also folds the null scores into the cache

    snv = sorted(d["delta_max"] for d in fresh if enformer_mod.is_snv(d["ref"], d["alt"]))
    indel = sorted(d["delta_max"] for d in fresh if not enformer_mod.is_snv(d["ref"], d["alt"]))
    data = {"model": config.INTERPRET_MODEL_HF,
            "af_threshold": config.INTERPRET_COMMON_AF,
            "n_snv": len(snv), "n_indel": len(indel), "snv": snv, "indel": indel}
    config.ENFORMER_CALIBRATION_FILE.write_text(json.dumps(data))
    enformer_mod._calibration = None              # invalidate the lazy cache
    log.info("interpret calibrate: null built — %d SNV + %d indel (AF≥%.2g) → %s",
             len(snv), len(indel), config.INTERPRET_COMMON_AF,
             config.ENFORMER_CALIBRATION_FILE.name)
    return config.ENFORMER_CALIBRATION_FILE


def interpret_report(snapshot: str, force: bool = False) -> Path:
    require_tools("samtools")
    if not config.REF_FASTA.exists():
        raise SystemExit(f"{config.REF_FASTA} missing. Run `setup` first.")

    out_dir = config.INTERPRET_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = reportpaths.family_report(f"interpret_{config.INTERPRET_MODEL}", snapshot)
    if summary.exists() and not force:
        log.info("interpret already reported (%s) — skipping (use --force).", summary)
        return summary

    # The full candidate set (rare carried non-coding, capped), independent of the cache, so the
    # report is stable and complete: we score only the *uncached* ones, then render *all* selected
    # variants from the cache. This also lets a re-run pick up a new calibration with no re-scoring.
    variants = _select_variants(snapshot, already=set())
    if not variants:
        log.info("interpret: no non-coding variants selected for %s.", snapshot)
        return _report(snapshot, [], 0, summary)

    cache = enformer_mod._load()
    to_score = [v for v in variants
                if (str(v["chrom"]), int(v["pos"]), v["ref"], v["alt"]) not in cache]
    if to_score:
        job = out_dir / f"job_{snapshot}.tsv.gz"
        out = out_dir / f"scores_{snapshot}.tsv"
        placed = _write_job(to_score, job)
        if placed:
            _run_remote(job, out)
            enformer_mod._table = None         # invalidate so the merge re-reads the cache
            _merge_scores(out, placed)
            enformer_mod._table = None
            cache = enformer_mod._load()

    # Join every selected variant to its cached score (some may stay unscored — edge/indel skips).
    scored = []
    for v in variants:
        en = cache.get((str(v["chrom"]), int(v["pos"]), v["ref"], v["alt"]))
        if en:
            scored.append({**v, "delta_max": en[0], "l2_center": en[1],
                           "top_track": en[2], "top_desc": en[3]})
    return _report(snapshot, scored, len(variants), summary)
