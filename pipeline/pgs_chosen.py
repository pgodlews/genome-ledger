"""Polygenic scores the family chooses from the whole PGS Catalog.

The curated panel (panels/complex_traits.tsv, prs.py) is one score per trait, picked by fixed
criteria. This module lets anyone pick any score from the Catalog instead, and keeps what it
produces apart from the panel's:

    working data   config.PRS_CHOSEN_DIR/<callset>/<PGS id>.json   one result each
    report         reports/md/genome/<person>/chosen/prs_<callset>.md   one per callset,
                   rebuilt from whatever results exist; the panel report is never touched

**Consent comes first, and it gates the data, not the display.** A score is downloaded and
computed for a person only once that person has consented to every trait the score is mapped
to (consent category `pgs:<trait id>`, in family/consent.tsv like every other grant). Nothing
is computed and then hidden: a result that was never consented to never exists on disk, so it
cannot reach a backup, a dashboard or a retriever later. Withdrawing consent deletes the
results that depended on it and rebuilds the report.

Consent is per trait, results are per score: agreeing to "age-related macular degeneration"
covers a newer AMD score too, while every result still names the exact score that produced it.
There is no age check. Guidelines advise against predictive testing of children for adult-onset
conditions with no treatment; the matrix says so every time, and the parent decides.

**Air-gapped machines.** Everything a score needs that does not depend on the family (the
scoring file, the 1000 Genomes reference distribution, the reference alleles at its sites) is
built by `fetch` on a connected machine with no genomes on it, written as a *bundle* with a
checksum manifest, carried across, and installed by `import_dir`. A bundle is the same set of
files prs.py caches, so scoring offline is the ordinary engine reading its ordinary cache.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import shutil
import urllib.request
import zlib
from pathlib import Path

from . import config, consent, pgs_catalog, prs, reportpaths
from .util import file_token, log, require_tools, sha256, today, write_report

CATEGORY_PREFIX = consent.PGS_PREFIX
# Bump when how a chosen result is computed changes: every result stamped with an older value
# is then stale and `pgs score --stale` recomputes it.
ENGINE = "1"
BUNDLE_FORMAT = 1
BUNDLE_FILE = "bundle.json"

CHILD_WARNING = (
    "If any of these people is a child: medical genetics guidelines (ACMG, ESHG) advise "
    "against predictive testing of children for adult-onset conditions that have no treatment "
    "or prevention, so that the child can decide as an adult. That is advice, not a rule here: "
    "the decision is the parent's.")
RELATIVES_NOTE = (
    "A polygenic score says a little about blood relatives too: a child's expected score is the "
    "average of the parents'. Much less than a single-gene result does, but not nothing.")


def category(efo: str) -> str:
    return CATEGORY_PREFIX + efo


# --- what a score is, from its scoring file or the catalog --------------------------------
def scoring_path(pid: str) -> Path:
    return config.PGS_DIR / f"{pid}_hmPOS_GRCh38.txt.gz"


def kg_dir(pid: str) -> Path:
    return config.PGS_DIR / f"{pid}_1kg"


def parse_header(lines) -> dict[str, str]:
    """`#key=value` lines at the top of a harmonized scoring file."""
    out = {}
    for line in lines:
        if not line.startswith("#"):
            break
        key, eq, value = line.lstrip("#").partition("=")
        if eq:
            out[key.strip()] = value.strip()
    return out


def file_header(path: Path) -> dict[str, str]:
    with gzip.open(path, "rt") as fh:
        head = []
        for line in fh:
            if not line.startswith("#"):
                break
            head.append(line.rstrip("\n"))
    return parse_header(head)


def _split_ids(text: str) -> list[str]:
    return [x for x in (t.strip() for t in re.split(r"[|,]", text or "")) if x]


def describe(pid: str) -> dict:
    """{pgs_id, traits, efo_ids, n_variants, hm_date} from the local scoring file when it is
    here, else from the catalog snapshot. Works before anything is downloaded: consent has to
    be decided before the download, so it cannot depend on it."""
    path = scoring_path(pid)
    if path.exists():
        h = file_header(path)
        efos = _split_ids(h.get("trait_efo", ""))
        labels = [x.strip() for x in (h.get("trait_mapped") or "").split("|") if x.strip()]
        if len(labels) != len(efos):
            labels = [x.strip() for x in re.split(r"[|,]", h.get("trait_mapped") or "")
                      if x.strip()]
        try:
            n = int(h.get("variants_number") or 0)
        except ValueError:
            n = 0
        if efos:
            return {"pgs_id": pid, "efo_ids": efos,
                    "traits": labels if len(labels) == len(efos) else efos,
                    "n_variants": n, "hm_date": h.get("HmPOS_date", ""),
                    "reported_trait": h.get("trait_reported", "")}
    try:
        s = pgs_catalog.load().summary(pid)
    except SystemExit:
        s = {}
    if not s:
        raise SystemExit(f"{pid}: not in the catalog snapshot and no scoring file here. Run "
                         "`run.py pgs sync`, or import a bundle that carries it.")
    return {"pgs_id": pid, "efo_ids": s["efo_ids"], "traits": s["traits"],
            "n_variants": s["n_variants"], "hm_date": "", "reported_trait": s["reported_trait"]}


def engine_problem(pid: str) -> str | None:
    """Why prs.py's additive engine cannot score this file, or None. Checked on the file
    itself before any reference data is built, and again before scoring: catalog metadata
    can lag the file, and per-dosage or haplotype weight columns appear only in the file."""
    path = scoring_path(pid)
    if not path.exists():
        return None
    h = file_header(path)
    why = pgs_catalog.unsupported(h.get("weight_type"))
    if why:
        return why
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                cols = set(line.rstrip("\n").split("\t"))
                break
        else:
            return "no column header"
    special = sorted(cols & {"dosage_0_weight", "dosage_1_weight", "dosage_2_weight",
                             "is_haplotype", "is_diplotype", "is_interaction"})
    if special:
        return f"uses {', '.join(special)}, which this engine does not model"
    return None


def trait_label(info: dict) -> str:
    return "; ".join(info["traits"]) or info.get("reported_trait") or info["pgs_id"]


def sensitive_traits() -> set[str]:
    """Trait ids the curated panel gates as severe and non-actionable (its S1 list). Used
    only to mark a column in the matrix; a chosen score always needs explicit consent."""
    out = set()
    try:
        import csv
        with config.COMPLEX_PANEL.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if (row.get("s1_gate") or "").strip() == "yes" and row.get("efo_id"):
                    out.add(row["efo_id"].strip())
    except OSError:
        pass
    return out


# --- people and consent -------------------------------------------------------------------
def people() -> list[dict]:
    """[{id, name, callsets}] for every person with at least one ingested callset."""
    from . import pedigree
    fam = pedigree.load(sequenced=pedigree.sequenced_ids())
    seen: dict[str, dict] = {}
    for cs in sorted(pedigree.sequenced_ids()):
        person = fam.person_of(cs)
        p = fam.people.get(person)
        seen.setdefault(person, {"id": person, "name": p.display_name if p else person,
                                 "callsets": fam.callsets(person)})
    return sorted(seen.values(), key=lambda p: p["id"])


def resolve_people(names: list[str] | None) -> list[dict]:
    everyone = people()
    if not names:
        return everyone
    by_key = {}
    for p in everyone:
        by_key[p["id"].lower()] = p
        by_key[p["name"].lower()] = p
        for cs in p["callsets"]:
            by_key[cs.lower()] = p
    out, unknown = [], []
    for n in names:
        p = by_key.get(n.strip().lower())
        (out.append(p) if p else unknown.append(n))
    if unknown:
        raise SystemExit(f"unknown person(s): {', '.join(unknown)} (known: "
                         f"{', '.join(p['id'] for p in everyone)})")
    return list({p["id"]: p for p in out}.values())


def may_score(person: str, efos: list[str]) -> bool:
    """Only with consent for every trait the score maps to. A score mapped to no trait at
    all cannot be consented to, so it cannot be scored."""
    return bool(efos) and all(consent.has_consented(person, category(e)) for e in efos)


def set_score_consent(person: str, pid: str, granted: bool) -> list[str]:
    """Grant or withdraw consent for every trait this score maps to; returns the trait ids.

    Withdrawing goes through consent.set_consent, which deletes every output that depended on
    the grant (gated_outputs below) and rebuilds the report."""
    efos = describe(pid)["efo_ids"]
    for efo in efos:
        if consent.has_consented(person, category(efo)) != granted:
            consent.set_consent(person, category(efo), granted)
    return efos


def matrix(pids: list[str], who: list[dict] | None = None) -> dict:
    """{people, scores, cells} for the consent matrix: cells[person][pid] is True only when
    the person may be scored with it."""
    who = who if who is not None else people()
    sens = sensitive_traits()
    scores = []
    for pid in pids:
        info = describe(pid)
        scores.append({**info, "label": trait_label(info),
                       "sensitive": bool(sens & set(info["efo_ids"]))})
    cells = {p["id"]: {s["pgs_id"]: may_score(p["id"], s["efo_ids"]) for s in scores}
             for p in who}
    return {"people": who, "scores": scores, "cells": cells}


# --- results on disk ----------------------------------------------------------------------
def result_path(callset: str, pid: str) -> Path:
    return config.PRS_CHOSEN_DIR / callset / f"{pid}.json"


def results(callset: str) -> list[dict]:
    d = config.PRS_CHOSEN_DIR / callset
    out = []
    for p in sorted(d.glob("PGS*.json")) if d.is_dir() else []:
        try:
            out.append(json.loads(p.read_text()))
        except (OSError, ValueError):
            log.warning("unreadable chosen-score result %s — ignored", p)
    return out


def scored_callsets() -> list[str]:
    root = config.PRS_CHOSEN_DIR
    return sorted(d.name for d in root.iterdir() if d.is_dir() and any(d.glob("PGS*.json"))) \
        if root.is_dir() else []


def report_outputs(callset: str) -> list[Path]:
    """The chosen report in every format tree: a rendering left behind keeps disclosing."""
    from .render import FORMATS
    md = reportpaths.chosen_report("prs", callset)
    return [md, *[reportpaths.derived(md, name, ext) for name, (ext, _r) in FORMATS.items()]]


def gated_outputs(sample: str, efo: str) -> list[Path]:
    """What a `pgs:<efo>` grant for this callset's person allows to exist: every result of a
    score mapped to that trait, and the report that lists them (rebuilt after a withdrawal
    from what remains)."""
    out = []
    for r in results(sample):
        if efo in (r.get("efo_ids") or []):
            out.append(result_path(sample, r["pgs_id"]))
    return [*out, *report_outputs(sample)] if out else out


def after_withdrawal(callsets) -> None:
    """Rebuild the reports a withdrawal just removed, from the results that remain."""
    for cs in sorted(set(callsets)):
        if results(cs) or reportpaths.chosen_report("prs", cs).exists():
            write_chosen_report(cs)


# --- readiness, fetching, bundles ----------------------------------------------------------
def _ref_files(pid: str, genomewide: bool) -> list[Path]:
    if genomewide:
        return [kg_dir(pid) / "scores_all.tsv", prs.ref_table(kg_dir(pid))]
    return [prs._kg_cache(pid)]


def _genomewide(pid: str) -> bool:
    return describe(pid)["n_variants"] > config.PRS_PLINK_MIN


def ready(pid: str) -> tuple[bool, str]:
    """Can this score be computed here without the network?"""
    if not scoring_path(pid).exists():
        return False, "no scoring file"
    gw = _genomewide(pid)
    if gw:
        if not (kg_dir(pid) / "scores_all.tsv").exists():
            return False, "no 1000 Genomes reference scores"
        pvars = list(kg_dir(pid).glob("chr*.pvar"))
        if not prs.ref_table(kg_dir(pid)).exists() and not pvars:
            return False, "no reference alleles"
        return True, "genome-wide, reference here"
    af = prs._kg_cache(pid)
    if not (af.exists() and af.stat().st_size > 0):
        return False, "no 1000 Genomes frequencies"
    return True, "reference here"


def _prune(pid: str) -> int:
    """Delete the reference panel's genotype intermediates once its scores and the REF table
    exist. They are only ever read to build those two files."""
    d = kg_dir(pid)
    if not ((d / "scores_all.tsv").exists() and prs.ref_table(d).exists()):
        return 0
    freed = 0
    for f in [*d.glob("chr*.pgen"), *d.glob("chr*.pvar"), *d.glob("chr*.psam"),
              *d.glob("chr*.log"), *d.glob("chr*.vcf.gz*")]:
        freed += f.stat().st_size
        f.unlink()
    return freed


fetch_errors: dict[str, str] = {}     # why the last fetch of each id failed, for callers
BATCH_MAX_ROWS = 60_000_000          # variant rows held in memory for one reference batch


def _bundle_done(pid: str, to: Path | None) -> bool:
    if not to:
        return False
    b, problems = verify_bundle(Path(to) / "bundles" / pid)
    return bool(b) and not problems


def _reference_batches(pids: list[str], to: Path | None, refresh: bool, size: int) -> None:
    """Build the reference data of the genome-wide scores in `pids` in batches of `size`
    (prs.reference_scores_batch), so each 1000G chromosome is read once per batch instead of
    once per score. The per-score loop in fetch then finds the caches and writes bundles."""
    batch: dict[str, list] = {}
    rows = 0

    def run(batch):
        import hashlib
        key = hashlib.sha256(",".join(sorted(batch)).encode()).hexdigest()[:12]
        wd = config.PGS_DIR / f".batch-{key}"
        try:
            done = prs.reference_scores_batch(batch, wd)
            log.info("batch %s: reference data for %d of %d score(s).", key, len(done),
                     len(batch))
        finally:
            shutil.rmtree(wd, ignore_errors=True)

    for pid in pids:
        if _bundle_done(pid, to) and not refresh:
            continue
        try:
            if refresh:
                for f in [scoring_path(pid), prs._kg_cache(pid)]:
                    f.unlink(missing_ok=True)
                shutil.rmtree(kg_dir(pid), ignore_errors=True)
            prs._scoring_file(pid)
            if engine_problem(pid):
                continue
            variants = prs._load_variants(scoring_path(pid))
        except Exception as e:  # noqa: BLE001 — the per-score loop reports it
            log.warning("%s: not batched (%s)", pid, e)
            continue
        if len(variants) <= config.PRS_PLINK_MIN:
            continue
        d = kg_dir(pid)
        if (d / "scores_all.tsv").exists() and prs.ref_table(d).exists():
            continue
        auto = prs._autosomal(variants)
        if batch and (len(batch) >= size or rows + len(auto) > BATCH_MAX_ROWS):
            run(batch)
            batch, rows = {}, 0
        batch[pid] = auto
        rows += len(auto)
    if batch:
        run(batch)


def fetch(pids: list[str], to: Path | None = None, refresh: bool = False,
          prune: bool = False, batch: int = 1, link: bool = False) -> list[str]:
    """Download each score and build its reference data here; with `to`, also write a bundle
    per score (and a copy of the catalog snapshot) for carrying to another machine.

    Needs the network and the reference tools (bcftools; plink2 for genome-wide scores), and
    nothing about the family. Returns the ids that are ready."""
    if config.PGS_OFFLINE:
        raise SystemExit("PGS_OFFLINE is set: fetch needs the network. Run it on a connected "
                         "machine and `pgs import` the bundles here.")
    require_tools("bcftools")
    if batch > 1:
        _reference_batches(pids, to, refresh, batch)
        refresh = False          # done above; the loop below must not discard what it built
    done = []
    for pid in pids:
        if to and not refresh:
            b, problems = verify_bundle(Path(to) / "bundles" / pid)
            if b and not problems:
                # Resumable: a long fetch interrupted (a spot VM reclaimed) picks up where it
                # stopped instead of rebuilding every reference panel it already finished.
                log.info("%s: bundle already written and intact — skipped.", pid)
                done.append(pid)
                continue
        try:
            if refresh:
                for f in [scoring_path(pid), prs._kg_cache(pid)]:
                    f.unlink(missing_ok=True)
                shutil.rmtree(kg_dir(pid), ignore_errors=True)
            prs._scoring_file(pid)
            problem = engine_problem(pid)
            if problem:
                raise ValueError(f"cannot be scored here: {problem}")
            variants = prs._load_variants(scoring_path(pid))
            if not variants:
                raise ValueError("no usable variants in the scoring file")
            auto = prs._autosomal(variants)
            if len(variants) > config.PRS_PLINK_MIN:
                log.info("%s: genome-wide (%d variants) — building the 1000 Genomes reference "
                         "distribution; this takes a while the first time.", pid, len(variants))
                if not prs._reference_scores(pid, auto):
                    raise RuntimeError("the 1000 Genomes reference scoring did not complete")
                if not prs.ref_table(kg_dir(pid)).exists() and not prs.export_ref_table(pid, auto):
                    raise RuntimeError("the reference panel's site list is incomplete")
                if prune:
                    freed = _prune(pid)
                    if freed:
                        log.info("%s: pruned %.1f GB of reference intermediates.", pid,
                                 freed / 1e9)
            elif not prs._reference_freqs(pid, auto):
                raise RuntimeError("no 1000 Genomes frequencies resolved")
            ok, why = ready(pid)
            if not ok:
                raise RuntimeError(why)
            if to:
                write_bundle(pid, Path(to), link=link)
            done.append(pid)
            log.info("%s ready.", pid)
        except Exception as e:  # noqa: BLE001 — one score must not stop the others
            fetch_errors[pid] = str(e)
            log.warning("%s: fetch failed (%s)", pid, e)
    if to:
        cur = pgs_catalog.current_dir()
        if cur is not None:
            dest = Path(to) / "catalog" / cur.name
            if not dest.exists():
                shutil.copytree(cur, dest)
    return done


def _bundle_files(pid: str) -> list[Path]:
    return [scoring_path(pid), *_ref_files(pid, _genomewide(pid))]


def write_bundle(pid: str, to: Path, link: bool = False) -> Path:
    """<to>/bundles/<pid>/: the score's files at their PGS_DIR-relative paths + bundle.json.

    `link` hard-links the files instead of copying them when `to` is on the same filesystem,
    for a machine building many bundles at once (the whole catalog's scoring files alone are
    ~90 GB). Off by default: a hard-linked bundle edited in place would change the cache."""
    info = describe(pid)
    d = Path(to) / "bundles" / pid
    shutil.rmtree(d, ignore_errors=True)
    files = {}
    for src in _bundle_files(pid):
        rel = src.relative_to(config.PGS_DIR).as_posix()
        dst = d / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        if link:
            try:
                os.link(src, dst)
            except OSError:
                shutil.copy2(src, dst)
        else:
            shutil.copy2(src, dst)
        files[rel] = sha256(dst)
    (d / BUNDLE_FILE).write_text(json.dumps({
        "format": BUNDLE_FORMAT, "pgs_id": pid, "created": today(), "engine": ENGINE,
        "efo_ids": info["efo_ids"], "traits": info["traits"],
        "n_variants": info["n_variants"], "hm_date": info["hm_date"],
        "kg_pop": config.KG_POP, "files": files}, indent=1) + "\n")
    log.info("Bundle written: %s", d)
    return d


def verify_bundle(d: Path) -> tuple[dict, list[str]]:
    try:
        b = json.loads((d / BUNDLE_FILE).read_text())
    except (OSError, ValueError):
        return {}, [f"no readable {BUNDLE_FILE}"]
    problems = []
    if b.get("format") != BUNDLE_FORMAT:
        problems.append(f"bundle format {b.get('format')} (this code reads {BUNDLE_FORMAT})")
    for rel, digest in (b.get("files") or {}).items():
        p = d / rel
        if ".." in Path(rel).parts or Path(rel).is_absolute():
            problems.append(f"{rel}: path escapes the bundle")
        elif not p.is_file():
            problems.append(f"{rel} missing")
        elif sha256(p) != digest:
            problems.append(f"{rel} does not match its checksum")
    if not (b.get("files") or {}):
        problems.append("the manifest lists no files")
    return b, problems


def import_dir(src: Path) -> dict:
    """Install a carried directory: catalog/<release>/ snapshots and bundles/<pid>/ bundles.
    Every file is checked against its manifest first; a damaged bundle is refused whole."""
    src = Path(src)
    out = {"catalog": [], "bundles": [], "refused": []}
    for snap in sorted((src / "catalog").glob("*")) if (src / "catalog").is_dir() else []:
        if (snap / pgs_catalog.MANIFEST).is_file():
            out["catalog"].append(pgs_catalog.install(snap).name)
    for d in sorted((src / "bundles").glob("*")) if (src / "bundles").is_dir() else []:
        b, problems = verify_bundle(d)
        if problems:
            out["refused"].append((d.name, problems))
            log.warning("%s: bundle refused — %s", d.name, "; ".join(problems))
            continue
        for rel in b["files"]:
            dst = config.PGS_DIR / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            if dst.exists() and sha256(dst) == b["files"][rel]:
                continue
            tmp = dst.with_name(f"{dst.name}.{os.getpid()}.part")
            shutil.copy2(d / rel, tmp)
            os.replace(tmp, dst)
        out["bundles"].append(b["pgs_id"])
        log.info("%s installed (%s).", b["pgs_id"], "; ".join(b.get("traits") or []))
    return out


# --- scoring ------------------------------------------------------------------------------
_tokens: dict[tuple, str] = {}


def _token(path: Path) -> str:
    """file_token, memoised per (path, size, mtime) for this process: a report rebuild checks
    every result's stamp, and the same multi-GB VCF must not be hashed once per score."""
    st = path.stat()
    key = (str(path), st.st_size, st.st_mtime_ns)
    if key not in _tokens:
        _tokens[key] = file_token(path)
    return _tokens[key]


def stamp(pid: str, callset: str) -> dict:
    """What a result depends on. Any change makes it stale."""
    from .call import sample_vcf
    from . import forcecall
    gw = _genomewide(pid)
    refs = {f.name: _token(f) for f in _ref_files(pid, gw) if f.exists()}
    vcf = sample_vcf(callset)
    fc = forcecall.forcecalled_vcf(callset)
    return {"engine": ENGINE, "hm_date": file_header(scoring_path(pid)).get("HmPOS_date", ""),
            "scoring_file": _token(scoring_path(pid)), "reference": refs,
            "kg_pop": config.KG_POP, "vcf": _token(vcf) if vcf.exists() else "",
            "forcecalled": _token(fc) if (not gw and fc.exists()) else ""}


def score_one(callset: str, pid: str) -> dict:
    """Compute one result with the panel's own engine, and save it. Callers check consent."""
    from .call import sample_vcf
    from . import forcecall
    vcf = sample_vcf(callset)
    if not vcf.exists():
        raise SystemExit(f"{vcf} missing — nothing to score for {callset}.")
    gt_vcfs = ([forcecall.forcecalled_vcf(callset), vcf] if forcecall.available(callset)
               else [vcf])
    problem = engine_problem(pid)
    if problem:
        raise ValueError(f"{pid} cannot be scored here: {problem}")
    info = describe(pid)
    out_dir = config.PRS_CHOSEN_DIR / callset
    out_dir.mkdir(parents=True, exist_ok=True)
    config.LOGS_DIR.mkdir(parents=True, exist_ok=True)     # plink2 logs its pgen build there
    r = prs._score_entry(callset, {"pgs_id": pid, "trait": trait_label(info)}, gt_vcfs, out_dir)
    result = {
        "pgs_id": pid, "traits": info["traits"], "efo_ids": info["efo_ids"],
        "callset": callset, "person": consent._person(callset),
        "raw": r["raw"], "used": r["used"], "total": r["total"],
        "pe": r["pe"], "pa": r["pa"], "method": r["method"],
        "scored": today(), "stamp": stamp(pid, callset),
        "consent": [category(e) for e in info["efo_ids"]],
    }
    path = result_path(callset, pid)
    tmp = path.with_name(f"{path.name}.{os.getpid()}.part")
    tmp.write_text(json.dumps(result, indent=1) + "\n")
    os.replace(tmp, path)
    return result


def score(pids: list[str], names: list[str] | None = None, all_callsets: bool = False,
          progress=None) -> dict:
    """Score each chosen score for each named person (default: everyone) who consented.

    Consent is checked per person BEFORE anything is downloaded; a score nobody may have is
    never fetched. Returns {scored: [(callset, pid)], skipped: [(person, pid, why)],
    failed: [(callset, pid, error)]}."""
    say = progress or (lambda msg: log.info("%s", msg))
    who = resolve_people(names)
    out = {"scored": [], "skipped": [], "failed": []}
    touched = set()
    for pid in pids:
        info = describe(pid)
        allowed = [p for p in who if may_score(p["id"], info["efo_ids"])]
        for p in who:
            if p not in allowed:
                out["skipped"].append((p["id"], pid, "no consent for " + trait_label(info)))
        if not allowed:
            say(f"{pid}: nobody has consented to {trait_label(info)} — nothing downloaded.")
            continue
        ok, why = ready(pid)
        if not ok:
            if config.PGS_OFFLINE:
                for p in allowed:
                    out["failed"].append((p["id"], pid, f"not available offline ({why})"))
                say(f"{pid}: {why} and PGS_OFFLINE is set — fetch it on a connected machine "
                    "and `pgs import` the bundle.")
                continue
            say(f"{pid}: fetching (consented: {', '.join(p['name'] for p in allowed)}) …")
            if pid not in fetch([pid]):
                why = fetch_errors.get(pid, "fetch failed")
                for p in allowed:
                    out["failed"].append((p["id"], pid, why))
                continue
        for p in allowed:
            for cs in (p["callsets"] if all_callsets else p["callsets"][:1]):
                try:
                    say(f"{pid} → {cs} …")
                    score_one(cs, pid)
                    out["scored"].append((cs, pid))
                    touched.add(cs)
                except (Exception, SystemExit) as e:  # noqa: BLE001 — isolate each result
                    out["failed"].append((cs, pid, str(e)))
                    log.warning("%s → %s failed: %s", pid, cs, e)
    for cs in sorted(touched):
        write_chosen_report(cs)
    return out


def status(callsets: list[str] | None = None) -> list[dict]:
    """Every result (or those of `callsets`) with whether it is current:
    [{callset, pgs_id, current, why}]."""
    out = []
    for cs in (callsets if callsets is not None else scored_callsets()):
        for r in results(cs):
            pid, why = r["pgs_id"], ""
            if not may_score(r.get("person") or cs, r.get("efo_ids") or []):
                why = "consent no longer recorded"
            elif not scoring_path(pid).exists():
                why = "scoring file gone"
            else:
                try:
                    now = stamp(pid, cs)
                    old = r.get("stamp") or {}
                    changed = [k for k in now if now[k] != old.get(k)]
                    why = f"changed: {', '.join(changed)}" if changed else ""
                except Exception as e:  # noqa: BLE001
                    why = f"cannot check ({e})"
            out.append({"callset": cs, "pgs_id": pid, "current": not why, "why": why})
    return out


def rescore_stale() -> dict:
    stale = [s for s in status() if not s["current"]]
    out = {"scored": [], "failed": []}
    touched = set()
    for s in stale:
        r = json.loads(result_path(s["callset"], s["pgs_id"]).read_text())
        if not may_score(r.get("person") or s["callset"], r.get("efo_ids") or []):
            continue
        try:
            score_one(s["callset"], s["pgs_id"])
            out["scored"].append((s["callset"], s["pgs_id"]))
            touched.add(s["callset"])
        except (Exception, SystemExit) as e:  # noqa: BLE001
            out["failed"].append((s["callset"], s["pgs_id"], str(e)))
    for cs in touched:
        write_chosen_report(cs)
    return out


# --- updates (online) ---------------------------------------------------------------------
def remote_header(pid: str, limit: int = 1 << 20) -> dict[str, str]:
    """The header of the Catalog's current harmonized file, reading only its first chunk."""
    req = urllib.request.Request(config.PGS_URL_TMPL.format(pid=pid),
                                 headers={"User-Agent": "genome-ledger"})
    dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
    text = ""
    with urllib.request.urlopen(req, timeout=60) as r:
        read = 0
        while read < limit:
            chunk = r.read(16384)
            if not chunk:
                break
            read += len(chunk)
            text += dec.decompress(chunk).decode("utf-8", "replace")
            lines = text.split("\n")
            if any(ln and not ln.startswith("#") for ln in lines[:-1]):
                break
    return parse_header(text.split("\n"))


def check_updates() -> list[str]:
    """What changed upstream for the scores already used here (needs the network for the
    re-harmonisation check; newer scores come from the local catalog snapshot)."""
    notes = []
    pids = sorted({r["pgs_id"] for cs in scored_callsets() for r in results(cs)})
    try:
        cat = pgs_catalog.load()
    except SystemExit:
        cat = None
    for pid in pids:
        local = file_header(scoring_path(pid)).get("HmPOS_date", "") \
            if scoring_path(pid).exists() else ""
        try:
            remote = remote_header(pid).get("HmPOS_date", "")
            if remote and remote != local:
                notes.append(f"{pid}: re-harmonised upstream ({local or '?'} → {remote}). "
                             f"Same model, new positions file: `pgs fetch {pid} --refresh`, "
                             "then `pgs score --stale`.")
        except Exception as e:  # noqa: BLE001 — offline or a transient error
            notes.append(f"{pid}: could not check upstream ({e}).")
        if cat is not None:
            newer = cat.newer(pid)
            if newer:
                notes.append(f"{pid}: {len(newer)} newer score(s) for the same trait since it "
                             f"was released: {', '.join(newer[:6])}"
                             f"{' …' if len(newer) > 6 else ''}. Not re-run: choosing another "
                             "model is a new choice (consent for the trait already covers it).")
    return notes


# --- the report ---------------------------------------------------------------------------
def _band(pe):
    return prs._band(pe) if pe is not None else "—"


def write_chosen_report(callset: str) -> Path | None:
    """Rebuild genome/<person>/chosen/prs_<callset>.md from the results on disk. With no
    results left, the report is removed in every format instead of written empty."""
    rows = results(callset)
    if not rows:
        for f in report_outputs(callset):
            f.unlink(missing_ok=True)
        return None
    try:
        cat = pgs_catalog.load()
    except SystemExit:
        cat = None
    current = {(s["callset"], s["pgs_id"]): s for s in status([callset])}
    table = ["| Trait | PGS ID | Raw score | EUR percentile | Band (EUR) | "
             "Variants (used / total) | Engine | Published performance | Training ancestry | "
             "Scored | State |",
             "|---|---|--:|--:|---|---|---|---|---|---|---|"]
    for r in sorted(rows, key=lambda r: ("; ".join(r["traits"]).lower(), r["pgs_id"])):
        s = cat.summary(r["pgs_id"]) if cat else {}
        perf = pgs_catalog.describe_performance(s.get("performance", {})) if s else "—"
        anc = (s.get("training_ancestry") or s.get("gwas_ancestry") or "—") if s else "—"
        st = current.get((callset, r["pgs_id"]), {})
        pct = f"{r['pe']:.1f}%" if r.get("pe") is not None else "—"
        table.append(
            f"| {'; '.join(r['traits'])} | {r['pgs_id']} | {r['raw']:.4f} | {pct} | "
            f"{_band(r.get('pe'))} | {r['used']} / {r['total']} | {r['method']} | {perf} | "
            f"{anc} | {r['scored']} | {'current' if st.get('current', True) else st.get('why')} |")
    lines = [
        f"# Polygenic scores you chose — {callset}",
        "",
        f"- **Chosen by the family from the whole PGS Catalog, not curated.** The curated panel "
        "is the separate **Polygenic risk scores** report; nothing here was checked against its "
        "selection criteria.",
        f"- {len(rows)} score(s), each computed only after this person consented to every "
        "trait it maps to. Withdrawing consent deletes the result.",
        f"- Percentiles are against the **{config.KG_POP}** 1000 Genomes reference, as in the "
        "panel report.",
        "",
        *table,
        "",
        "## How to read this",
        "- **Published performance** is the median over European evaluation cohorts of the "
        "score's AUROC or C-index **without** age, sex or clinical covariates (ancestry "
        "components allowed). When no such evaluation exists it says so: a figure from a model "
        "with age and sex overstates what the score alone can tell you.",
        "- The **EUR percentile** is relative genetic load against the European reference, not "
        "absolute risk. A score trained or evaluated in another ancestry may be miscalibrated "
        "here, whatever percentile it returns.",
        "- **State** is `current` while the scoring file, the reference data, the genotypes and "
        "the engine are unchanged since the score was computed; otherwise it names what "
        "changed. `run.py pgs score --stale` recomputes those.",
        "",
        "## Caveats",
        f"- {RELATIVES_NOTE}",
        "- Scores are computed by the same engine as the panel: small scores per variant with an "
        "analytic percentile, genome-wide scores with plink2 and an empirical percentile.",
        "- Research-grade, not diagnostic.",
        "",
    ]
    md = reportpaths.chosen_report("prs", callset)
    write_report(md, lines)
    log.info("Chosen-scores report: %s", md)
    return md
