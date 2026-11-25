"""Keep a bundle store in step with the PGS Catalog: `run.py pgs update --store DIR`.

Meant for a nightly timer. Each run costs one request to the Catalog's REST API; nothing else
happens unless the Catalog has published a new release since the local snapshot (about monthly,
usually 5-65 new scores, occasionally a dump of a thousand or more). On a new release:

  1. take the new catalog snapshot (pgs_catalog.sync);
  2. list the scores the engine can compute that the store does not hold yet, minus those
     already known to be uncomputable (STORE/refused.tsv: haplotype, interaction or per-dosage
     rows, found by reading the files), so they are not retried every month;
  3. build their bundles into STORE/<release>/ in batches (pgs_chosen.fetch), unless there are
     more than `max_new`: then build nothing and say so, because a dump that size is a job for
     a bigger machine;
  4. for the scores the family has actually been scored with, ask whether the Catalog
     re-harmonised their files (same model, new positions); rebuild those bundles and
     recompute the results that depended on them (consent still holds, so this is in scope);
  5. send one notification with what changed (notify.notify, which holds it during quiet hours).

The store is a directory of release folders, each `pgs import`-able:
    STORE/2026-11-24/bundles/PGS.../   STORE/2026-11-24/catalog/2026-11-24/   STORE/refused.tsv
"""

from __future__ import annotations

import json
import os
from pathlib import Path

from . import config, notify, pgs_catalog, pgs_chosen
from .util import log

REFUSED = "refused.tsv"


def bundled(store: Path) -> set[str]:
    """Score ids with an intact-looking bundle anywhere in the store (bundle.json present)."""
    out = set()
    for rel in Path(store).iterdir() if Path(store).is_dir() else []:
        b = rel / "bundles"
        if b.is_dir():
            out |= {d.name for d in b.iterdir() if (d / pgs_chosen.BUNDLE_FILE).is_file()}
    return out


def refused(store: Path) -> dict[str, str]:
    try:
        lines = (Path(store) / REFUSED).read_text().splitlines()
    except OSError:
        return {}
    out = {}
    for line in lines:
        pid, _tab, why = line.partition("\t")
        if pid.startswith("PGS"):
            out[pid] = why
    return out


def _write_refused(store: Path, rows: dict[str, str]) -> None:
    p = Path(store) / REFUSED
    tmp = p.with_name(f"{p.name}.{os.getpid()}.part")
    tmp.write_text("".join(f"{k}\t{v}\n" for k, v in sorted(rows.items())))
    os.replace(tmp, p)


def remote_release() -> dict:
    info = json.loads(pgs_catalog._get(config.PGS_REST_URL + "info"))
    return info.get("latest_release") or {}


def update(store: Path, max_new: int = 300, batch: int = 40, dry_run: bool = False,
           force: bool = False) -> dict:
    """One update run. Returns what it found and did; notifies only when something changed."""
    store = Path(store)
    latest = remote_release()
    release = latest.get("date") or ""
    cur = pgs_catalog.current_dir()
    previous = cur.name if cur else ""
    out = {"release": release, "previous": previous, "changed": bool(release) and release != previous,
           "built": [], "failed": [], "refused": [], "too_many": False, "reharmonised": []}
    if not out["changed"] and not force:
        log.info("PGS Catalog release %s unchanged — nothing to do.", previous or "?")
        return out
    if not store.is_dir():
        raise SystemExit(f"bundle store {store} is not reachable")
    snap = pgs_catalog.sync() if not dry_run else cur
    cat = pgs_catalog.load(snap)
    known_bad = refused(store)
    missing = sorted(set(pgs_catalog.select(cat)) - bundled(store) - set(known_bad))
    out["missing"] = missing
    log.info("PGS Catalog %s: %d computable score(s) not in the store yet.", release, len(missing))

    if len(missing) > max_new:
        out["too_many"] = True
    elif missing and not dry_run:
        to = store / release
        done = pgs_chosen.fetch(missing, to=to, batch=batch, prune=True)
        out["built"] = done
        for pid in missing:
            if pid in done:
                continue
            why = pgs_chosen.fetch_errors.get(pid, "")
            if "does not model" in why or "not on an additive" in why:
                known_bad[pid] = why
                out["refused"].append(pid)
            else:
                out["failed"].append((pid, why))
        if out["refused"]:
            _write_refused(store, known_bad)

    # Re-harmonised files among the scores the family was actually scored with.
    used = sorted({r["pgs_id"] for cs in pgs_chosen.scored_callsets()
                   for r in pgs_chosen.results(cs)})
    for pid in used:
        try:
            local = pgs_chosen.file_header(pgs_chosen.scoring_path(pid)).get("HmPOS_date", "") \
                if pgs_chosen.scoring_path(pid).exists() else ""
            remote = pgs_chosen.remote_header(pid).get("HmPOS_date", "")
        except Exception as e:  # noqa: BLE001 — offline or transient; next release retries
            log.warning("%s: could not check upstream (%s)", pid, e)
            continue
        if remote and local and remote != local:
            out["reharmonised"].append(pid)
    if out["reharmonised"] and not dry_run:
        pgs_chosen.fetch(out["reharmonised"], to=store / release, refresh=True, batch=batch,
                         prune=True)
        out["rescored"] = pgs_chosen.rescore_stale()

    msg = summary(out)
    out["message"] = msg
    if not dry_run:
        notify.notify(msg)
    return out


def summary(out: dict) -> str:
    parts = [f"PGS Catalog {out['release']} (was {out['previous'] or 'none'}):"]
    if out["too_many"]:
        parts.append(f"{len(out.get('missing', []))} new scores, too many to build here; "
                     "nothing built. Run a bigger batch build (pgs fetch --batch) on demand.")
    else:
        parts.append(f"{len(out['built'])} new score bundle(s) built")
        if out["refused"]:
            parts.append(f"{len(out['refused'])} refused (haplotype/interaction/per-dosage)")
        if out["failed"]:
            parts.append(f"{len(out['failed'])} FAILED, will retry next release: "
                         + ", ".join(p for p, _w in out["failed"][:8]))
    if out["reharmonised"]:
        n = len((out.get("rescored") or {}).get("scored", []))
        parts.append(f"re-harmonised upstream: {', '.join(out['reharmonised'])}; "
                     f"{n} result(s) recomputed")
    return "; ".join(parts[:1]) + " " + "; ".join(parts[1:]) + "."
