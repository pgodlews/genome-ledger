"""Bring the scores you chose up to date with the PGS Catalog: `run.py pgs update`.

Run it when you want to; it is not a timer job. It only ever touches scores someone here chose:

  1. If the Catalog has published a release newer than the local snapshot, take a new snapshot,
     so `pgs search`, `pgs show` and the picker see the new scores.
  2. For every score someone here was scored with, ask whether the Catalog re-harmonised its
     file (the same model with its positions re-mapped: `#HmPOS_date` in the header). If it did,
     download it again, rebuild its reference data and recompute the results that used it.
     Consent still holds: it is the same model for the same trait.
  3. List newer scores for the same traits. They are not run: switching models is a new choice.

It never downloads or builds a score nobody chose. On an air-gapped machine, do the network part
on a connected one (`pgs sync --to DIR`, `pgs fetch <ids> --refresh --to DIR`) and `pgs import`
the result, then `pgs score --stale`.
"""

from __future__ import annotations

import json

from . import config, pgs_catalog, pgs_chosen
from .util import log


def remote_release() -> dict:
    info = json.loads(pgs_catalog._get(config.PGS_REST_URL + "info"))
    return info.get("latest_release") or {}


def used_scores() -> list[str]:
    """Every score id with a result for someone here."""
    return sorted({r["pgs_id"] for cs in pgs_chosen.scored_callsets()
                   for r in pgs_chosen.results(cs)})


def update(dry_run: bool = False) -> dict:
    """One on-demand update. Returns what it found and did; with `dry_run`, only what it found."""
    if config.PGS_OFFLINE:
        raise SystemExit("PGS_OFFLINE is set: pgs update needs the network. On a connected "
                         "machine run `pgs sync --to DIR` and `pgs fetch <ids> --refresh --to DIR`, "
                         "then `pgs import DIR` and `pgs score --stale` here.")
    cur = pgs_catalog.current_dir()
    previous = cur.name if cur else ""
    release = remote_release().get("date") or ""
    out = {"release": release, "previous": previous,
           "new_release": bool(release) and release != previous,
           "used": used_scores(), "reharmonised": [], "unchecked": [], "rescored": {},
           "newer": {}, "dry_run": dry_run}
    if out["new_release"] and not dry_run:
        pgs_catalog.sync()

    for pid in out["used"]:
        path = pgs_chosen.scoring_path(pid)
        local = pgs_chosen.file_header(path).get("HmPOS_date", "") if path.exists() else ""
        try:
            remote = pgs_chosen.remote_header(pid).get("HmPOS_date", "")
        except Exception as e:  # noqa: BLE001 — offline or transient: say so, carry on
            log.warning("%s: could not check upstream (%s)", pid, e)
            out["unchecked"].append(pid)
            continue
        if remote and remote != local:
            out["reharmonised"].append(pid)
    if out["reharmonised"] and not dry_run:
        pgs_chosen.fetch(out["reharmonised"], refresh=True)
        out["rescored"] = pgs_chosen.rescore_stale()

    try:
        cat = pgs_catalog.load()
    except SystemExit:
        cat = None
    for pid in out["used"] if cat is not None else []:
        newer = cat.newer(pid)
        if newer:
            out["newer"][pid] = newer
    out["message"] = summary(out)
    return out


def summary(out: dict) -> str:
    lines = []
    if out["new_release"]:
        verb = "is available" if out["dry_run"] else "snapshot taken"
        lines.append(f"PGS Catalog release {out['release']} {verb} (was {out['previous'] or 'none'}).")
    else:
        lines.append(f"PGS Catalog release {out['previous'] or '?'} is current.")
    if not out["used"]:
        lines.append("No chosen scores have results yet, so there is nothing else to update.")
        return "\n".join(lines)
    if out["reharmonised"]:
        if out["dry_run"]:
            lines.append(f"Re-harmonised upstream, would be rebuilt and recomputed: "
                         f"{', '.join(out['reharmonised'])}.")
        else:
            r = out["rescored"] or {}
            lines.append(f"Re-harmonised upstream and rebuilt: {', '.join(out['reharmonised'])}; "
                         f"{len(r.get('scored', []))} result(s) recomputed"
                         + (f", {len(r['failed'])} failed" if r.get("failed") else "") + ".")
    else:
        lines.append(f"None of the {len(out['used'])} chosen score(s) was re-harmonised.")
    if out["unchecked"]:
        lines.append(f"Could not check upstream: {', '.join(out['unchecked'])}.")
    for pid, newer in out["newer"].items():
        lines.append(f"{pid}: {len(newer)} newer score(s) for the same trait: "
                     f"{', '.join(newer[:6])}{' …' if len(newer) > 6 else ''}. Not run: choosing "
                     "another model is a new choice.")
    return "\n".join(lines)
