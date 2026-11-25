"""`run.py pgs …`: search the PGS Catalog offline, record consent, fetch, carry and score.

    pgs sync                        download the Catalog's metadata (connected machine)
    pgs categories | search | show  browse it, offline
    pgs consent PGS… [--grant/--revoke WHO]   the family × score consent matrix
    pgs score PGS… [--people WHO]   score whoever consented (fetches on demand when online)
    pgs fetch PGS… --to DIR         build bundles to carry to an air-gapped machine
    pgs import DIR                  install carried bundles / catalog snapshots
    pgs status | check | report     what is current, what changed upstream, rebuild reports
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import config, pgs_catalog, pgs_chosen


def add_parser(sub) -> None:
    sp = sub.add_parser("pgs", help="polygenic scores you choose from the whole PGS Catalog "
                                    "(search offline, consent per trait, bundles for "
                                    "air-gapped machines)")
    ps = sp.add_subparsers(dest="pgs_cmd", required=True)

    p = ps.add_parser("sync", help="download the Catalog's metadata (needs the network)")
    p.add_argument("--to", type=Path, help="write the snapshot under DIR/catalog/ to carry it")
    p.add_argument("--force", action="store_true")

    ps.add_parser("categories", help="trait categories with trait and score counts")

    p = ps.add_parser("search", help="find traits, or list a trait's scores")
    p.add_argument("term", nargs="?", default="")
    p.add_argument("--category", default="")
    p.add_argument("--trait", metavar="EFO", help="list the scores for this trait id")
    p.add_argument("--limit", type=int, default=40)

    p = ps.add_parser("list", help="score ids by category and evidence, one per line (for "
                                   "`fetch --from-file`); a summary goes to stderr")
    p.add_argument("--category", action="append", default=[],
                   help="a trait category (repeatable); 'disease' = every disease, disorder and "
                        "cancer category")
    p.add_argument("--evidence", default="none", choices=pgs_catalog.EVIDENCE,
                   help="european: a usable European evaluation; covariate-free: one without "
                        "age/sex/clinical covariates")
    p.add_argument("--engine", default="all", choices=["all", "small", "genome-wide"])
    p.add_argument("--include-unsupported", action="store_true")

    p = ps.add_parser("show", help="one score: what it is and every published evaluation")
    p.add_argument("pgs_id")

    p = ps.add_parser("consent", help="show the family × score consent matrix; grant/revoke")
    p.add_argument("pgs_ids", nargs="+")
    p.add_argument("--grant", metavar="WHO", help="comma-separated people, or 'all'")
    p.add_argument("--revoke", metavar="WHO", help="comma-separated people, or 'all'")

    p = ps.add_parser("score", help="score consented people (fetches on demand unless "
                                    "PGS_OFFLINE=1)")
    p.add_argument("pgs_ids", nargs="*")
    p.add_argument("--people", metavar="WHO", help="comma-separated people (default: everyone "
                                                   "who consented)")
    p.add_argument("--all-callsets", action="store_true",
                   help="every callset of each person, not only the primary one")
    p.add_argument("--stale", action="store_true", help="recompute results that are no longer "
                                                        "current")

    p = ps.add_parser("fetch", help="download scores and build their reference data; with "
                                    "--to, write bundles to carry (needs the network)")
    p.add_argument("pgs_ids", nargs="*")
    p.add_argument("--from-file", type=Path, help="PGS ids, one per line (# comments allowed)")
    p.add_argument("--shard", metavar="I/N", help="take every N-th id starting at I (1-based), "
                                                  "to split a long list across machines")
    p.add_argument("--to", type=Path, help="write bundles under DIR/bundles/")
    p.add_argument("--refresh", action="store_true", help="discard the cached copy first")
    p.add_argument("--batch", type=int, default=1, metavar="N",
                   help="build the 1000 Genomes reference of up to N genome-wide scores in one "
                        "pass over each chromosome (much faster for long lists; same results)")
    p.add_argument("--link", action="store_true",
                   help="hard-link bundle files to the cache instead of copying (same filesystem; "
                        "for a machine building many bundles)")
    p.add_argument("--prune", action="store_true",
                   help="delete the reference panel's genotype intermediates once a "
                        "genome-wide score's bundle data exists (they can run to GBs per score)")

    p = ps.add_parser("import", help="install a carried directory (bundles/ and catalog/)")
    p.add_argument("dir", type=Path)

    ps.add_parser("status", help="every chosen result and whether it is current")
    ps.add_parser("check", help="upstream changes: re-harmonised files, newer scores (network)")
    p = ps.add_parser("report", help="rebuild the chosen-scores reports")
    p.add_argument("callsets", nargs="*")


def _who(text: str | None) -> list[str] | None:
    if not text or text.strip().lower() == "all":
        return None
    return [x for x in (t.strip() for t in text.split(",")) if x]


def _ids(args) -> list[str]:
    ids = list(args.pgs_ids or [])
    if getattr(args, "from_file", None):
        for line in args.from_file.read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if line:
                ids.append(line)
    ids = [i.strip().upper() for i in ids if i.strip()]
    ids = list(dict.fromkeys(ids))
    if getattr(args, "shard", None):
        try:
            i, n = (int(x) for x in args.shard.split("/"))
            if not 1 <= i <= n:
                raise ValueError
        except ValueError:
            raise SystemExit("--shard takes I/N with 1 <= I <= N, e.g. 3/50")
        ids = ids[i - 1::n]
    bad = [i for i in ids if not i.startswith("PGS") or not i[3:].isdigit()]
    if bad:
        raise SystemExit(f"not PGS ids: {', '.join(bad)}")
    return ids


def print_matrix(m: dict) -> None:
    scores, who = m["scores"], m["people"]
    print("Scores:")
    for i, s in enumerate(scores, 1):
        gw = "genome-wide" if s["n_variants"] > config.PRS_PLINK_MIN else "small"
        flag = "  ⚠ severe, no established intervention" if s["sensitive"] else ""
        print(f"  [{i}] {s['pgs_id']}  {s['label']}  ({s['n_variants']:,} variants, {gw}){flag}")
    print()
    width = max([len(p["name"]) for p in who] + [6])
    print(" " * (width + 4) + "  ".join(f"[{i}]" for i in range(1, len(scores) + 1)))
    for n, p in enumerate(who, 1):
        cells = "  ".join(" ✓ " if m["cells"][p["id"]][s["pgs_id"]] else " · "
                          for s in scores)
        print(f"{n:>2}. {p['name']:<{width}}  {cells}")
    print("\n✓ = consented to every trait the score maps to; · = not (nothing is computed).")
    print(f"\n{pgs_chosen.CHILD_WARNING}\n{pgs_chosen.RELATIVES_NOTE}\n")


def _interactive(m: dict) -> None:
    who = m["people"]
    for i, s in enumerate(m["scores"], 1):
        ans = input(f"[{i}] {s['pgs_id']} {s['label']} — tick people (e.g. 1,3), untick "
                    "with - (e.g. -2), 'all', Enter = no change: ").strip().lower()
        if not ans:
            continue
        picks = range(1, len(who) + 1) if ans == "all" else []
        if ans != "all":
            picks = [x.strip() for x in ans.split(",") if x.strip()]
        for x in picks:
            neg = isinstance(x, str) and x.startswith("-")
            try:
                n = int(str(x).lstrip("-"))
                person = who[n - 1]
            except (ValueError, IndexError):
                print(f"  ignored {x!r}")
                continue
            pgs_chosen.set_score_consent(person["id"], s["pgs_id"], not neg)
            print(f"  {person['name']}: {'revoked' if neg else 'granted'}")


def _print_scores(rows: list[dict]) -> None:
    print(f"{'PGS ID':<10} {'variants':>10} {'engine':<7} {'released':<10} {'licence':<22} "
          "published performance | training ancestry")
    for s in rows:
        ready = "  [here]" if pgs_chosen.ready(s["pgs_id"])[0] else ""
        if s.get("unsupported"):
            ready += f"  [cannot score: {s['unsupported']}]"
        print(f"{s['pgs_id']:<10} {s['n_variants']:>10,} {s['engine']:<7} "
              f"{s['released'][:10]:<10} {s['licence']:<22} "
              f"{pgs_catalog.describe_performance(s['performance'])} | "
              f"{s['training_ancestry'] or s['gwas_ancestry'] or '—'}{ready}")


def run(args, refresh_pages=None) -> int:
    cmd = args.pgs_cmd
    if cmd == "sync":
        d = pgs_catalog.sync(args.to / "catalog" if args.to else None, force=args.force)
        print(d)
        return 0
    if cmd == "categories":
        for c in pgs_catalog.load().categories():
            print(f"{c['traits']:>4} traits {c['scores']:>5} scores  {c['label']}")
        return 0
    if cmd == "search":
        cat = pgs_catalog.load()
        efo = args.trait
        traits = [] if efo else cat.traits(args.term, args.category)
        if not efo and len(traits) == 1:
            efo = traits[0]["efo"]
        if not efo and args.term:
            exact = [t for t in traits if t["label"].lower() == args.term.strip().lower()]
            efo = exact[0]["efo"] if len(exact) == 1 else None
        if efo:
            print(f"{cat.trait_label.get(efo, efo)} ({efo}) — Catalog release {cat.release}\n")
            _print_scores(cat.scores_for(efo))
            return 0
        for t in traits[:args.limit]:
            print(f"{t['efo']:<16} {t['scores']:>4} score(s)  {t['label']}"
                  f"{'  [' + ', '.join(t['categories']) + ']' if t['categories'] else ''}")
        if len(traits) > args.limit:
            print(f"… {len(traits) - args.limit} more; narrow the term or raise --limit")
        if not traits:
            print("no matching trait")
        return 0
    if cmd == "list":
        cat = pgs_catalog.load()
        cats = []
        for c in args.category:
            cats += list(pgs_catalog.DISEASE_CATEGORIES) if c.lower() == "disease" else [c]
        ids = pgs_catalog.select(cat, cats or None, args.evidence, args.engine,
                                 args.include_unsupported)
        n = [cat.summary(p)["n_variants"] for p in ids]
        gw = sum(1 for x in n if x > config.PRS_PLINK_MIN)
        print("\n".join(ids))
        print(f"{len(ids)} scores ({gw} genome-wide, {len(ids) - gw} small), "
              f"{sum(n) / 1e6:,.0f}M variant rows — Catalog release {cat.release}",
              file=sys.stderr)
        return 0
    if cmd == "show":
        cat = pgs_catalog.load()
        s = cat.summary(args.pgs_id.upper())
        if not s:
            raise SystemExit(f"{args.pgs_id} is not in Catalog release {cat.release}")
        pub = s["publication"]
        print(f"{s['pgs_id']}  {s['name']}  — {s['reported_trait']}")
        print(f"  traits: {', '.join(f'{t} ({e})' for t, e in zip(s['traits'], s['efo_ids']))}")
        print(f"  {s['n_variants']:,} variants ({s['engine']}), method {s['method']}, "
              f"weights {s['weight_type']}, licence {s['licence']}, released {s['released']}")
        if s["unsupported"]:
            print(f"  ⚠ cannot be scored by this engine: {s['unsupported']}")
        print(f"  {pub['first_author']} {pub['year']} {pub['journal']} {pub['doi']}")
        print(f"  GWAS ancestry: {s['gwas_ancestry'] or '—'}; "
              f"training: {s['training_ancestry'] or '—'}")
        print(f"  published performance: {pgs_catalog.describe_performance(s['performance'])}")
        print(f"  here: {'yes' if pgs_chosen.ready(s['pgs_id'])[0] else 'no'}")
        print("\n  evaluations:")
        for e in cat.evaluations(s["pgs_id"]):
            m = ", ".join(f"{k}={e[k]}" for k in ("auroc", "c_index", "or", "hr") if e[k])
            print(f"   - {e['ancestry'] or '?'} n={e['n']} cases={e['cases']} "
                  f"[{e['covariate_class']}: {e['covariates'][:50] or 'none'}] {m} "
                  f"{e['other'][:60]} — {e['publication']}")
        newer = cat.newer(s["pgs_id"])
        if newer:
            print(f"\n  newer scores for the same trait: {', '.join(newer[:10])}")
        return 0
    if cmd == "consent":
        ids = _ids(args)
        for flag, granted in ((args.grant, True), (args.revoke, False)):
            if flag is None:
                continue
            for person in pgs_chosen.resolve_people(_who(flag)):
                for pid in ids:
                    pgs_chosen.set_score_consent(person["id"], pid, granted)
                if not granted and refresh_pages:
                    refresh_pages(person["id"])
        m = pgs_chosen.matrix(ids)
        print_matrix(m)
        if args.grant is None and args.revoke is None and sys.stdin.isatty():
            _interactive(m)
            print_matrix(pgs_chosen.matrix(ids))
        return 0
    if cmd == "score":
        if args.stale:
            out = pgs_chosen.rescore_stale()
            print(f"recomputed {len(out['scored'])}, failed {len(out['failed'])}")
            for cs, pid, err in out["failed"]:
                print(f"  {pid} → {cs}: {err}")
            return 1 if out["failed"] else 0
        ids = _ids(args)
        if not ids:
            raise SystemExit("name at least one PGS id (or --stale)")
        out = pgs_chosen.score(ids, _who(args.people), args.all_callsets, progress=print)
        for cs, pid in out["scored"]:
            print(f"scored {pid} → {cs}")
        for who, pid, why in out["skipped"]:
            print(f"skipped {pid} for {who}: {why}")
        for cs, pid, err in out["failed"]:
            print(f"FAILED {pid} → {cs}: {err}")
        if refresh_pages:
            for cs in {cs for cs, _p in out["scored"]}:
                refresh_pages(cs)
        return 1 if out["failed"] else 0
    if cmd == "fetch":
        ids = _ids(args)
        if not ids:
            raise SystemExit("name PGS ids, or --from-file")
        done = pgs_chosen.fetch(ids, args.to, args.refresh, args.prune, args.batch, args.link)
        print(f"{len(done)} of {len(ids)} ready" + (f"; bundles in {args.to / 'bundles'}"
                                                    if args.to else ""))
        missing = [i for i in ids if i not in done]
        if missing:
            print("failed: " + ", ".join(missing))
        return 0 if not missing else 1
    if cmd == "import":
        out = pgs_chosen.import_dir(args.dir)
        print(f"catalog: {', '.join(out['catalog']) or '—'}; bundles installed: "
              f"{len(out['bundles'])}; refused: {len(out['refused'])}")
        for pid, problems in out["refused"]:
            print(f"  {pid}: {'; '.join(problems)}")
        return 1 if out["refused"] else 0
    if cmd == "status":
        rows = pgs_chosen.status()
        if not rows:
            print("no chosen scores computed yet")
        for r in rows:
            print(f"{r['callset']:<20} {r['pgs_id']:<10} "
                  f"{'current' if r['current'] else r['why']}")
        return 0
    if cmd == "check":
        notes = pgs_chosen.check_updates()
        print("\n".join(notes) if notes else "nothing changed upstream")
        return 0
    if cmd == "report":
        for cs in args.callsets or pgs_chosen.scored_callsets():
            print(pgs_chosen.write_chosen_report(cs) or f"{cs}: no results, report removed")
        return 0
    return 2


def main(argv=None) -> int:   # for `python -m pipeline.pgs_cli` in tests
    p = argparse.ArgumentParser()
    add_parser(p.add_subparsers(dest="cmd"))
    return run(p.parse_args(argv))
