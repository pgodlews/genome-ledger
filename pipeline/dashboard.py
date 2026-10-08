"""What a dashboard shows at a glance: the genome map, the headline numbers, the findings.

The navigation pages used to be lists of links. These summaries are read from the reports
themselves (the markdown is the source of truth, and a dashboard that re-derived a number
could disagree with the report it links to), plus the variant rows for the genome map and
the findings list, which need positions and descent that no report table carries.

Everything degrades by omission: a section whose report or data is missing is left out, so
a machine holding reports synced from elsewhere — no database — still gets a correct page.
"""

from __future__ import annotations

import re
from pathlib import Path

from . import config, reportpaths as rp
from .util import disease_name, log

# GRCh38 chromosome lengths in Mb — the widths of the genome map.
CHROM_MB = {"1": 248, "2": 242, "3": 198, "4": 190, "5": 182, "6": 171, "7": 159, "8": 145,
            "9": 138, "10": 134, "11": 135, "12": 133, "13": 114, "14": 107, "15": 102,
            "16": 90, "17": 83, "18": 80, "19": 59, "20": 64, "21": 47, "22": 51, "X": 156}


# --- reading the reports ---------------------------------------------------------
def sections(md_path: Path) -> list[tuple[str, list[str]]]:
    """[(heading, lines)] for a report: the text before the first `##`/`###` heading is
    filed under "", front matter dropped. [] when the file is missing."""
    from .render import split_front_matter
    try:
        _fields, body = split_front_matter(Path(md_path).read_text())
    except OSError:
        return []
    out: list[tuple[str, list[str]]] = [("", [])]
    for line in body.splitlines():
        if line.startswith(("## ", "### ")):
            out.append((line.lstrip("# ").strip(), []))
        else:
            out[-1][1].append(line)
    return out


def table_rows(lines: list[str]) -> list[list[str]]:
    """Body rows of the first markdown table in `lines`, as stripped cells. Placeholder
    rows (`| _none_ | | |`) are not rows."""
    rows = [ln for ln in lines if ln.strip().startswith("|") and ln.strip().endswith("|")]
    out = []
    for ln in rows[2:]:
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        if cells and cells[0] and not cells[0].startswith("_none_") and any(cells):
            out.append(cells)
    return out


def _rows_under(md: Path, prefix: str) -> list[list[str]]:
    for heading, lines in sections(md):
        if heading.lower().startswith(prefix.lower()):
            return table_rows(lines)
    return []


def panel_counts(callset: str, snap: str) -> tuple[list[list[str]], list[list[str]]]:
    """(disease-risk rows, carrier rows) from a callset's monogenic panel report."""
    md = rp.snapshot_report("panel", callset, snap)
    return _rows_under(md, "Disease risk"), _rows_under(md, "Carrier status")


def prs_percentiles(callset: str) -> list[tuple[str, float]]:
    """[(trait, percentile)] from a callset's PRS report, highest first."""
    out = []
    for _heading, lines in sections(rp.genome_report("prs", callset)):
        for cells in table_rows(lines):
            if len(cells) > 3:
                m = re.fullmatch(r"(\d+(?:\.\d+)?)%", cells[3])
                if m:
                    out.append((cells[0], float(m.group(1))))
    return sorted(out, key=lambda t: -t[1])


def delta_summary(snap: str) -> dict | None:
    """The latest DELTA report for a snapshot: counts by tier, the earlier snapshot it
    compares against, its rows, and its path."""
    fam_dir = rp.snapshot_md_dir(snap) / rp.FAMILY
    deltas = sorted(fam_dir.glob("delta_*.md")) if fam_dir.is_dir() else []
    if not deltas:
        return None
    md = deltas[-1]
    text = md.read_text()
    def n(pat: str) -> int:
        m = re.search(pat, text)
        return int(m.group(1)) if m else 0
    rows = [r for _h, lines in sections(md) for r in table_rows(lines) if len(r) >= 6]
    prev = md.stem[len("delta_"):].split("_to_")[0]
    return {"md": md, "prev": prev, "rows": rows,
            "actionable": n(r"(\d+) actionable now"), "monitor": n(r"(\d+) monitor"),
            "log": n(r"(\d+) log-only"), "total": n(r"Changes triaged:\*\* (\d+)")}


def findings(person: str, snap: str, fam) -> list[dict] | None:
    """The actionable variants this person carries, with where each came from. None when
    this machine's database has no rows for the snapshot (reports synced from elsewhere)."""
    from . import inheritance
    try:
        if not inheritance.loaded(snap):
            return None
        variants = inheritance.collect(snap, fam)
    except (Exception, SystemExit) as e:  # noqa: BLE001 — a dashboard never fails on this
        log.debug("dashboard: no variant rows for %s (%s)", snap, e)
        return None
    out = []
    for v in variants:
        c = v["carriers"].get(person)
        if not c:
            continue
        chrom, _, rest = v["label"].partition(":")
        pos = int(re.match(r"\d+", rest).group(0)) if re.match(r"\d+", rest) else 0
        out.append({**v, "zyg": c["zyg"], "chrom": chrom, "pos": pos,
                    "origin": inheritance.ORIGIN_LABEL[inheritance.origin(v, person, fam)],
                    "low_ab": c.get("low_ab")})
    return sorted(out, key=lambda f: (f["zyg"] != "HOM", f["gene"]))


# --- rendering -------------------------------------------------------------------
def _e(s) -> str:
    from .render import _esc
    return _esc(str(s))


def _href(md: Path, out_dir: Path) -> str:
    import os
    from .render import _attr
    return _attr(os.path.relpath(rp.html_for(md), out_dir))


def genome_map(found: list[dict]) -> str:
    """Every chromosome at its true relative length, each finding pinned where it sits."""
    if not found:
        return ""
    by_chrom: dict[str, list[dict]] = {}
    for f in found:
        by_chrom.setdefault(f["chrom"], []).append(f)
    bars = []
    for chrom, mb in CHROM_MB.items():
        pins = ""
        for f in by_chrom.get(chrom, []):
            left = max(0.0, min(100.0, f["pos"] / (mb * 1e6) * 100))
            cls = " hom" if f["zyg"] == "HOM" else ""
            pins += (f'<span class="lab{cls}" style="left:{left:.1f}%">{_e(f["gene"])}</span>'
                     f'<span class="pin{cls}" style="left:{left:.1f}%"></span>')
        bars.append(f'<div class="chr" style="flex-grow:{mb}" title="chr {chrom}">{pins}</div>')
    return ('<section class="card" aria-label="Genome map">'
            '<div style="display:flex;align-items:baseline;justify-content:space-between">'
            '<h2>Where the findings sit</h2><div class="legend2">'
            '<span class="hom">two copies</span><span>one copy</span></div></div>'
            f'<div class="gmap">{"".join(bars)}</div>'
            '<div class="gaxis"><span>chr 1</span><span>chr 6</span><span>chr 12</span>'
            '<span>chr 22 · X</span></div></section>')


def _tile(label: str, value, sub: str, cls: str = "") -> str:
    return (f'<div class="card"><p class="k">{_e(label)}</p>'
            f'<p class="big {cls}">{value}</p><p class="v">{_e(sub)}</p></div>')


def stat_tiles(callset: str, snap: str | None) -> str:
    tiles = []
    if snap:
        risk, carrier = panel_counts(callset, snap)
        if rp.snapshot_report("panel", callset, snap).exists():
            genes = lambda rows: " · ".join(dict.fromkeys(r[2] for r in rows if len(r) > 2))
            tiles.append(_tile("Disease-risk findings", len(risk), genes(risk) or "none",
                               "c-hom" if risk else "c-mut"))
            tiles.append(_tile("Carrier status", len(carrier), genes(carrier) or "none",
                               "c-het" if carrier else "c-mut"))
        d = delta_summary(snap)
        if d:
            mine = [r for r in d["rows"] if r[1] == callset]
            tiles.append(_tile(f"Changed since {d['prev'][:10]}", len(mine),
                               "see What changed" if mine else "nothing for this person",
                               "c-acc" if mine else "c-mut"))
    prs = prs_percentiles(callset)
    if prs:
        trait, pct = prs[0]
        tiles.append(_tile("Highest polygenic score",
                           f'{pct:.0f}<span style="font-size:1.3rem" class="c-mut">th</span>',
                           trait))
    if not tiles:
        return ""
    return f'<section class="grid g{min(4, max(2, len(tiles)))}">{"".join(tiles)}</section>'


def findings_card(found: list[dict] | None, has_tree: bool) -> str:
    if found is None:
        return ""
    items = []
    for f in found:
        hom = f["zyg"] == "HOM"
        sig = f["sig"].replace("_", " ")
        low = " Low allele fraction — review before acting on it." if f["low_ab"] else ""
        items.append(
            f'<div class="finding"><div class="z{" hom" if hom else ""}">'
            f'{"2×" if hom else "1×"}</div><div>'
            f'<p class="n">{_e(disease_name(f["disease"]) or f["gene"])}</p>'
            f'<p class="m">{_e(f["gene"])} · {_e(f["label"])} · {_e(sig)}</p>'
            f'<p class="s">{"Two copies" if hom else "One copy"}, {_e(f["origin"])}.{low}</p>'
            '</div></div>')
    link = ('<a href="inheritance.html" target="_top" style="font-size:.85rem">'
            'Trace through the family</a>') if has_tree else ""
    body = "".join(items) or '<p class="v">No pathogenic or likely-pathogenic variants.</p>'
    return ('<section class="card" aria-label="Findings">'
            '<div style="display:flex;align-items:baseline;justify-content:space-between">'
            f'<h2>Findings</h2>{link}</div>{body}</section>')


def track(pct: float) -> str:
    cls = "hi" if pct >= 90 else "up" if pct >= 75 else ""
    return (f'<span class="track"><i class="{cls}" style="left:{max(0, min(100, pct)):.0f}%">'
            '</i></span>')


def prs_card(callset: str, out_dir: Path) -> str:
    prs = prs_percentiles(callset)
    if not prs:
        return ""
    rows = "".join(f'<span>{_e(t[:1].upper() + t[1:])}</span>{track(p)}<b>{p:.0f}</b>'
                   for t, p in prs[:8])
    href = _href(rp.genome_report("prs", callset), out_dir)
    return ('<section class="card" aria-label="Polygenic scores">'
            '<div style="display:flex;align-items:baseline;justify-content:space-between">'
            f'<h2>Polygenic scores</h2><a href="{href}" style="font-size:.85rem">'
            f'{"All " + str(len(prs)) if len(prs) > 8 else "Full report"}</a></div>'
            f'<div class="prs">{rows}</div>'
            '<p class="v" style="margin-top:.9rem">Percentile against the reference '
            'population. Relative genetic load, not absolute risk.</p></section>')


def genome_cards(callset: str, out_dir: Path) -> str:
    """Three small cards: drug response, HLA, maternal line — each only if reported."""
    cards = []
    pgx = rp.genome_report("pgx", callset)
    rows = _rows_under(pgx, "Actionable")
    if pgx.exists():
        chips = "".join(f"<span>{_e(r[0])} {_e(r[1])}</span>" for r in rows[:4] if len(r) > 1)
        cards.append(f'<a class="card" href="{_href(pgx, out_dir)}"><p class="k">Drug response'
                     f'</p><p class="t">{len(rows)} non-normal call{"" if len(rows) == 1 else "s"}'
                     f'</p><div class="chips">{chips}</div></a>')
    hla = rp.genome_report("hla", callset)
    for heading, lines in sections(hla):
        if heading.lower().startswith("celiac"):
            text = " ".join(ln.strip() for ln in lines if ln.strip())
            m = re.search(r"\*\*(.+?)\*\*", text)
            first = re.sub(r"\*\*", "", text.split(". ")[0]).strip().rstrip(".")
            cards.append(f'<a class="card" href="{_href(hla, out_dir)}"><p class="k">HLA</p>'
                         f'<p class="t">Celiac: {_e(m.group(1) if m else "see report")}</p>'
                         f'<p class="v">{_e(first)}.</p></a>')
            break
    anc = rp.genome_report("ancestry", callset)
    for _heading, lines in sections(anc):
        for ln in lines:
            m = re.search(r"mtDNA haplogroup\*\* \(maternal line\): \*\*(.+?)\*\*", ln)
            if m:
                hg = m.group(1).split(" (")[0]
                cards.append(f'<a class="card" href="{_href(anc, out_dir)}"><p class="k">'
                             f'Maternal line</p><p class="t">mtDNA {_e(hg)}</p>'
                             '<p class="v">Shared along the maternal line.</p></a>')
                break
    if not cards:
        return ""
    return f'<section class="grid g{min(3, max(2, len(cards)))}">{"".join(cards)}</section>'


def overview(person: str, callset: str, snap: str | None, fam, out_dir: Path) -> str:
    """The at-a-glance part of a person's page."""
    found = findings(person, snap, fam) if snap else None
    two = "".join(x for x in (findings_card(found, (out_dir / "inheritance.html").exists()),
                              prs_card(callset, out_dir)) if x)
    n_cols = (1 if found is not None else 0) + (1 if prs_percentiles(callset) else 0)
    return "".join((
        genome_map(found or []),
        '<div style="height:1rem"></div>' if found else "",
        stat_tiles(callset, snap),
        f'<div class="grid g{max(1, n_cols)}">{two}</div>' if two else "",
        genome_cards(callset, out_dir),
    ))


def people_cards(fam, snap: str | None, base: str) -> str:
    """One card per person for the landing page: findings at a glance, link to their page."""
    from .render import _attr
    cards = []
    for pid, p in fam.people.items():
        rel = f'<span class="rel" data-rel-of="{_attr(pid)}"></span>'
        if not fam.is_sequenced(pid):
            cards.append(f'<div class="card" data-person="{_attr(pid)}"><p class="nm">'
                         f'<span class="unseq">{_e(p.display_name)}</span>{rel}</p>'
                         '<p class="v">not sequenced</p></div>')
            continue
        badges = genes = ""
        if snap:
            risk, carrier = panel_counts(pid, snap)
            if risk:
                badges += f'<span class="badge hom">{len(risk)} disease-risk</span> '
            if carrier:
                badges += f'<span class="badge het">{len(carrier)} carrier</span>'
            genes = " · ".join(dict.fromkeys(r[2] for r in risk + carrier if len(r) > 2))
        n = len(fam.callsets(pid))
        extra = f"{n} callsets" if n > 1 else ""
        cards.append(
            f'<a class="card" data-person="{_attr(pid)}" '
            f'href="{_attr(base)}index_{_attr(pid)}.html"><p class="nm">'
            f'<span>{_e(p.display_name)}</span>{rel}</p><div>{badges or "&nbsp;"}</div>'
            f'<p class="g">{_e(" · ".join(x for x in (genes, extra) if x)) or "&nbsp;"}</p></a>')
    return f'<section class="grid g3 people">{"".join(cards)}</section>'


def changes_card(snap: str | None, base: str, out_dir: Path) -> str:
    d = delta_summary(snap) if snap else None
    if not d:
        return ""
    href = _href(d["md"], out_dir)
    cell = lambda n, label, cls: (f'<div><p class="big {cls if n else "c-mut"}">{n}</p>'
                                  f'<p class="v">{label}</p></div>')
    return (f'<a class="card" href="{href}"><p class="k">What changed · '
            f'{_e(d["prev"][:10])} → {_e(snap[:10])}</p>'
            '<div class="grid g3" style="margin:0">'
            + cell(d["actionable"], "actionable now", "c-hom")
            + cell(d["monitor"], "to monitor", "c-het")
            + cell(d["log"], "logged only", "c-acc") + "</div></a>")


def hero(fam) -> tuple[str, str]:
    """(headline, sentence) for the landing page: how many genomes over how many generations."""
    from . import inheritance
    words = ["No", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine",
             "Ten", "Eleven", "Twelve"]
    w = lambda n: words[n] if n < len(words) else str(n)
    n = sum(1 for p in fam.people if fam.is_sequenced(p))
    try:
        gens = len(set(inheritance.generations(fam).values()))
    except Exception:  # noqa: BLE001
        gens = 1
    head = f"{w(n)} genome{'' if n == 1 else 's'}"
    if gens > 1:
        head += f",<br>{w(gens).lower()} generations"
    return head + ".", "Re-read against the latest variant database."
