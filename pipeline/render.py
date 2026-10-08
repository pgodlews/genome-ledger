"""Render the pipeline's markdown reports to styled, self-contained HTML.

Markdown stays the source of truth; this produces a sibling `.html` for each report —
a single file with all CSS embedded (no external assets, no web fonts), so it is
portable, emailable, archivable, prints cleanly, and stays readable for decades.

Zero dependencies on purpose: the report markdown is a small, fully machine-generated
subset (headings, GFM tables, bold, inline code, bullets, full-line italic notes, a few
emoji), so a stdlib-only converter is total and rot-free. It is NOT a general markdown
engine — it handles exactly what the report modules emit.
"""

from __future__ import annotations

import html as _html
import logging
import os
import re
from pathlib import Path

from . import reportpaths, theme

log = logging.getLogger("genome")

# --- inline formatting -------------------------------------------------------
_CODE = re.compile(r"`([^`]+)`")
_BOLD = re.compile(r"\*\*([^*]+)\*\*")


def _esc(s: str) -> str:
    """Escape &, <, > — load-bearing: variant cells like `T>G` and transition arrows
    must not corrupt the HTML. Run BEFORE the bold/code regexes."""
    return _html.escape(s, quote=False)


def _attr(s: str) -> str:
    """Escape a value going into a quoted HTML attribute — quotes included.

    `_esc` deliberately leaves quotes alone (report bodies are full of them and they are
    harmless in text), so it is the wrong tool for an `href` or an `id`. Sample ids reach
    both, and they come from family.tsv, which is hand-edited."""
    return _html.escape(s, quote=True)


def _inline(text: str) -> str:
    t = _esc(text)
    # Stash code spans BEFORE bold: _BOLD would otherwise also fire inside the inserted
    # <code>…</code> HTML (`**x**` in backticks must stay literal).
    codes: list[str] = []

    def _stash(m: re.Match) -> str:
        codes.append(m.group(1))
        return f"\x00{len(codes) - 1}\x00"

    t = _CODE.sub(_stash, t)
    t = _BOLD.sub(r"<strong>\1</strong>", t)  # **bold** (single * star-alleles are safe)
    return re.sub(r"\x00(\d+)\x00",
                  lambda m: f"<code>{codes[int(m.group(1))]}</code>", t)
    # Note: no inline _italic_ — underscores appear inside ClinVar terms
    # (Likely_pathogenic, Conflicting_classifications_…). Full-line _…_ notes are
    # handled at the block level instead (see _convert_body).


# --- GFM tables --------------------------------------------------------------
def _cells(line: str) -> list[str]:
    s = line.strip()
    s = s[1:] if s.startswith("|") else s
    s = s[:-1] if s.endswith("|") else s
    return [c.strip() for c in s.split("|")]


def _is_separator(cells: list[str]) -> bool:
    return bool(cells) and all(re.fullmatch(r":?-{2,}:?", c or "-") for c in cells) \
        and any("-" in c for c in cells)


def _align(sep_cell: str) -> str:
    left, right = sep_cell.startswith(":"), sep_cell.endswith(":")
    return "center" if left and right else "right" if right else "left" if left else ""


_PCT_CELL = re.compile(r"(\d+(?:\.\d+)?)%")


def _table(rows: list[str]) -> str:
    header = _cells(rows[0])
    aligns, body_start = [""] * len(header), 1
    if len(rows) > 1 and _is_separator(_cells(rows[1])):
        aligns = [_align(c) for c in _cells(rows[1])]
        body_start = 2
    # A column of percentiles reads faster as a position on a track than as a number.
    pct_cols = {i for i, h in enumerate(header) if "percentile" in h.lower()}

    def _style(i: int) -> str:
        a = aligns[i] if i < len(aligns) else ""
        return f' style="text-align:{a}"' if a and i not in pct_cols else ""

    def _cell(i: int, c: str) -> str:
        m = _PCT_CELL.fullmatch(c.strip()) if i in pct_cols else None
        if not m:
            return _inline(c)
        pct = float(m.group(1))
        return f'<span class="pct"><b>{_esc(c.strip())}</b>{dashboard_track(pct)}</span>'

    out = ['<div class="tablewrap">', "<table>", "<thead><tr>"]
    out += [f"<th{_style(i)}>{_inline(h)}</th>" for i, h in enumerate(header)]
    out.append("</tr></thead><tbody>")
    for r in rows[body_start:]:
        cs = _cells(r)
        out.append("<tr>")
        out += [f"<td{_style(i)}>{_cell(i, c)}</td>" for i, c in enumerate(cs)]
        out.append("</tr>")
    out += ["</tbody></table>", "</div>"]
    return "\n".join(out)


def dashboard_track(pct: float) -> str:
    cls = "hi" if pct >= 90 else "up" if pct >= 75 else ""
    return (f'<span class="track"><i class="{cls}" style="left:{max(0.0, min(100.0, pct)):.0f}%">'
            "</i></span>")


# --- block converter ---------------------------------------------------------
def split_front_matter(md: str) -> tuple[dict[str, str], str]:
    """(fields, body) for a report that opens with a `---` YAML block; ({}, md) otherwise.
    Only the flat `key: value` form the pipeline itself writes (util.front_matter)."""
    if not md.startswith("---\n"):
        return {}, md
    end = md.find("\n---\n", 4)
    if end < 0:
        return {}, md
    fields = {}
    for line in md[4:end].splitlines():
        key, sep, value = line.partition(":")
        if sep:
            value = value.strip()
            # util.front_matter writes the title as a JSON string; hand back the string.
            if len(value) > 1 and value[0] == value[-1] == '"':
                try:
                    import json
                    value = json.loads(value)
                except ValueError:
                    pass
            fields[key.strip()] = value
    return fields, md[end + 5:]


def _meta_line(fields: dict[str, str]) -> str:
    """The front matter, shown as one quiet line under the title in the HTML."""
    bits = [fields.get("callset") or fields.get("scope", ""),
            f"snapshot {fields['snapshot']}" if "snapshot" in fields else "",
            f"ClinVar {fields['clinvar_date']}" if "clinvar_date" in fields else "",
            f"generated {fields['generated']}" if "generated" in fields else ""]
    bits = [b for b in bits if b]
    return f'<p class="note">{_esc(" · ".join(bits))}</p>' if bits else ""


def _convert_body(md: str) -> str:
    fields, md = split_front_matter(md)
    meta_line = _meta_line(fields)
    lines = md.splitlines()
    out: list[str] = []
    para: list[str] = []
    bullets: list[str] = []

    def flush_para() -> None:
        if para:
            text = " ".join(para).strip()
            if len(text) > 1 and text.startswith("_") and text.endswith("_"):
                out.append(f'<p class="note">{_inline(text[1:-1])}</p>')
            else:
                out.append(f"<p>{_inline(text)}</p>")
            para.clear()

    def flush_bullets() -> None:
        if bullets:
            out.append("<ul>")
            out.extend(f"<li>{_inline(b)}</li>" for b in bullets)
            out.append("</ul>")
            bullets.clear()

    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if s.startswith("|") and s.endswith("|"):
            flush_para(); flush_bullets()
            tbl = [lines[i]]
            i += 1
            while i < len(lines) and lines[i].strip().startswith("|"):
                tbl.append(lines[i]); i += 1
            out.append(_table(tbl))
            continue
        if not s:
            flush_para(); flush_bullets()
        elif s.startswith("### "):
            flush_para(); flush_bullets(); out.append(f"<h3>{_inline(s[4:])}</h3>")
        elif s.startswith("## "):
            flush_para(); flush_bullets(); out.append(f"<h2>{_inline(s[3:])}</h2>")
        elif s.startswith("# "):
            flush_para(); flush_bullets(); out.append(f"<h1>{_inline(s[2:])}</h1>")
            if meta_line:
                out.append(meta_line); meta_line = ""
        elif s == "---":
            flush_para(); flush_bullets(); out.append("<hr>")
        elif s.startswith("> "):
            # A markdown blockquote is how reports carry a notice (withheld traits, demo
            # data). Consecutive quoted lines are one callout.
            flush_para(); flush_bullets()
            quoted = [s[2:]]
            while i + 1 < len(lines) and lines[i + 1].strip().startswith("> "):
                i += 1
                quoted.append(lines[i].strip()[2:])
            out.append(f'<div class="disclaimer">{_inline(" ".join(quoted))}</div>')
        elif s.startswith("- "):
            flush_para(); bullets.append(s[2:])
        else:
            flush_bullets(); para.append(s)
        i += 1
    flush_para(); flush_bullets()
    return "\n".join(out)


# --- document shell ----------------------------------------------------------
# The stylesheet, theme switch and logo live in theme.py; these names are kept because
# inheritance.py and older callers import them from here.
_STYLE = theme.STYLE


# --- two-pane shell (navigation tree + content) -------------------------------
# Only the index pages use this. The report pages stay single-column and fully
# self-contained on purpose (see the module docstring): a report that is emailed on its
# own must not carry a navigation tree of links to files that did not travel with it.
_SHELL_STYLE = theme.SHELL_STYLE


PANE = "report"   # name of the reading pane the navigation tree's links open in


def _shell(title: str, nav_html: str, body_html: str, pane: bool = False) -> str:
    """A navigation page: tree pane on the left, content on the right.

    With `pane`, the right-hand side is a frame that starts out showing `body_html` and
    then shows whichever report is clicked in the tree. The reports are complete,
    self-contained pages (they have to survive being emailed alone), so opening one used to
    replace the whole window and lose the tree — every report read cost a trip back.
    A named frame needs no script: the starting content rides in `srcdoc`, whose relative
    links resolve against this page exactly as they did when it was inline."""
    from . import config
    if config.DEMO_NOTICE:
        body_html = (f'<div class="disclaimer">{_inline(config.DEMO_NOTICE)}</div>\n'
                     + body_html)
    head = ("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
            f"<title>{_esc(title)}</title>\n<style>{_STYLE}{_SHELL_STYLE}</style>\n"
            f"{theme.SCRIPT}\n</head>\n")
    # The theme switch sits at the foot of the tree on every navigation page.
    nav_html = nav_html.replace("</nav>", theme.SWITCHER + "</nav>", 1)
    if pane:
        start = f"{head}<body>\n<main class=\"content\">\n{body_html}\n</main>\n</body>\n</html>\n"
        return (f"{head}<body class=\"shell\">\n{nav_html}\n"
                f'<iframe class="pane" name="{PANE}" title="Report" srcdoc="{_attr(start)}">'
                "</iframe>\n</body>\n</html>\n")
    return (f"{head}<body class=\"shell\">\n{nav_html}\n<main class=\"content\">\n"
            f"{body_html}\n</main>\n</body>\n</html>\n")


def _wrap(title: str, body_html: str) -> str:
    """A report page: one column, self-contained, with its own theme switch — hidden when
    the report is shown inside a dashboard's reading pane, which has one already."""
    return (
        "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        f"<title>{_esc(title)}</title>\n<style>{_STYLE}</style>\n{theme.SCRIPT}\n</head>\n"
        f"<body>\n<main class=\"report\">\n<div class=\"topbar\">{theme.SWITCHER}</div>\n"
        f"{body_html}\n</main>\n</body>\n</html>\n"
    )


def to_html(md_text: str, title: str) -> str:
    """Full self-contained HTML document for a markdown report."""
    return _wrap(title, _convert_body(md_text))


def ensure_fonts() -> Path:
    """Copy the bundled fonts (and their licences) into reports/html/assets/fonts/ and
    return that folder. Only missing or changed files are written."""
    import shutil
    src = Path(__file__).resolve().parent / "fonts"
    dest = reportpaths.html_root() / "assets" / "fonts"
    if src.is_dir():
        dest.mkdir(parents=True, exist_ok=True)
        for f in src.iterdir():
            target = dest / f.name
            if f.is_file() and (not target.exists() or target.stat().st_size != f.stat().st_size):
                shutil.copyfile(f, target)
    return dest


def with_fonts(html: str, out_path: Path) -> str:
    """Point a page at the bundled fonts. For a page inside reports/html/ the placeholder in
    its stylesheet becomes @font-face rules with the right relative path; for a page written
    anywhere else it is left as the CSS comment it is, and system fonts are used."""
    try:
        out_path.resolve().relative_to(reportpaths.html_root().resolve())
    except ValueError:
        return html
    prefix = os.path.relpath(ensure_fonts(), out_path.resolve().parent).replace(os.sep, "/")
    return html.replace(theme.FONTS_TOKEN, theme.font_faces(prefix))


# --- formats -----------------------------------------------------------------
# A format is a named renderer: it takes one markdown report and returns the bytes of its
# rendering, which land at reports/<name>/<same relative path>.<ext>. Markdown is the source
# of truth and is never a registered format. See docs/report-formats.md.
FORMATS: dict[str, tuple[str, object]] = {}


def register_format(name: str, ext: str, renderer) -> None:
    """Register `renderer(body_md, fields, title) -> str | bytes` as format `name`.

    `body_md` is the report without its front matter, `fields` the front matter as a dict,
    `title` the report's first heading."""
    FORMATS[name] = (ext, renderer)


def render_formats(md_path: Path) -> list[Path]:
    """Render one markdown report into every registered format; return the files written.
    One format failing is logged and does not stop the others."""
    md_path = Path(md_path)
    text = md_path.read_text()
    fields, body = split_front_matter(text)
    m = re.search(r"^# (.+)$", body, re.M)
    title = m.group(1).strip() if m else md_path.stem
    written = []
    for name, (ext, renderer) in FORMATS.items():
        try:
            out = reportpaths.derived(md_path, name, ext)
            data = renderer(body, fields, title)
            out.parent.mkdir(parents=True, exist_ok=True)
            if isinstance(data, bytes):
                out.write_bytes(data)
            else:
                out.write_text(with_fonts(data, out) if name == "html" else data,
                               encoding="utf-8")
            written.append(out)
        except Exception as e:  # noqa: BLE001 — a derived format never blocks a report
            logging.getLogger("genome").warning("%s render failed for %s: %s",
                                                name, md_path.name, e)
    return written


def _render_html(body: str, fields: dict, title: str) -> str:
    return _wrap(title, (lambda h, ml: h.replace("</h1>", "</h1>\n" + ml, 1) if ml else h)(
        _convert_body(body), _meta_line(fields)))


register_format("html", "html", _render_html)


def write_html_for(md_path: Path) -> Path:
    """Render a report `.md` to HTML and return the html path: the same relative path under
    reports/html/ for a report in the reports tree, a sibling file otherwise."""
    md_path = Path(md_path)
    text = md_path.read_text()
    m = re.search(r"^# (.+)$", text, re.M)
    title = m.group(1).strip() if m else md_path.stem
    html_path = reportpaths.html_for(md_path)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.write_text(with_fonts(to_html(text, title), html_path), encoding="utf-8")
    return html_path


# --- per-callset report discovery -------------------------------------------
def current_snapshot() -> str | None:
    """The snapshot whose pages represent "now".

    The newest of what the database holds and what the reports directory holds, because
    those disagree on a machine that received report files from the box that ran the scan.
    Only this snapshot's page may present the undated per-sample reports as its own."""
    from . import reportpaths as rp
    from .snapshot import latest_snapshot
    return max([x for x in [latest_snapshot(), *rp.md_snapshots(), *_html_snapshots()] if x],
               default=None)


def _html_snapshots() -> list[str]:
    """Dated directories under reports/html/ — every snapshot that has pages."""
    from . import reportpaths as rp
    root = rp.html_root()
    return sorted(d.name for d in root.iterdir()
                  if d.is_dir() and d.name not in (rp.GENOME, "assets")) if root.is_dir() else []


# Reports that are genuinely per-snapshot live in the dated directory; everything else is
# a per-sample stage dir with no date in it at all. On an older snapshot's page those are
# NOT that snapshot's reports — an older snapshot's page was listing the PGx and de-novo reports
# of a genome sequenced after that snapshot was taken — so they move under their own heading.
STALE_GROUP = "Current (not from this snapshot)"


def _discover(sample: str, snap: str | None = None,
              undated_group: str | None = None) -> list[tuple[str, str, Path]]:
    """(group, label, md_path) for each existing report of one CALLSET, in display
    order: the snapshot's FULL/panel/family/delta, then the per-sample engines.

    Keyed on the callset id as it appears in samples.tsv (`Adam`, `Adam_vendor`),
    because that is what the report filenames are keyed on. `person_reports()` groups
    these by person.

    `snap` defaults to the latest snapshot but is explicit so an older — or, when this
    machine's database is behind the one that ran the scan, a *newer* — dated report
    directory can still be given its dashboards."""
    from . import config, reportpaths as rp
    from .snapshot import latest_snapshot

    items: list[tuple[str, str, Path]] = []
    snap = snap or latest_snapshot()
    if snap:
        for label, kind in (("Actionable findings (FULL)", "full"),
                            ("Monogenic disease panel", "panel"),
                            ("ACMG-style classification", "acmg")):
            md = rp.snapshot_report(kind, sample, snap)
            if md.exists():
                items.append(("Health", label, md))
        fam_dir = rp.snapshot_md_dir(snap) / rp.FAMILY
        if (fam_dir / "panel_family.md").exists():
            items.append(("Health", "Family reproductive-risk panel",
                          fam_dir / "panel_family.md"))
        deltas = sorted(fam_dir.glob("delta_*.md"))
        if deltas:
            items.append(("What changed", "Latest changes (DELTA)", deltas[-1]))
        # Interpret (GPU regulatory-effect) is snapshot-scoped like DELTA, not per-sample —
        # one report per (model, snapshot), not per person. Prefer the current
        # INTERPRET_MODEL's file; fall back to a glob since both enformer_ and borzoi_
        # reports exist in practice (the model can be switched between runs).
        interp = rp.family_report(f"interpret_{config.INTERPRET_MODEL}", snap)
        if not interp.exists():
            found = sorted(fam_dir.glob("interpret_*.md"))
            interp = found[0] if found else interp
        if interp.exists():
            items.append(("Health", "Regulatory effect (Enformer/Borzoi)", interp))

    for group, label, kind in GENOME_REPORTS:
        md = rp.genome_report(kind, sample)
        if md.exists():
            items.append((undated_group or group, label, md))
    return items


# The per-genome reports, in dashboard order: (group, label, kind). The kind is the file
# name prefix a stage passes to reportpaths.genome_report, so a new stage is listed here once.
GENOME_REPORTS = (
    ("Health", "Pharmacogenomics (drug response)", "pgx"),
    ("Health", "HLA risk (B27, celiac)", "hla"),
    ("Health", "Repeat expansions", "repeats"),
    ("Health", "Structural variants (SV/CNV)", "sv"),
    ("Health", "Mitochondrial variants", "mito"),
    ("Health", "Callability & coverage", "callable"),
    ("Health", "Phasing (cis/trans)", "phase"),
    ("Health", "Trio de-novo variants", "denovo"),
    ("Health", "Segregation & compound het", "segregation"),
    ("Health", "Force-call QC", "forcecall"),
    ("Health", "HLA typing (4-digit alleles)", "hla_typing"),
    ("Wellness", "Traits & wellness", "traits"),
    ("Wellness", "Polygenic risk scores", "prs"),
    ("Ancestry", "Ancestry — haplogroups", "ancestry"),
)


def person_reports(person: str, fam=None, snap: str | None = None,
                   undated_group: str | None = None
                   ) -> list[tuple[str, list[tuple[str, str, Path]]]]:
    """[(callset_id, items)] for every ingested callset of `person`, primary callset first.

    `Adam` and `Adam_vendor` are one person called twice, so their reports belong on
    one page — the bare id's reports as the person's own view, the suffixed ones beside
    them as alternate callsets (`add, never replace`)."""
    from . import pedigree
    fam = fam or pedigree.load()
    ids = fam.callsets(person) or [person]
    return [(cid, _discover(cid, snap, undated_group)) for cid in ids]


def person_report_paths(sample: str) -> list[Path]:
    """All of one CALLSET's report .md files, existing only."""
    return [md for _g, _l, md in _discover(sample)]


# --- navigation tree ---------------------------------------------------------
_GROUPS = ("Health", "Wellness", "Ancestry", "What changed", STALE_GROUP)


def _links(items, out_dir: Path, current: Path | None = None, target: str = "") -> str:
    """<li> links to report HTML, relative to the page being written. `target` names the
    frame the link opens in (the navigation tree aims its links at the reading pane)."""
    li = []
    tgt = f' target="{_attr(target)}"' if target else ""
    for label, md in items:
        href = os.path.relpath(reportpaths.html_for(md), out_dir)
        cls = ' class="here"' if current and md == current else ""
        li.append(f'<li{cls}><a href="{_attr(href)}"{tgt}>{_esc(label)}</a></li>')
    return "".join(li)


EGO_KEY = "genome-ego"   # browser-side: whose perspective relationships are shown from


def perspective_control(fam) -> str:
    """The "Viewing as" picker for the navigation tree: choose whose point of view family
    relationships are labelled from.

    The pages are static files, so the choice lives in the browser (localStorage) and a
    small script fills in the labels: every element carrying `data-rel-of="<person>"`
    receives that person's relationship to the chosen one. With scripting off the picker is
    inert and the pages simply show names, as they did before it existed."""
    picks = ['<li><a href="#" data-ego-pick="">No one (names only)</a></li>']
    picks += [f'<li><a href="#" data-ego-pick="{_attr(pid)}">{_esc(p.display_name)}</a></li>'
              for pid, p in fam.people.items()]
    return ('<details class="persp"><summary>Viewing as: <b data-ego-name>no one</b>'
            f'</summary><ul class="fam">{"".join(picks)}</ul></details>')


_PERSPECTIVE_SCRIPT = """
<script>
(function(){
const REL = __REL__, NAMES = __NAMES__, KEY = __KEY__;
function get(){ try { return localStorage.getItem(KEY); } catch (e) { return null; } }
function set(v){ try { v ? localStorage.setItem(KEY, v) : localStorage.removeItem(KEY); }
                 catch (e) {} }
function apply(){
  let ego = get();
  if (!ego || !(ego in REL)) ego = null;
  document.querySelectorAll('[data-rel-of]').forEach(el => {
    const p = el.dataset.relOf;
    el.textContent = ego ? (p === ego ? 'you' : (REL[ego][p] || '')) : '';
  });
  document.querySelectorAll('[data-ego-name]').forEach(el =>
    el.textContent = ego ? NAMES[ego] : 'no one');
  document.querySelectorAll('[data-ego-pick]').forEach(el =>
    el.parentElement.classList.toggle('here', el.dataset.egoPick === (ego || '')));
  document.querySelectorAll('[data-person]').forEach(el =>
    el.classList.toggle('ego', el.dataset.person === ego));
}
document.querySelectorAll('[data-ego-pick]').forEach(el =>
  el.addEventListener('click', e => {
    e.preventDefault(); set(el.dataset.egoPick); apply();
    const d = el.closest('details'); if (d) d.open = false;
  }));
apply();
})();
</script>
"""


def _js(value) -> str:
    """JSON for inlining in a <script>: `</` cannot be allowed to close the element."""
    import json
    return json.dumps(value).replace("</", "<\\/")


def perspective_script(fam) -> str:
    """The script behind `perspective_control`, with every pairwise relationship computed
    here (people x people is tiny) so the browser only looks labels up."""
    rel = {ego: {other: r for other, r in fam.others(ego)} for ego in fam.people}
    names = {pid: p.display_name for pid, p in fam.people.items()}
    return (_PERSPECTIVE_SCRIPT.replace("__REL__", _js(rel))
            .replace("__NAMES__", _js(names)).replace("__KEY__", _js(EGO_KEY)))


def snapshot_switcher(out_dir: Path, page: str) -> str:
    """A collapsed "Snapshots" section for the navigation tree: the same page in every other
    dated report directory beside `out_dir`, newest first.

    Each snapshot is one pinned ClinVar release, so this is how a reader moves between runs
    — the same person, as the database stood on another date. Built from what is on disk
    when the page is written, and only where that directory actually has this page (someone
    sequenced later has no page in an earlier snapshot). A static page cannot learn about
    snapshots made after it, so every dashboard build refreshes the older directories' pages
    too (`refresh_other_snapshots`)."""
    from .snapshot import load_manifest
    parent = out_dir.parent
    if not parent.is_dir():
        return ""
    names = sorted((d.name for d in parent.iterdir()
                    if d.is_dir() and (d / "index.html").exists() or d == out_dir),
                   reverse=True)
    if len(names) < 2:
        return ""
    rows = []
    for i, name in enumerate(names):
        try:
            clinvar = load_manifest(name).get("clinvar_date")
        except Exception:  # noqa: BLE001 — a report dir synced without its manifest
            clinvar = None
        note = " · ".join(x for x in (f"ClinVar {clinvar}" if clinvar else "",
                                      "latest" if i == 0 else "") if x)
        rel = f'<span class="rel">{_esc(note)}</span>' if note else ""
        if name == out_dir.name:
            rows.append(f'<li class="here">{_esc(name)}{rel}</li>')
        elif (parent / name / page).exists():
            rows.append(f'<li><a href="../{_attr(name)}/{_attr(page)}">{_esc(name)}</a>'
                        f"{rel}</li>")
        else:
            rows.append(f'<li><span class="unseq">{_esc(name)}</span>'
                        '<span class="rel">not in this snapshot</span></li>')
    return (f'<details><summary>Snapshots ({len(names)})</summary>'
            f'<ul class="fam">{"".join(rows)}</ul></details>')


def refresh_other_snapshots(snap: str | None) -> int:
    """Rewrite the navigation pages of every OTHER dated directory, so their snapshot lists
    include the one just built. Only pages that already exist are rewritten — this never
    gives an old snapshot a page for someone it did not have. Returns directories touched."""
    root = reportpaths.html_root()
    if not root.is_dir():
        return 0
    done = 0
    for d in sorted(p for p in root.iterdir() if p.is_dir() and p.name != snap):
        if not (d / "index.html").exists():
            continue
        try:
            for page in sorted(d.glob("index_*.html")):
                build_index(page.stem[len("index_"):], out_dir=d, snap=d.name)
            build_snapshot_index(d.name, out_dir=d)
            if (d / "inheritance.html").exists():
                from . import inheritance
                inheritance.build_page(d.name, out_dir=d)
            done += 1
        except (Exception, SystemExit) as e:  # noqa: BLE001 — never fail the scan over nav
            log.warning("snapshot list not refreshed for %s: %s", d.name, e)
    return done


def _tree(person: str, out_dir: Path, fam, snap: str | None = None,
          undated_group: str | None = None) -> str:
    """The left navigation pane: expandable <details> sections, no JavaScript.

    Report links open in the reading pane beside the tree (`target=PANE`), so the tree
    stays put while reports are read one after another. Links to another person, or to
    the inheritance page, are whole pages with their own tree and replace this one.

    Plain <details>/<summary> rather than a scripted tree — it expands and collapses with
    the browser's own machinery, so the page works identically over file:// and behind
    `run.py serve`, and degrades to a nested list in anything that cannot style it."""
    display = fam.people[person].display_name if person in fam.people else person
    out = [f'<nav class="tree" aria-label="Reports"><p class="brand">{theme.LOGO}'
           f'<a href="index.html">Genome reports</a></p>',
           f'<p class="who">{_esc(display)}</p>']

    by_callset = person_reports(person, fam, snap, undated_group)
    primary = by_callset[0] if by_callset else (person, [])
    others = by_callset[1:]

    for gname in _GROUPS:
        links = [(label, md) for g, label, md in primary[1] if g == gname]
        if not links:
            continue
        out.append(f'<details open><summary>{_esc(gname)}</summary>'
                   f'<ul>{_links(links, out_dir, target=PANE)}</ul></details>')

    # Alternate callsets: same person, called by another pipeline or against another
    # reference. Nested one level deeper so they read as a subordinate view, never as a
    # separate person (which is exactly what the old family.tsv turned them into).
    if others:
        inner = []
        for cid, items in others:
            links = [(label, md) for _g, label, md in items]
            if not links:
                continue
            inner.append(f'<details><summary>{_esc(_callset_label(cid, person))}</summary>'
                         f'<ul>{_links(links, out_dir, target=PANE)}</ul></details>')
        if inner:
            out.append('<details><summary>Other callsets</summary>'
                       f'<div class="sub">{"".join(inner)}</div></details>')

    fam_li = []
    for other_id, rel in fam.others(person):
        dn = _esc(fam.people[other_id].display_name)
        if fam.is_sequenced(other_id):
            fam_li.append(f'<li><a href="index_{_attr(other_id)}.html">{dn}</a>'
                          f'<span class="rel">{_esc(rel)}</span></li>')
        else:
            fam_li.append(f'<li><span class="unseq">{dn}</span>'
                          f'<span class="rel">{_esc(rel)}, not sequenced</span></li>')
    if fam_li:
        out.append(f'<details open><summary>Family</summary><ul class="fam">'
                   f'{"".join(fam_li)}</ul></details>')

    # Only when the file exists: the inheritance page is built from the variant rows, and
    # a machine holding report files synced from elsewhere has the reports without the
    # rows. Linking it unconditionally put a 404 in the nav of every page in that dir.
    if (out_dir / "inheritance.html").exists():
        out.append('<details open><summary>Family-wide</summary><ul>'
                   '<li><a href="inheritance.html">Variant inheritance tree</a></li>'
                   "</ul></details>")
    out.append(snapshot_switcher(out_dir, f"index_{person}.html"))
    out.append("</nav>")
    # This page IS this person's view (the Family list above is relative to them), so
    # opening it makes them the perspective the family-wide pages use.
    if person in fam.people:
        out.append(f"<script>try{{localStorage.setItem({_js(EGO_KEY)},{_js(person)})}}"
                   "catch(e){}</script>")
    return "".join(out)


_SUFFIX_NAMES = {"vendor": "vendor pipeline", "t2t": "T2T-CHM13 (lifted to GRCh38)"}


def _callset_label(callset: str, person: str) -> str:
    """Human label for an alternate callset: `Adam_vendor` -> 'vendor pipeline'."""
    suffix = callset[len(person):].lstrip("_")
    return _SUFFIX_NAMES.get(suffix, suffix or callset)


# --- per-person dashboard ----------------------------------------------------
def build_index(person: str, out_dir: Path | None = None,
                snap: str | None = None) -> Path:
    """One landing page per PERSON, written beside that snapshot's reports.

    Lives in `reports/html/<snapshot>/` rather than at the html root: the page is a view of
    one snapshot's report set, so it belongs with the files it links, and the root then
    holds only dated directories plus the landing `index.html`. Every link is relative, so
    the whole dated directory stays portable if copied or archived."""
    from . import config
    from . import pedigree
    from .snapshot import latest_snapshot

    fam = pedigree.load()
    person = fam.person_of(person)
    snap = snap or latest_snapshot()
    out_dir = out_dir or (reportpaths.html_root() / snap if snap
                          else reportpaths.html_root())
    out = out_dir / f"index_{person}.html"
    out.parent.mkdir(parents=True, exist_ok=True)

    p = fam.people.get(person)
    title = p.display_name if p else person
    # Only the current snapshot may show the undated per-sample reports as its own.
    undated_group = None if snap == current_snapshot() else STALE_GROUP
    by_callset = person_reports(person, fam, snap, undated_group)
    n = sum(len(items) for _c, items in by_callset)

    from . import dashboard
    from .snapshot import load_manifest
    try:
        clinvar = load_manifest(snap).get("clinvar_date") if snap else None
    except Exception:  # noqa: BLE001 — a report dir synced without its manifest
        clinvar = None
    meta = " · ".join(x for x in (
        f"Snapshot {snap}" if snap else "", f"ClinVar {clinvar}" if clinvar else "",
        f"{n} report(s)", f"{len(by_callset)} callset(s)") if x)
    body = [f'<p class="meta">{_esc(meta)}</p>', f"<h1>{_esc(title)}</h1>"]
    if snap and undated_group:
        body.append('<p class="disclaimer">This is an <b>older</b> snapshot. Only the '
                    "FULL / panel / ACMG / DELTA reports below are from it; the "
                    f"per-sample reports under <i>{_esc(STALE_GROUP)}</i> have no "
                    "snapshot of their own and show today's results.</p>")
    primary = by_callset[0][0] if by_callset else person
    try:
        body.append(dashboard.overview(person, primary, snap, fam, out_dir))
    except (Exception, SystemExit) as e:  # noqa: BLE001 — the links below still work
        log.warning("dashboard summary skipped for %s: %s", person, e)
    for cid, items in by_callset:
        if not items:
            continue
        which = ("All reports" if cid == person
                 else f"Alternate callset — {_callset_label(cid, person)}")
        body.append(f"<h2>{_esc(which)} <code>{_esc(cid)}</code></h2>")
        for gname in _GROUPS:
            links = [(label, md) for g, label, md in items if g == gname]
            if links:
                body.append(f"<h3>{_esc(gname)}</h3><ul>{_links(links, out_dir)}</ul>")
    if not n:
        body.append("<p>No reports found yet for this person.</p>")
    body.append('<p class="foot">Research-grade, not a clinical diagnosis. Confirm any '
                "actionable finding with a qualified clinician.</p>")

    out.write_text(with_fonts(
        _shell(f"Genome reports — {title}", _tree(person, out_dir, fam, snap, undated_group),
               "\n".join(body), pane=True), out), encoding="utf-8")
    log.info("Dashboard: %s (%d report link(s))", out, n)
    return out


def _roster(fam, base: str) -> str:
    """<li> rows for the whole family: a link per sequenced person, a plain name for the
    rest. `base` prefixes the dated directory when linking in from the reports root."""
    rows = []
    for pid, p in fam.people.items():
        dn = _esc(p.display_name)
        if fam.is_sequenced(pid):
            cs = fam.callsets(pid)
            extra = f'<span class="rel">{len(cs)} callsets</span>' if len(cs) > 1 else ""
            rows.append(f'<li><a href="{_attr(base)}index_{_attr(pid)}.html">{dn}</a>'
                        f'<span class="rel" data-rel-of="{_attr(pid)}"></span>{extra}</li>')
        else:
            rows.append(f'<li><span class="unseq">{dn}</span>'
                        f'<span class="rel" data-rel-of="{_attr(pid)}"></span>'
                        f'<span class="rel">not sequenced</span></li>')
    return "".join(rows)


def _landing_page(fam, snap: str | None, base: str, title: str, subtitle: str,
                  older: list[str] | None = None, has_tree: bool = True,
                  switcher: str = "", out_dir: Path | None = None) -> str:
    from . import dashboard
    rows = _roster(fam, base)
    nav = ['<nav class="tree" aria-label="Reports">',
           f'<p class="brand">{theme.LOGO}{_esc(title)}</p>',
           perspective_control(fam),
           f'<details open><summary>Family</summary><ul class="fam">{rows}</ul></details>']
    if snap and has_tree:
        nav.append('<details open><summary>Family-wide</summary><ul>'
                   f'<li><a href="{_attr(base)}inheritance.html">Variant inheritance tree</a>'
                   "</li></ul></details>")
    nav.append(switcher)
    nav.append("</nav>")

    # Hrefs to reports are relative to the page: the root landing page sits one level above
    # the snapshot directory its links point into (`base`), the snapshot's own page inside it.
    if out_dir is None:
        out_dir = reportpaths.html_root() if base else (
            reportpaths.html_root() / snap if snap else reportpaths.html_root())
    headline, sentence = dashboard.hero(fam)
    try:
        changes = dashboard.changes_card(snap, base, out_dir)
    except (Exception, SystemExit):  # noqa: BLE001 — the people list still works
        changes = ""
    body = [f'<p class="meta">{subtitle}</p>',
            f'<div class="grid g2" style="align-items:end">'
            f'<div><h1 class="hero">{headline}</h1><p class="lede">{_esc(sentence)}</p></div>'
            f"{changes}</div>",
            "<h2>People</h2>"]
    try:
        body.append(dashboard.people_cards(fam, snap, base))
    except (Exception, SystemExit):  # noqa: BLE001
        body.append(f'<ul class="fam">{rows}</ul>')
    if older:
        body.append('<h2>Other snapshots</h2><div class="snaps">'
                    + (f'<span class="cur">{_esc(snap)}</span>' if snap else "")
                    + "".join(f'<i></i><a href="{_attr(d)}/index.html">{_esc(d)}</a>'
                              for d in older) + "</div>")
    body.append('<p class="foot">Research-grade, not a clinical diagnosis. Confirm any '
                "actionable finding with a qualified clinician.</p>")
    return _shell(title, "".join(nav), "\n".join(body) + perspective_script(fam))


def build_landing(out_path: Path | None = None) -> Path:
    """`reports/html/index.html` — the one stable entry point.

    The per-person pages move with each snapshot; this does not, so it is what a bookmark
    or a LAN link should point at. It opens the newest dated directory that actually has
    pages, and lists the rest.

    Newest-with-pages rather than `latest_snapshot()`, because those two can disagree:
    this is a *reports* entry point, and a machine whose database is behind the one that
    ran the scan (reports synced from it, snapshots not) would otherwise send every link
    to a stale directory while a newer one sits beside it."""
    from . import config
    from . import pedigree
    from .snapshot import latest_snapshot

    fam = pedigree.load()
    root = reportpaths.html_root()
    out = out_path or (root / "index.html")
    out.parent.mkdir(parents=True, exist_ok=True)
    dated = sorted((d.name for d in root.iterdir()
                    if d.is_dir() and (d / "index.html").exists()), reverse=True)
    db_snap = latest_snapshot()
    snap = dated[0] if dated else db_snap
    older = [d for d in dated if d != snap]
    subtitle = ("No snapshot yet — run <code>scan</code> first." if not snap
                else f"Showing <code>{_esc(snap)}</code>."
                + (f" This machine's database is at <code>{_esc(db_snap)}</code>, so the "
                   "newer report set was built elsewhere." if db_snap and db_snap != snap
                   else ""))
    has_tree = bool(snap) and (root / snap / "inheritance.html").exists()
    out.write_text(with_fonts(_landing_page(fam, snap, f"{snap}/" if snap else "",
                                            "Genome reports", subtitle, older, has_tree),
                              out), encoding="utf-8")
    log.info("Landing page: %s", out)
    return out


def build_snapshot_index(snap: str | None = None, out_dir: Path | None = None) -> Path:
    """`reports/html/<snapshot>/index.html` — the dated directory's own front page.

    What the per-person pages' "Genome reports" link goes back to, and what makes a copied
    dated directory self-sufficient: its links never leave the directory."""
    from . import config
    from . import pedigree
    from .snapshot import latest_snapshot

    fam = pedigree.load()
    snap = snap or latest_snapshot()
    out_dir = out_dir or (reportpaths.html_root() / snap if snap
                          else reportpaths.html_root())
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "index.html"
    subtitle = (f"Report set for snapshot <code>{_esc(snap)}</code>." if snap
                else "No snapshot yet — run <code>scan</code> first.")
    out.write_text(with_fonts(
        _landing_page(fam, snap, "", "Genome reports", subtitle,
                      has_tree=(out_dir / "inheritance.html").exists(),
                      switcher=snapshot_switcher(out_dir, "index.html")), out),
        encoding="utf-8")
    log.info("Snapshot index: %s", out)
    return out


# Built-in formats beyond HTML. Imported last, because a format module imports names from
# this one; importing it here (not from run.py) means the format exists for every caller
# that writes a report — the CLI, the tests, a notebook.
from . import formats_json  # noqa: E402,F401
