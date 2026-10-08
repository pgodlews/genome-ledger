"""JSON rendering of reports — a registered format beside HTML.

The markdown in reports/md/ is the source of truth; this renders each report to
reports/json/<same relative path>.json with its front-matter fields, its title, and
its content in a structured form: an ordered list of sections, each holding its text
lines and its tables as lists of row objects keyed by the column headers — a
retriever never has to parse a markdown table itself. Placement, consent-driven
deletion and failure isolation come from render.render_formats; this module only
answers "what is this report, as JSON". Importing it registers the format, and
render.py imports it, so the format exists wherever reports are written.
See docs/report-formats.md, "Adding a format".

    {
      "title": "...", <every front-matter field>,
      "summary": ["lines before the first heading"],
      "tables":  [tables before the first heading],
      "sections": [{"heading": "...", "level": 2 | 3, "text": [...], "tables": [...]}]
    }

Sections are a list, not an object keyed by heading: a report may repeat a heading (the
de-novo report has "How the candidates were filtered" twice), and their order is content.
`level` 3 is a subsection of the nearest level-2 section before it; a report may also open
with level-3 sections and no level-2 parent at all (the polygenic report groups its tables
that way), which is why a table is structured wherever it stands — under a `##`, under a
`###`, or before any heading.
"""

from __future__ import annotations

import json

from .render import _cells, _is_separator, register_format


def _table_object(rows: list[str]) -> dict:
    """One markdown table as {"headers": [...], "rows": [row objects]} — each row an
    object keyed by the table's column headers. The |---| separator row is markdown
    table syntax, not data. Cell text is kept verbatim: markdown emphasis stays
    markdown (the JSON carries the source of truth's content, not its HTML)."""
    headers = _cells(rows[0])
    body_start = 2 if len(rows) > 1 and _is_separator(_cells(rows[1])) else 1
    out = []
    for r in rows[body_start:]:
        cells = _cells(r)
        cells += [""] * (len(headers) - len(cells))   # a short row pads, never drops a key
        out.append(dict(zip(headers, cells)))
    return {"headers": headers, "rows": out}


def to_json(body: str, fields: dict, title: str) -> str:
    """body: the markdown without its front matter. fields: the front matter.
    title: the report's first heading. Return str or bytes."""
    top = {"text": [], "tables": []}      # everything before the first heading
    sections: list[dict] = []
    current = top
    lines = body.splitlines()
    i = 0
    while i < len(lines):
        line, s = lines[i], lines[i].strip()
        if s.startswith("|") and s.endswith("|"):
            tbl = [line]
            i += 1
            while i < len(lines) and lines[i].strip().startswith("|"):
                tbl.append(lines[i])
                i += 1
            current["tables"].append(_table_object(tbl))
            continue
        if s.startswith(("## ", "### ")):
            level = 3 if s.startswith("### ") else 2
            current = {"heading": s[level + 1:].strip(), "level": level,
                       "text": [], "tables": []}
            sections.append(current)
        elif s.startswith("# "):
            pass                    # the title itself — carried as its own key
        elif s:
            current["text"].append(s)
        i += 1
    doc = {**fields, "title": title, "summary": top["text"], "tables": top["tables"],
           "sections": sections}
    return json.dumps(doc, ensure_ascii=False, indent=2) + "\n"


register_format("json", "json", to_json)
