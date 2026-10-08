# Reports: layout, the file contract, and adding a format

Every report the pipeline writes is a markdown file. Markdown is the source of truth;
anything else (HTML today) is rendered from it. This page describes where reports live, what
every report file looks like, and how to add another output format.

## Layout

One tree per format, the same relative path in each (`pipeline/reportpaths.py` owns this):

```
$GENOMES_ROOT/reports/
  md/                           the source of truth
    <snapshot>/                   the LATEST snapshot only
      _family/                      delta_<prev>_to_<snapshot>.md, panel_family.md,
                                    inheritance.md, interpret_<model>.md, validation_….md
      <person>/                     full_<callset>.md, panel_<callset>.md, acmg_<callset>.md
    genome/                       reports that depend on the genome, not on a ClinVar release
      <person>/                     pgx_… prs_… hla_… hla_typing_… repeats_… sv_… mito_…
                                    callable_… phase_… denovo_… segregation_… forcecall_…
                                    traits_… ancestry_…   (each `<kind>_<callset>.md`)
  md-history/                   older snapshots' markdown, same shape as md/<snapshot>/
    <snapshot>/…
  html/                         rendered reports, plus the navigation pages
    index.html                    the stable entry point
    <snapshot>/                   index.html, index_<person>.html, inheritance.html,
                                  _family/…, <person>/…
    genome/<person>/…
  json/                         every report as structured JSON, same relative paths
```

Three rules explain the shape:

- **`md/` holds only the latest snapshot.** A snapshot is one pinned ClinVar release, so an
  older snapshot's report states classifications that have since been superseded. `md/` is
  meant to be handed whole to a search or retrieval tool, and such a tool will happily
  quote "uncertain significance" from last quarter for a variant that is pathogenic today.
  When a newer snapshot is built, the previous one's directory moves to `md-history/`.
  HTML keeps every snapshot, because there a person picks which one to read.
- **Reports are filed per person, named per callset.** `Adam` and `Adam_vendor` are one
  person called twice; both callsets' reports sit in `Adam/`. Handing a tool `md/genome/Adam`
  and `md/<snapshot>/Adam` gives it one person and nobody else.
- **Stage directories hold working data only.** `pgx/`, `prs/`, `repeats/` and the rest keep
  VCFs, tool output and caches — no reports. Only `reports/html/` is web-servable.

Per-genome reports are not dated: they change when the genome or the tool changes, not when
ClinVar does. On an older snapshot's dashboard they are grouped under *Current (not from
this snapshot)* for that reason.

## The file contract

Every report in the tree has the same shape, enforced in one place (`util.write_report`):

```markdown
---
title: "Polygenic risk scores — Adam"
report: prs                 # the kind: full, panel, acmg, delta, pgx, prs, hla, …
scope: genome               # genome | snapshot | family
person: Adam                # absent for family-scope reports
callset: Adam               # which of the person's callsets
snapshot: 2026-09-06        # snapshot and family scope only
clinvar_date: 2026-09-06    # the ClinVar release that snapshot pinned
generated: 2026-09-07
---

# Polygenic risk scores — Adam

- summary bullets: the headline numbers

## <findings sections, as many as the report needs>

## Caveats
- what this result cannot tell you
```

- **Front matter** is derived from where the file is filed, so a writer cannot get it wrong.
  It is flat `key: value` YAML. `demo: true` appears on output of `run.py demo`.
- **Title** is the first line of the body: `# <Report name> — <callset>`.
- **Summary** comes first, as bullets directly under the title.
- **One closing section**, always last, under one of three headings:
  `## Caveats`, `## Method & caveats` (when it also explains how the result was derived), or
  `## Assumptions & Known Limitations` (the auto-generated block of pinned versions). A
  report written without one gets a minimal `## Caveats` appended, and `## Notes` is filed
  as `## Caveats`. A reader, or a tool splitting on headings, always finds a result's limits
  in the same place.
- **Notices** (withheld traits, demo data) are blockquotes: lines starting with `> `.

## Adding a report

1. Ask `reportpaths` for the path — never build one by hand:
   `reportpaths.genome_report("<kind>", sample)`,
   `reportpaths.snapshot_report("<kind>", sample, snapshot)`, or
   `reportpaths.family_report("<name>", snapshot)`.
2. Write it with `util.write_report(path, lines)`. That adds the front matter, enforces the
   closing section, and renders every registered format.
3. For a per-genome report, add one line to `GENOME_REPORTS` in `pipeline/render.py` so the
   dashboards list it. A test fails if a kind is written but not listed, or listed but
   never written.

## Adding a format

A format is a named renderer. It receives one report and returns its rendering; the
pipeline decides where the file goes — `reports/<format>/<same relative path>.<ext>` — and
runs it for every report, on every write and on `run.py render`.

```python
# pipeline/formats_text.py
from .render import register_format

def to_text(body: str, fields: dict, title: str) -> str:
    """body: the markdown without its front matter. fields: the front matter, as a dict of
    strings. title: the report's first heading. Return str or bytes."""
    keep = [ln for ln in body.splitlines() if not ln.startswith("# ")]
    header = f"{title}\n{fields.get('person', 'family')} · {fields.get('generated', '')}\n"
    return header + "\n".join(keep).strip() + "\n"

register_format("text", "txt", to_text)
```

Three things the example relies on, each of which has bitten someone:

- **Registration is an import.** A format exists once its module has been imported. Add
  `from . import formats_text` beside the other built-in format at the bottom of
  `pipeline/render.py`; then the CLI, the tests and any other caller that writes a report
  all have it. Importing it only from `run.py` leaves the tests without it.
- **`title` is the clean one.** `fields` also carries a `title`; build your output so the
  argument wins (`{**fields, "title": title}`, not the other way round).
- **`body` starts with the title line and may have content before the first `##`.** The
  summary bullets, a notice, even a table can sit there, and some reports group their
  tables under `###` headings with no `##` above them. Handle a table wherever it stands.

### Built-in: JSON

`pipeline/formats_json.py` writes `reports/json/<same path>.json` for every report, so a
tool never has to parse markdown:

```json
{
  "title": "Polygenic risk scores — Adam",
  "report": "prs", "scope": "genome", "person": "Adam", "callset": "Adam",
  "generated": "2026-09-07",
  "summary": ["lines before the first heading: the headline bullets and any notice"],
  "tables": ["tables before the first heading, same shape as below"],
  "sections": [
    {"heading": "Preventive treatment exists", "level": 3,
     "text": ["non-table lines, markdown markers kept"],
     "tables": [{"headers": ["Trait", "PGS ID", "EUR percentile"],
                 "rows": [{"Trait": "…", "PGS ID": "…", "EUR percentile": "91.0%"}]}]}
  ]
}
```

- Front-matter fields are at the top level, as strings (`"demo": "true"`).
- `sections` is an ordered list, because a report may repeat a heading and the order is
  content. `level` is 2 for `##` and 3 for `###`; a level-3 section belongs to the nearest
  level-2 section before it, if there is one.
- Every table is a list of row objects keyed by its column headers, wherever it stands.
  Cell text is verbatim markdown (`**Pathogenic**`); the `|---|` separator row is syntax,
  not data.

What you get for free, and what to keep in mind:

- **Placement and mirroring.** `reports/json/genome/Adam/prs_Adam.json` appears beside the
  other trees with no path code of your own.
- **Consent.** When a person withdraws consent, the report is deleted from every registered
  format's tree, not only from `md/` — a rendering left behind keeps disclosing.
- **Failure isolation.** A renderer that raises is logged and skipped; the markdown and the
  other formats are still written.
- **Render from the markdown, not from the database.** A format that re-queries the data
  can disagree with the report it sits beside. If a format needs something the markdown
  does not carry, add it to the report (or its front matter) first.
- **Old layout.** `run.py migrate-reports` treats every registered format's folder as part
  of the new layout, so a format needs no entry there.
- **History.** Formats mirror `md-history/` as well, into the same `<format>/<snapshot>/`
  path. If a format is meant for retrieval, skip reports whose `fields["snapshot"]` is not
  the latest, or point the tool at `md/` instead.
- **Navigation pages** (`index*.html`, the interactive inheritance tree) belong to the HTML
  format only. They are built by `render`, not from a markdown file, and other formats do
  not need an equivalent.
- **Binary formats** (PDF) return `bytes`. Add any new dependency as an optional one, so
  the core install stays light.

## Moving from the old layout

Before this layout, per-snapshot reports sat in `reports/<snapshot>/` beside their HTML,
and each per-genome stage kept its report in its own directory. One command re-files them:

```bash
uv run python run.py migrate-reports            # dry run: lists every move
uv run python run.py migrate-reports --apply    # move the markdown, delete the old HTML
uv run python run.py render                     # rebuild reports/html/
```

Markdown is moved, never regenerated — several reports cost hours of compute, and a stage
decides whether to re-run by whether its report exists.
