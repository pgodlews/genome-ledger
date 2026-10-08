"""Family-wide variant inheritance: who carries which variant, and where it came from.

One page for the whole family rather than one per person. For each ClinVar-actionable
variant carried by anyone, it answers the question a per-person report cannot: *whose* is
it — did it come down the father's line, the mother's line, both, or neither (a de-novo
candidate)?

Two deliberate choices about the evidence:

* **Carrier status is per PERSON, not per callset.** `Adam` and `Adam_vendor` are one
  genome called twice, so a variant called by either is a variant this person carries — the
  same reading `denovo.parent_callsets()` and `segregation._parent_calls()` already use.
  Callset disagreement is not hidden: it is reported per variant, because a site only one
  caller sees is exactly the site worth doubting.
* **"De-novo" here means "neither parent's callsets carry it"**, which is the *candidate*
  definition, not the confirmed one. `pipeline/denovo.py` goes further and force-calls the
  parents plus pileups their raw reads; its report is the authority and this page links to
  it. Labelling the weaker inference as if it were the stronger one would be the single
  most misleading thing this page could do.
"""

from __future__ import annotations

from pathlib import Path

from . import reportpaths, theme
from . import allele_balance, clinsig, config
from .load import TABLE, connect
from .snapshot import latest_snapshot
from .util import log

# Where a person's copy of a variant came from. Ordered by how much the page can actually
# claim: a resolved parental origin, then the two states that are inferences about absence.
FROM_FATHER, FROM_MOTHER, FROM_BOTH = "father", "mother", "both"
DE_NOVO, UNRESOLVED, FOUNDER = "de_novo", "unresolved", "founder"

ORIGIN_LABEL = {
    FROM_FATHER: "from father",
    FROM_MOTHER: "from mother",
    FROM_BOTH: "from both parents",
    DE_NOVO: "de-novo candidate",
    UNRESOLVED: "unresolved",
    FOUNDER: "origin not traceable",
}


# --- evidence ----------------------------------------------------------------
def loaded(snap: str) -> bool:
    """Does this database hold variant rows for `snap`?

    Reports can be copied between machines while the database that produced them stays
    behind, so a dated report directory is not evidence that the rows are here. Without
    this check `collect()` returns nothing and the page renders as "no pathogenic variant
    in this family" — an empty result that looks exactly like a clean one."""
    con = connect()
    try:
        return bool(con.execute(
            f"SELECT 1 FROM {TABLE} WHERE snapshot_id = ? LIMIT 1", [snap]).fetchone())
    finally:
        con.close()


def _rows(con, snap: str):
    """Every actionable-significance call in the snapshot, one row per (callset, variant)."""
    return con.execute(f"""
        SELECT sample_id, chrom, pos, ref, alt, gene, consequence,
               clinvar_sig, clinvar_disease, gnomad_af, zygosity
        FROM {TABLE}
        WHERE snapshot_id = ? AND {clinsig.sql_actionable()}
        ORDER BY chrom, pos, ref, alt, sample_id
    """, [snap]).fetchall()


def collect(snap: str, fam) -> list[dict]:
    """[{key, gene, sig, disease, af, csq, carriers}] for each actionable variant.

    `carriers` maps a PERSON id to {"zyg": "HOM"|"HET", "callsets": {callset: zyg},
    "low_ab": (alt, total) | None} — a person is a carrier when any of their callsets calls
    them one, and the per-callset detail rides along so the page can flag single-caller
    sites. `low_ab` is the weakest below-threshold allele balance among those calls (see
    allele_balance): a low-fraction het is what sequencer slippage and paralog reads look
    like, and without it an artifact in a child reads as a confident de-novo."""
    con = connect()
    try:
        rows = _rows(con, snap)
    finally:
        con.close()

    keys_by_sid: dict[str, list] = {}
    for sid, chrom, pos, ref, alt, *_ in rows:
        keys_by_sid.setdefault(sid, []).append((chrom, int(pos), ref, alt))
    ab = {sid: allele_balance.fractions(sid, keys) for sid, keys in keys_by_sid.items()}

    by_key: dict[tuple, dict] = {}
    for sid, chrom, pos, ref, alt, gene, csq, sig, disease, af, zyg in rows:
        key = (chrom, int(pos), ref, alt)
        v = by_key.setdefault(key, {
            "key": key, "label": f"{chrom}:{pos} {ref}>{alt}", "gene": gene or "?",
            "sig": sig or "", "disease": disease or "", "af": af,
            "csq": (csq or "").split("&")[0], "carriers": {},
        })
        person = fam.person_of(sid)
        z = (zyg or "").upper()
        c = v["carriers"].setdefault(person, {"zyg": z, "callsets": {}, "low_ab": None})
        c["callsets"][sid] = z
        frac = ab.get(sid, {}).get(key)
        if allele_balance.is_low(z, frac) and (
                c["low_ab"] is None or frac[0] / frac[1] < c["low_ab"][0] / c["low_ab"][1]):
            c["low_ab"] = frac
        # HOM outranks HET: one caller seeing both alleles is the stronger observation, and
        # a person read as HET when they are HOM is read as a carrier when they are affected.
        if z == "HOM":
            c["zyg"] = "HOM"
    return sorted(by_key.values(), key=lambda v: (v["gene"], v["key"]))


def origin(variant: dict, person: str, fam) -> str:
    """Where `person`'s copy of `variant` came from, as far as the callsets can say."""
    carriers = variant["carriers"]
    if person not in carriers:
        return ""
    p = fam.people.get(person)
    father, mother = (p.father, p.mother) if p else (None, None)
    if not father or not mother:
        return FOUNDER
    # An absent parent callset cannot testify. Without both parents sequenced, "neither
    # parent carries it" is a statement about the data, not about the family.
    if not (fam.is_sequenced(father) and fam.is_sequenced(mother)):
        return UNRESOLVED
    f, m = father in carriers, mother in carriers
    if f and m:
        return FROM_BOTH
    if f:
        return FROM_FATHER
    if m:
        return FROM_MOTHER
    return DE_NOVO


def discordant(variant: dict, fam) -> list[str]:
    """People whose callsets disagree about this variant — called by some, not by all.

    A site one caller sees and another does not is the page's own confidence signal, and it
    is the reason a de-novo candidate most often evaporates under `pipeline/denovo.py`."""
    out = []
    for person, c in variant["carriers"].items():
        if len(c["callsets"]) < len(fam.callsets(person)):
            out.append(person)
    return out


# --- tree layout --------------------------------------------------------------
def generations(fam) -> dict[str, int]:
    """person -> row in the drawing. Children sit one row below their parents; a founder
    who married in is pulled down to their partner's row so couples draw side by side."""
    gen: dict[str, int] = {}

    def depth(pid: str, seen: frozenset = frozenset()) -> int:
        if pid in gen:
            return gen[pid]
        if pid in seen:                       # cycle guard; validate() reports the cycle
            return 0
        parents = fam.parents(pid)
        return 1 + max((depth(x, seen | {pid}) for x in parents), default=-1)

    for pid in fam.people:
        gen[pid] = depth(pid)
    for pid in fam.people:                    # married-in founders join their partner's row
        partner = fam.partner(pid)
        if not fam.parents(pid) and partner and gen.get(partner, 0) > gen[pid]:
            gen[pid] = gen[partner]
    return gen


def layout(fam) -> tuple[dict[str, tuple[int, int]], int, int]:
    """person -> (x, y) centre in SVG units, plus the drawing's width and height."""
    gen = generations(fam)
    rows: dict[int, list[str]] = {}
    for pid in fam.people:                    # file order, then partners pulled adjacent
        rows.setdefault(gen[pid], []).append(pid)
    for g, members in rows.items():
        ordered: list[str] = []
        for pid in members:
            if pid in ordered:
                continue
            ordered.append(pid)
            partner = fam.partner(pid)
            if partner in members and partner not in ordered:
                ordered.append(partner)
        rows[g] = ordered

    widest = max((len(m) for m in rows.values()), default=1)
    width = max(widest * (NODE_W + NODE_GAP) + NODE_GAP, 520)
    pos: dict[str, tuple[int, int]] = {}
    for g, members in rows.items():
        span = len(members) * (NODE_W + NODE_GAP) - NODE_GAP
        x0 = (width - span) // 2
        for i, pid in enumerate(members):
            pos[pid] = (x0 + i * (NODE_W + NODE_GAP) + NODE_W // 2,
                        ROW_TOP + g * ROW_H + NODE_H // 2)
    height = ROW_TOP + (max(rows, default=0) + 1) * ROW_H
    return pos, width, height


NODE_W, NODE_H, NODE_GAP = 132, 52, 40
ROW_H, ROW_TOP = 132, 24


# --- drawing ------------------------------------------------------------------
def _svg(fam, pos, width, height) -> str:
    """The family tree itself: couple bars, descent lines, one node per person.

    Every element that changes with the selected variant carries a stable id
    (`n-<person>`, `e-<child>`) so the page's script can restyle it without redrawing —
    the geometry is computed once, here, and never in the browser."""
    from .render import _attr, _esc

    parts = [f'<svg id="tree" viewBox="0 0 {width} {height}" width="100%" '
             f'height="{height}" role="img" aria-label="Family tree">']

    drawn_couples = set()
    for pid in fam.people:
        partner = fam.partner(pid)
        pair = tuple(sorted((pid, partner or "")))
        if not partner or partner not in pos or pair in drawn_couples:
            continue
        drawn_couples.add(pair)
        (x1, y1), (x2, y2) = pos[pid], pos[partner]
        parts.append(f'<line class="couple" x1="{min(x1, x2) + NODE_W // 2}" y1="{y1}" '
                     f'x2="{max(x1, x2) - NODE_W // 2}" y2="{y2}" />')

    for pid in fam.people:
        parents = [x for x in fam.parents(pid) if x in pos]
        if not parents or pid not in pos:
            continue
        # Descent hangs from the midpoint of the couple bar when both parents are drawn,
        # and straight down from the single known parent otherwise.
        px = sum(pos[p][0] for p in parents) // len(parents)
        py = max(pos[p][1] for p in parents) + NODE_H // 2
        cx, cy = pos[pid][0], pos[pid][1] - NODE_H // 2
        mid = (py + cy) // 2
        parts.append(f'<path id="e-{_attr(pid)}" class="edge" fill="none" '
                     f'd="M {px} {py} V {mid} H {cx} V {cy}" />')

    for pid, p in fam.people.items():
        if pid not in pos:
            continue
        x, y = pos[pid]
        shape = ("rect" if p.sex == "male" else "ellipse")
        if shape == "rect":
            node = (f'<rect x="{x - NODE_W // 2}" y="{y - NODE_H // 2}" width="{NODE_W}" '
                    f'height="{NODE_H}" rx="8" />')
        else:
            node = (f'<ellipse cx="{x}" cy="{y}" rx="{NODE_W // 2}" ry="{NODE_H // 2}" />')
        sub = "" if fam.is_sequenced(pid) else (
            f'<text class="sub" x="{x}" y="{y + 16}" text-anchor="middle">not sequenced</text>')
        # The count rides on a pip pinned to the node's top-right corner rather than
        # floating inside it: an ellipse has no corner to tuck text into, and a number
        # sitting on the outline reads as a badge instead of as part of the name.
        bx, by = x + NODE_W // 2 - 6, y - NODE_H // 2 + 4
        parts.append(
            f'<g id="n-{_attr(pid)}" class="node" data-person="{_attr(pid)}">'
            f'{node}<text class="name" x="{x}" y="{y + (0 if sub else 5)}" '
            f'text-anchor="middle">{_esc(p.display_name)}</text>{sub}'
            f'<text class="relab" data-rel-of="{_attr(pid)}" x="{x}" '
            f'y="{y - 12 if sub else y + 19}" text-anchor="middle"></text>'
            f'<circle class="pip" cx="{bx}" cy="{by}" r="10" />'
            f'<text class="badge" x="{bx}" y="{by + 4}" text-anchor="middle"></text></g>')

    parts.append("</svg>")
    return "".join(parts)


_TREE_STYLE = """
.legend { display:flex; flex-wrap:wrap; gap:.5rem 1.1rem; margin:.6rem 0 1rem;
  font-size:.86rem; color:var(--fg2); align-items:center; }
.legend span::before { content:""; display:inline-block; width:.8rem; height:.8rem;
  border-radius:4px; margin-right:.4rem; vertical-align:-1px; }
.legend .hom::before { background:var(--hom); } .legend .het::before { background:var(--het); }
.legend .non::before { border:1.5px solid var(--tick); }
.legend .unk::before { border:1.5px dashed var(--tick); }
.legend .dn::before { border:2px solid var(--dn); }
.legend span:last-child::before { display:none; }
.treewrap { border:1px solid var(--line); border-radius:20px; padding:1rem;
  background:var(--panel); overflow-x:auto; }
#tree text { font-family:var(--sans); font-size:13px; font-weight:600; fill:var(--fg); }
#tree text.sub { font-size:10px; font-weight:400; fill:var(--muted); }
#tree text.relab { font-size:10px; font-weight:400; fill:var(--muted); }
#tree .node rect, #tree .node ellipse { fill:var(--panel); stroke:var(--edge);
  stroke-width:1.5; }
#tree .node.some rect, #tree .node.some ellipse { stroke:var(--accent); stroke-width:2.5; }
#tree .node.hom rect, #tree .node.hom ellipse { fill:var(--hom); stroke:var(--hom); }
#tree .node.het rect, #tree .node.het ellipse { fill:var(--het); stroke:var(--het); }
#tree .node.hom text.name, #tree .node.hom text.relab { fill:var(--on-hom); }
#tree .node.het text.name, #tree .node.het text.relab { fill:var(--on-het); }
#tree .node.unk rect, #tree .node.unk ellipse { fill:var(--raised); stroke-dasharray:4 3; }
#tree .node.dn rect, #tree .node.dn ellipse { stroke:var(--dn); stroke-width:3.5; }
#tree .node.ego rect, #tree .node.ego ellipse { stroke:var(--fg); stroke-width:3;
  stroke-dasharray:none; }
#tree text.badge { font-size:11px; font-weight:700; fill:var(--bg); }
#tree .pip { fill:none; stroke:none; }
#tree .node.hom .pip, #tree .node.het .pip, #tree .node.some .pip { fill:var(--fg);
  stroke:var(--panel); stroke-width:2; }
#tree .edge { stroke:var(--edge); stroke-width:2; }
#tree .edge.on { stroke:var(--accent); stroke-width:3.5; }
#tree .edge.dn { stroke:var(--dn); stroke-width:3.5; stroke-dasharray:6 4; }
#tree .couple { stroke:var(--edge); stroke-width:2; }
table.vars tbody tr { cursor:pointer; }
table.vars td:first-child, table.vars td:nth-child(2) { white-space:nowrap; }
table.vars td.gene { font-family:var(--mono); font-weight:500; }
table.vars tbody tr.sel { background:var(--sel); box-shadow:inset 3px 0 0 var(--accent); }
#overview { min-height:36px; padding:0 .9rem; border:1px solid var(--track);
  border-radius:10px; background:var(--raised); color:var(--fg); cursor:pointer;
  font:inherit; font-size:.86rem; font-weight:600; }
#sel { font-family:var(--display); font-size:1.15rem; font-weight:600; color:var(--fg); }
.pill { display:inline-block; padding:.05rem .5rem; border-radius:10px; font-size:.76rem;
  border:1px solid var(--line); color:var(--muted); }
.pill.dn { border-color:var(--dn); color:var(--dn); }
.pill.common { border-color:var(--het-text); color:var(--het-text); }
.pill.low { border-color:var(--hom-text); color:var(--hom-text); }
"""

# Restyles the tree for the selected variant. The drawing never changes shape — only the
# classes on nodes and descent edges do — so this stays a few lines and the page works as a
# static tree (everyone's carrier badge) before any click.
_TREE_SCRIPT = """
<script>
const V = __VARIANTS__, PEOPLE = __PEOPLE__;
const rows = [...document.querySelectorAll('table.vars tbody tr')];
const cap = document.getElementById('sel');

// ClinVar concatenates every submitted condition with '|'; the first real one is enough
// for a one-line caption, and 'not provided' is not a condition.
function disease(d){
  const first = (d||'').split('|').filter(x=>x && x!=='not_provided')[0];
  return first ? first.replace(/_/g,' ') : '';
}
function clear(){
  document.querySelectorAll('#tree .node').forEach(n=>
    n.classList.remove('hom','het','some','unk','dn'));
  document.querySelectorAll('#tree .edge').forEach(e=>e.classList.remove('on','dn'));
}
function badge(p, text){
  const g = document.getElementById('n-'+p);
  if (g) g.querySelector('.badge').textContent = text;
  return g;
}
function overview(){
  clear(); rows.forEach(r=>r.classList.remove('sel'));
  cap.textContent = 'All actionable variants \u2014 each person is badged with how many they carry.';
  for (const p of PEOPLE){
    const n = V.filter(v=>v.carriers[p]).length;
    const g = badge(p, n ? n : '');
    if (g && n) g.classList.add('some');
  }
}
function show(i){
  clear(); rows.forEach((r,j)=>r.classList.toggle('sel', i===j));
  const v = V[i], d = disease(v.disease);
  cap.textContent = v.gene+' '+v.label+' \u2014 '+v.sig.replace(/_/g,' ')+(d ? ' \u00b7 '+d : '');
  for (const p of PEOPLE){
    const c = v.carriers[p];
    const g = badge(p, c ? (c.zyg==='HOM' ? '2' : '1') : '');
    if (!g) continue;
    // classList.add('') throws (empty token) and would abort this loop at the first
    // non-carrier, leaving the rest of the tree showing the previous selection.
    const cls = c ? (c.zyg==='HOM' ? 'hom' : 'het') : (v.unknown.includes(p) ? 'unk' : '');
    if (cls) g.classList.add(cls);
    if (!c) continue;
    const e = document.getElementById('e-'+p);
    if (!e) continue;
    if (v.origin[p]==='de_novo'){ g.classList.add('dn'); e.classList.add('dn'); }
    else if (['father','mother','both'].includes(v.origin[p])) e.classList.add('on');
  }
}
rows.forEach((r,i)=>r.addEventListener('click',()=>show(i)));
document.getElementById('overview').addEventListener('click',overview);
overview();
</script>
"""


def _table(variants: list[dict], fam) -> str:
    """One row per variant: gene, change, classification, who carries it, where it came
    from. Clicking a row drives the tree above; the table is complete without JavaScript."""
    from .render import _attr, _esc

    head = ("<table class=\"vars\"><thead><tr><th>Gene</th><th>Variant</th>"
            "<th>Classification</th><th>gnomAD AF</th><th>Carriers</th>"
            "<th>Inheritance</th></tr></thead><tbody>")
    body = []
    for v in variants:
        common = v["af"] is not None and v["af"] >= config.RARE_AF_THRESHOLD
        af = "—" if v["af"] is None else f"{v['af']:.3g}"
        who, inh = [], []
        for pid in fam.people:
            c = v["carriers"].get(pid)
            if not c:
                continue
            dn = fam.people[pid].display_name
            # "(1/2)" — called by one of this person's two callsets. Terser than naming the
            # callers and it is the number that matters: partial support is the page's own
            # confidence signal, and the T2T-lifted calls disagree often enough that
            # spelling it out per name buried the table.
            total = len(fam.callsets(pid))
            frac = (f' <span class="pill" title="called by {len(c["callsets"])} of '
                    f'{total} callsets">{len(c["callsets"])}/{total}</span>'
                    if total > 1 and len(c["callsets"]) < total else "")
            low = ""
            if c.get("low_ab"):
                a, n = c["low_ab"]
                low = (f' <span class="pill low" title="only {a} of {n} reads carry the '
                       f'variant; a real het sits near 50% — likely sequencing or mapping '
                       f'artifact, review">low AF {a}/{n}</span>')
            who.append(f"{_esc(dn)}{' (hom)' if c['zyg'] == 'HOM' else ''}{frac}{low}")
            o = origin(v, pid, fam)
            if o and o != FOUNDER:
                cls = ' class="pill dn"' if o == DE_NOVO else ' class="pill"'
                inh.append(f"{_esc(dn)} <span{cls}>{_esc(ORIGIN_LABEL[o])}</span>")
        flags = ' <span class="pill common">common</span>' if common else ""
        body.append(
            f'<tr><td class="gene">{_esc(v["gene"])}</td>'
            f'<td><code>{_esc(v["label"])}</code>{flags}</td>'
            f'<td>{_esc(clinsig.cell(v["sig"]).replace("_", " "))}</td><td>{_esc(af)}</td>'
            f'<td>{", ".join(who) or "—"}</td><td>{"<br>".join(inh) or "—"}</td></tr>')
    return ('<div class="tablewrap">' + head + "".join(body)
            + "</tbody></table></div>")


def write_markdown(snap: str, variants: list[dict], fam) -> Path:
    """The same findings as the interactive page, as a markdown report: one row per variant,
    who carries it and where each carrier got it. The tree is a picture; this is the text a
    reader — or a retriever — can quote."""
    from .util import disease_name, write_report

    def who(v: dict) -> str:
        out = []
        for p, c in v["carriers"].items():
            name = fam.people[p].display_name if p in fam.people else p
            zyg = "homozygous" if c["zyg"] == "HOM" else "heterozygous"
            low = ", low allele fraction" if c.get("low_ab") else ""
            out.append(f"{name} ({zyg}, {ORIGIN_LABEL[origin(v, p, fam)]}{low})")
        return "; ".join(out)

    unseq = [fam.people[p].display_name for p in fam.people if not fam.is_sequenced(p)]
    lines = [
        "# Variant inheritance — family",
        "",
        f"- **Variants:** {len(variants)} ClinVar pathogenic / likely-pathogenic variant(s) "
        "carried by someone in the family",
        f"- **People sequenced:** {sum(1 for p in fam.people if fam.is_sequenced(p))} of "
        f"{len(fam.people)}" + (f" (not sequenced: {', '.join(unseq)})" if unseq else ""),
        "",
        "## Who carries what, and from whom",
        "",
        "| Gene | Variant | ClinVar | Condition | Carriers |",
        "|------|---------|---------|-----------|----------|",
        *([f"| {v['gene']} | {v['label']} | {clinsig.cell(v['sig']).replace('_', ' ')} | "
           f"{disease_name(v['disease'])} | {who(v)} |" for v in variants]
          or ["| _none_ | | | | |"]),
        "",
        "## Method & caveats",
        "- Carrier status is per person: a variant called by any of a person's callsets "
        "counts, since alternate callsets are the same genome called twice.",
        "- Origin is inferred from which sequenced parent carries the variant, not from "
        "read-backed phasing. A parent must be sequenced for anything to be said about "
        "descent; otherwise the origin is unresolved.",
        "- A de-novo candidate is an argument from absence — the commonest cause is a variant "
        "the parent carries that their caller missed. The `denovo` report settles it.",
        "- Research-grade, not a clinical diagnosis.",
        "",
    ]
    return write_report(reportpaths.family_report("inheritance", snap), lines)


def build_page(snap: str | None = None, out_dir: Path | None = None) -> Path:
    """Write `reports/html/<snapshot>/inheritance.html` — the family-wide inheritance tree —
    and its markdown counterpart."""
    import json

    from . import pedigree
    from .render import (_SHELL_STYLE, _STYLE, _attr, _esc, _inline, perspective_control,
                         perspective_script, snapshot_switcher)

    fam = pedigree.load(sequenced=pedigree.sequenced_ids())
    snap = snap or latest_snapshot()
    if not snap:
        raise SystemExit("No snapshot yet. Run `snapshot` first.")
    if not loaded(snap):
        raise SystemExit(
            f"No variant rows for snapshot {snap} in {config.DUCKDB_FILE.name} — this "
            "page is built from the database, not from the markdown reports beside it. "
            "Build it on the machine that ran the scan, or load that snapshot here first.")
    out_dir = out_dir or (reportpaths.html_root() / snap)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "inheritance.html"

    variants = collect(snap, fam)
    write_markdown(snap, variants, fam)
    pos, width, height = layout(fam)
    people = [p for p in fam.people]

    # Everything the script needs, computed here in Python. `unknown` is the set of people
    # whose carrier status for this variant is not established either way — nobody with a
    # callset is ever "unknown", so an unsequenced relative reads as a gap in the data
    # rather than as a non-carrier.
    payload = [{
        "label": v["label"], "gene": v["gene"], "sig": v["sig"], "disease": v["disease"],
        "carriers": {p: {"zyg": c["zyg"]} for p, c in v["carriers"].items()},
        "origin": {p: origin(v, p, fam) for p in v["carriers"]},
        "unknown": [p for p in people if not fam.is_sequenced(p)],
    } for v in variants]

    nav = ('<nav class="tree" aria-label="Reports"><p class="brand">' + theme.LOGO +
           '<a href="index.html">Genome reports</a></p>'
           '<p class="who">Family-wide</p>'
           + perspective_control(fam) +
           '<details open><summary>Family-wide</summary><ul>'
           '<li class="here"><a href="inheritance.html">Variant inheritance tree</a></li>'
           "</ul></details>"
           '<details open><summary>People</summary><ul class="fam">'
           + "".join(
               f'<li><a href="index_{_attr(p)}.html">{_esc(fam.people[p].display_name)}</a>'
               f'<span class="rel" data-rel-of="{_attr(p)}"></span></li>'
               for p in people if fam.is_sequenced(p))
           + "</ul></details>" + snapshot_switcher(out_dir, "inheritance.html") + "</nav>")

    n_dn = sum(1 for v in variants for p in v["carriers"] if origin(v, p, fam) == DE_NOVO)
    n_disc = sum(1 for v in variants if discordant(v, fam))
    n_low = sum(1 for v in variants for c in v["carriers"].values() if c.get("low_ab"))
    body = [
        f'<p class="meta">Snapshot {_esc(snap)}</p>',
        "<h1>Who carries what, and from whom</h1>",
        (f'<div class="disclaimer">{_inline(config.DEMO_NOTICE)}</div>'
         if config.DEMO_NOTICE else ""),
        f'<p class="lede" style="max-width:none">'
        f"{len(variants)} ClinVar pathogenic / likely-pathogenic variant(s) carried by "
        f"someone in this family. Click a row to trace it through the tree."
        + (f" {n_disc} of them are called by only some of a carrier's callsets "
           f"(the <code>n/m</code> marker) — treat those as the least certain rows."
           if n_disc else "")
        + f" A <span class=\"pill low\">low AF</span> marker means fewer than "
        f"{config.LOW_VAF_THRESHOLD:.0%} of that person's reads carry the variant, where a "
        f"real heterozygote sits near 50% — the typical signature of a sequencing or "
        f"mapping artifact"
        + (f" ({n_low} carrier call(s) here)" if n_low else " (none here)")
        + ". Such a call is shown, not hidden; review it before acting on it.</p>",
        f'<p id="sel" class="note"></p>',
        '<div class="legend"><span class="hom">homozygous</span>'
        '<span class="het">heterozygous</span><span class="non">not carried</span>'
        '<span class="unk">not sequenced</span>'
        '<span class="dn">de-novo candidate</span>'
        '<span><button id="overview" type="button">Show all variants</button></span></div>',
        f'<div class="treewrap">{_svg(fam, pos, width, height)}</div>',
        "<h2>Variants</h2>",
        _table(variants, fam),
    ]
    if n_dn:
        body.append(
            f"<h2>About the {n_dn} de-novo candidate(s)</h2>"
            "<p><b>Candidate, not confirmed.</b> A variant is marked de-novo here when "
            "neither parent's callsets contain it. That is an argument from absence, and "
            "the commonest reason for it is a variant the parent genuinely carries that "
            "their caller missed — not a new mutation. The <code>denovo</code> stage settles "
            "it properly by force-calling both parents and piling up their raw reads; read "
            "that person's <b>Trio de-novo variants</b> report for the verdict.</p>")
    body.append(
        "<h2>How this is derived</h2><ul>"
        "<li>Carrier status is per <b>person</b>: a variant called by any of that person's "
        "callsets counts, since alternate callsets are the same genome called twice.</li>"
        "<li>A parent must be sequenced for the page to say anything about descent; "
        "otherwise the variant is listed as <i>unresolved</i>.</li>"
        "<li>The list is gated to ClinVar pathogenic / likely-pathogenic only. Common "
        f"variants (gnomAD AF ≥ {config.RARE_AF_THRESHOLD}) are shown but flagged.</li>"
        "<li>Which parent transmitted a variant is inferred from who carries it, not from "
        "read-backed phasing; the <code>phase</code> and <code>segregation</code> reports "
        "resolve cis/trans where it matters.</li></ul>")

    body.append('<p class="foot">Research-grade, not a clinical diagnosis. Confirm any '
                "actionable finding with a qualified clinician.</p>")
    script = (_TREE_SCRIPT.replace("__VARIANTS__", json.dumps(payload))
              .replace("__PEOPLE__", json.dumps(people))) + perspective_script(fam)
    html = (
        "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        "<title>Variant inheritance</title>\n"
        f"<style>{_STYLE}{_SHELL_STYLE}{_TREE_STYLE}</style>\n{theme.SCRIPT}\n</head>\n"
        f"<body class=\"shell\">\n{nav.replace('</nav>', theme.SWITCHER + '</nav>', 1)}\n"
        "<main class=\"content\">\n"
        + "\n".join(body) + f"\n{script}\n</main>\n</body>\n</html>\n")
    from .render import with_fonts
    out.write_text(with_fonts(html, out), encoding="utf-8")
    log.info("Inheritance tree: %s (%d variant(s), %d de-novo candidate(s))",
             out, len(variants), n_dn)
    return out
