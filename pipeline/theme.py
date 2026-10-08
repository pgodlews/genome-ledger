"""The look of every HTML page: colour tokens, light / dark / auto, typography, components.

One stylesheet, inlined into each page, because a report has to stay a single self-contained
file. Nothing is fetched: the font stacks name Bricolage Grotesque and IBM Plex first and
fall back to system fonts, so a machine without them still gets a clean page and no page
ever makes a network request about someone's genome.

Colour carries meaning and is used for nothing else:
  --hom   two copies of a variant (or a result in the top band)
  --het   one copy (or an elevated result)
  --accent  navigation, selection, and neutral data marks

Theme: `data-theme` on <html> is "light", "dark" or "auto". With no attribute at all (script
disabled) the page follows the system setting. The choice is kept in localStorage and shared
by every page, including a report shown in the dashboard's reading pane.
"""

from __future__ import annotations

THEME_KEY = "genome-theme"

_DARK = ("--bg:#0B0E14;--nav:#0F131C;--panel:#121723;--raised:#19202E;--line:#1E2534;"
         "--track:#283044;--edge:#3A4560;--tick:#4A556C;--fg:#EAEFF7;--fg2:#C9D2E3;"
         "--muted:#9BA7BD;--hom:#FF7566;--hom-text:#FF8A7D;--on-hom:#2A0B07;--het:#F5C451;"
         "--het-text:#F5C451;--on-het:#2A2005;--warn-bg:#1E1A10;--warn-line:#3A3320;"
         "--warn-fg:#F3E3B3;--sel:#162029;--zone:#3A2622;--accent:#6FE3C1;--dn:#B79CFF;"
         "color-scheme:dark;")
_LIGHT = ("--bg:#F3F5F9;--nav:#FFFFFF;--panel:#FFFFFF;--raised:#EDF0F5;--line:#DCE1EA;"
          "--track:#D5DBE6;--edge:#AEB8C9;--tick:#8793A8;--fg:#121723;--fg2:#37415A;"
          "--muted:#5B667D;--hom:#FF7566;--hom-text:#B3261E;--on-hom:#2A0B07;--het:#F5C451;"
          "--het-text:#8A5A00;--on-het:#2A2005;--warn-bg:#FFF6DD;--warn-line:#EBD9A0;"
          "--warn-fg:#5A4300;--sel:#E3F4EF;--zone:#F7D3CD;--accent:#0B7A63;--dn:#6A3FB5;"
          "color-scheme:light;")

_SANS = "'IBM Plex Sans',system-ui,-apple-system,'Segoe UI',Roboto,Helvetica,Arial,sans-serif"
_DISPLAY = "'Bricolage Grotesque','IBM Plex Sans',system-ui,-apple-system,'Segoe UI',sans-serif"
_MONO = "'IBM Plex Mono',ui-monospace,'SF Mono',Menlo,Consolas,monospace"

# Bundled fonts (pipeline/fonts/, SIL Open Font License — the licence files sit beside them).
# `render` copies them to reports/html/assets/fonts/ and each page points at that folder by
# a relative path, so nothing is fetched from the network and an air-gapped machine gets
# the intended typography. A report opened without the folder falls back to system fonts.
FONTS_TOKEN = "/*@fonts@*/"
_FACES = (("Bricolage Grotesque", "200 800", "BricolageGrotesque-Variable.woff2"),
          ("IBM Plex Sans", "400", "IBMPlexSans-Regular.woff2"),
          ("IBM Plex Sans", "500", "IBMPlexSans-Medium.woff2"),
          ("IBM Plex Sans", "600", "IBMPlexSans-SemiBold.woff2"),
          ("IBM Plex Mono", "400", "IBMPlexMono-Regular.woff2"),
          ("IBM Plex Mono", "500", "IBMPlexMono-Medium.woff2"))


def font_faces(prefix: str) -> str:
    """@font-face rules for the bundled fonts, `prefix` being the relative path from the
    page to the fonts folder. No double quotes and no quoted url(): the same text is placed
    inside a `srcdoc` attribute for the reading pane's first page."""
    return "".join(
        f"@font-face{{font-family:'{family}';font-weight:{weight};font-style:normal;"
        f"font-display:swap;src:url({prefix}/{name})}}" for family, weight, name in _FACES)


# --- every page ------------------------------------------------------------------
STYLE = FONTS_TOKEN + f"""
:root,:root[data-theme="dark"]{{{_DARK}}}
:root[data-theme="light"]{{{_LIGHT}}}
@media (prefers-color-scheme: light){{:root:not([data-theme]),:root[data-theme="auto"]{{{_LIGHT}}}}}
:root{{--sans:{_SANS};--display:{_DISPLAY};--mono:{_MONO};}}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--fg); font-family:var(--sans);
  font-size:15px; line-height:1.55; }}
a {{ color:var(--accent); }}
main.report {{ max-width:960px; margin:0 auto; padding:2rem 2rem 4rem; }}
h1 {{ font-family:var(--display); font-weight:700; font-size:2.9rem; line-height:1.04;
  letter-spacing:-.03em; margin:.2rem 0 .8rem; }}
h2 {{ font-family:var(--display); font-weight:600; font-size:1.4rem; letter-spacing:-.01em;
  margin:2.2rem 0 .8rem; }}
h3 {{ font-family:var(--display); font-weight:600; font-size:1.08rem; margin:1.6rem 0 .5rem;
  color:var(--fg2); }}
code {{ font-family:var(--mono); font-size:.88em; background:var(--raised);
  padding:.1em .4em; border-radius:6px; }}
hr {{ border:0; border-top:1px solid var(--line); margin:1.8rem 0; }}
ul {{ padding-left:1.25rem; }} li {{ margin:.25rem 0; }}
p.note {{ color:var(--muted); }}
h1 + p.note {{ font-family:var(--mono); font-size:.8rem; margin:-.2rem 0 1.4rem; }}
.disclaimer {{ background:var(--warn-bg); border:1px solid var(--warn-line);
  color:var(--warn-fg); padding:.9rem 1.1rem; border-radius:14px; margin:1rem 0 1.4rem;
  font-size:.93rem; }}
.disclaimer strong, .disclaimer b {{ color:var(--het-text); }}
.tablewrap {{ overflow-x:auto; margin:.7rem 0 1.5rem; border:1px solid var(--line);
  border-radius:16px; background:var(--panel); }}
table {{ border-collapse:collapse; width:100%; font-size:.92rem; }}
th,td {{ padding:.7rem .95rem; text-align:left; vertical-align:top;
  border-top:1px solid var(--line); }}
thead th {{ border-top:0; font-weight:600; font-size:.74rem; letter-spacing:.07em;
  text-transform:uppercase; color:var(--muted); white-space:nowrap; }}
tbody tr:hover {{ background:var(--raised); }}
.pct {{ display:inline-flex; align-items:center; gap:.7rem; white-space:nowrap; }}
.pct b {{ font-family:var(--mono); font-weight:500; min-width:3.4em; text-align:right; }}
.track {{ position:relative; display:inline-block; width:150px; height:8px;
  border-radius:4px; background:var(--track); }}
.track::before {{ content:""; position:absolute; left:90%; right:0; top:0; bottom:0;
  border-radius:0 4px 4px 0; background:var(--zone); }}
.track::after {{ content:""; position:absolute; left:50%; top:-3px; width:1px; height:14px;
  background:var(--tick); }}
.track i {{ position:absolute; top:-4px; width:16px; height:16px; margin-left:-8px;
  border-radius:50%; background:var(--accent); z-index:1; }}
.track i.hi {{ background:var(--hom); }} .track i.up {{ background:var(--het); }}
.seg {{ display:inline-flex; padding:3px; border-radius:12px; border:1px solid var(--track);
  background:var(--raised); }}
.seg button {{ flex:1 1 0; min-height:36px; padding:0 .9rem; border:0; border-radius:9px;
  background:transparent; color:var(--fg2); font:inherit; font-size:.82rem; font-weight:600;
  cursor:pointer; }}
.seg button[aria-pressed="true"] {{ background:var(--fg); color:var(--bg); }}
.topbar {{ display:flex; justify-content:flex-end; margin-bottom:.6rem; }}
.framed .topbar {{ display:none; }}
@media print {{
  :root,:root[data-theme]{{{_LIGHT}}}
  body {{ background:#fff; }} main.report {{ max-width:none; padding:0; }}
  .topbar,.seg {{ display:none; }} .tablewrap {{ border:0; }}
}}
"""

# --- navigation pages (tree on the left, content on the right) ---------------------
SHELL_STYLE = """
body.shell { display:flex; align-items:stretch; min-height:100vh; }
nav.tree { flex:0 0 17rem; padding:1.6rem 1.1rem 2rem; border-right:1px solid var(--line);
  background:var(--nav); overflow-y:auto; max-height:100vh; position:sticky; top:0;
  font-size:.92rem; display:flex; flex-direction:column; gap:.2rem; }
nav.tree a { color:var(--fg2); text-decoration:none; }
nav.tree a:hover { color:var(--fg); }
nav.tree .brand { display:flex; align-items:center; gap:.6rem; margin:0 .6rem 1.1rem;
  font-family:var(--display); font-weight:700; font-size:1.12rem; letter-spacing:-.01em;
  color:var(--fg); }
nav.tree .brand a { color:var(--fg); }
nav.tree .brand svg { color:var(--accent); flex:0 0 auto; }
nav.tree .who { margin:-.7rem .6rem 1rem; color:var(--muted); font-size:.85rem; }
nav.tree details { margin:.15rem 0 .7rem; }
nav.tree summary { cursor:pointer; padding:.3rem .75rem; list-style:none; user-select:none;
  font-size:.7rem; letter-spacing:.12em; text-transform:uppercase; color:var(--muted); }
nav.tree summary::-webkit-details-marker { display:none; }
nav.tree summary::after { content:"\\25BE"; float:right; transition:transform .12s ease; }
nav.tree details:not([open]) > summary::after { transform:rotate(-90deg); }
nav.tree ul { list-style:none; margin:.2rem 0 0; padding:0; }
nav.tree li { margin:0; border-radius:10px; }
nav.tree li > a, nav.tree li > span.unseq { display:inline-block; padding:.5rem .75rem; }
nav.tree li:hover { background:var(--raised); }
nav.tree li.here { background:var(--raised); color:var(--fg); font-weight:600;
  padding:.5rem .75rem; }
nav.tree li.here > a { padding:0; color:var(--fg); }
nav.tree ul.fam li { display:flex; align-items:center; justify-content:space-between;
  gap:.5rem; padding-right:.75rem; }
nav.tree .sub { padding-left:.6rem; }
nav.tree details.persp { margin:0 0 1rem; border:1px solid var(--track); border-radius:12px;
  background:var(--raised); }
nav.tree details.persp > summary { padding:.7rem .85rem; text-transform:none;
  letter-spacing:0; font-size:.85rem; }
nav.tree details.persp > summary b { color:var(--fg); }
nav.tree details.persp ul { padding:0 .3rem .4rem; }
nav.tree .seg { margin-top:auto; display:flex; }
.rel { color:var(--muted); font-size:.8rem; margin-left:.45rem; white-space:nowrap; }
.unseq { color:var(--muted); }
iframe.pane { flex:1 1 auto; min-width:0; border:0; height:100vh; position:sticky; top:0;
  background:var(--bg); }
main.content { flex:1 1 auto; min-width:0; max-width:1120px; margin:0 auto;
  padding:2.2rem 2.6rem 3.5rem; }
main.content h3 { font-family:var(--sans); font-size:.74rem; margin:1.2rem 0 .3rem;
  color:var(--muted); text-transform:uppercase; letter-spacing:.1em; }
main.content ul.fam { list-style:none; padding-left:0; }
.meta { font-family:var(--mono); font-size:.8rem; color:var(--muted); margin:0 0 .3rem; }
.hero { font-size:4.4rem; line-height:.98; letter-spacing:-.04em; margin:.2rem 0 .9rem; }
.lede { font-size:1.08rem; color:var(--fg2); max-width:36rem; margin:0 0 1.6rem; }
.card { display:block; padding:1.3rem 1.5rem; border-radius:20px; background:var(--panel);
  border:1px solid var(--line); color:var(--fg); text-decoration:none; }
a.card:hover { border-color:var(--edge); }
.card h2 { margin:0 0 .9rem; font-size:1.25rem; }
.card .k { font-size:.82rem; color:var(--muted); margin:0 0 .35rem; }
.card .big { font-family:var(--display); font-weight:700; font-size:3.2rem; line-height:1;
  margin:0 0 .35rem; }
.card .v { font-size:.85rem; color:var(--fg2); margin:0; }
.card .t { font-family:var(--display); font-weight:600; font-size:1.35rem;
  letter-spacing:-.01em; margin:0 0 .5rem; }
.c-hom { color:var(--hom-text); } .c-het { color:var(--het-text); }
.c-acc { color:var(--accent); } .c-mut { color:var(--muted); }
.grid { display:grid; gap:1rem; margin:0 0 1rem; }
.g4 { grid-template-columns:repeat(4,minmax(0,1fr)); }
.g3 { grid-template-columns:repeat(3,minmax(0,1fr)); }
.g2 { grid-template-columns:repeat(2,minmax(0,1fr)); }
.chips { display:flex; flex-wrap:wrap; gap:.4rem; font-family:var(--mono); font-size:.76rem; }
.chips span { padding:.25rem .55rem; border-radius:8px; background:var(--raised);
  color:var(--fg2); }
.badge { display:inline-block; padding:.3rem .6rem; border-radius:8px; font-size:.8rem;
  font-weight:600; }
.badge.hom { background:var(--hom); color:var(--on-hom); }
.badge.het { background:var(--het); color:var(--on-het); }
.legend2 { display:flex; gap:1.1rem; font-size:.82rem; color:var(--muted); }
.legend2 span::before { content:""; display:inline-block; width:10px; height:10px;
  border-radius:50%; margin-right:.4rem; background:var(--het); }
.legend2 span.hom::before { background:var(--hom); }
.gmap { display:flex; gap:4px; align-items:center; padding-top:2rem; }
.gmap .chr { height:14px; border-radius:7px; background:var(--track); position:relative;
  min-width:6px; }
.gmap .pin { position:absolute; top:-5px; width:24px; height:24px; margin-left:-12px;
  border-radius:50%; background:var(--het); border:4px solid var(--panel); }
.gmap .pin.hom { background:var(--hom); }
.gmap .lab { position:absolute; top:-2rem; transform:translateX(-50%); font-family:var(--mono);
  font-size:.76rem; font-weight:500; color:var(--het-text); white-space:nowrap; }
.gmap .lab.hom { color:var(--hom-text); }
.gaxis { display:flex; justify-content:space-between; font-family:var(--mono);
  font-size:.7rem; color:var(--muted); margin-top:.9rem; }
.finding { display:flex; gap:.9rem; align-items:flex-start; padding:.9rem 0;
  border-top:1px solid var(--line); }
.finding:first-of-type { border-top:0; padding-top:0; }
.finding .z { flex:0 0 44px; height:44px; border-radius:12px; display:flex;
  align-items:center; justify-content:center; font-family:var(--mono); font-weight:500;
  font-size:.82rem; background:var(--het); color:var(--on-het); }
.finding .z.hom { background:var(--hom); color:var(--on-hom); }
.finding .n { font-weight:600; font-size:1.02rem; margin:0; }
.finding .m { font-family:var(--mono); font-size:.76rem; color:var(--muted); margin:.15rem 0; }
.finding .s { font-size:.9rem; color:var(--fg2); margin:0; }
.prs { display:grid; grid-template-columns:minmax(0,11rem) minmax(0,1fr) 2.6rem;
  gap:.75rem; align-items:center; font-size:.9rem; }
.prs .track { width:100%; } .prs b { font-family:var(--mono); font-weight:500;
  text-align:right; }
.people a.card.me { border:2px solid var(--accent); background:var(--sel); }
.people .nm { display:flex; align-items:baseline; justify-content:space-between;
  font-family:var(--display); font-weight:600; font-size:1.6rem; letter-spacing:-.02em;
  margin:0 0 .9rem; }
.people .nm .rel { font-family:var(--sans); font-weight:400; }
.people .g { font-family:var(--mono); font-size:.76rem; color:var(--muted);
  margin:.9rem 0 0; }
.snaps { display:flex; align-items:center; flex-wrap:wrap; }
.snaps a, .snaps span.cur { padding:.8rem 1.1rem; border-radius:14px;
  border:1px solid var(--track); color:var(--fg); text-decoration:none;
  font-family:var(--mono); font-size:.88rem; }
.snaps span.cur { border:2px solid var(--accent); background:var(--sel); }
.snaps i { flex:0 0 3rem; height:2px; background:var(--track); }
.foot { font-size:.82rem; color:var(--muted); border-top:1px solid var(--line);
  padding-top:1.1rem; margin-top:1.6rem; }
@media (max-width: 900px) {
  .g4 { grid-template-columns:repeat(2,minmax(0,1fr)); }
  .g3,.g2 { grid-template-columns:minmax(0,1fr); }
  .hero { font-size:2.8rem; }
}
@media (max-width: 720px) {
  body.shell { display:block; }
  nav.tree { flex:none; width:auto; max-height:none; position:static;
    border-right:0; border-bottom:1px solid var(--line); }
  iframe.pane { display:block; width:100%; position:static; }
  main.content { padding:1.4rem 1.1rem 2.5rem; }
}
@media print { nav.tree { display:none; } }
"""

# Applied in <head>, before first paint, so a page never flashes the wrong theme. The same
# script wires the Light / Dark / Auto buttons and keeps a report in the reading pane in
# step with the page around it.
SCRIPT = """<script>
(function(){
var K=%r, R=document.documentElement;
function get(){try{return localStorage.getItem(K)||'auto'}catch(e){return 'auto'}}
function apply(t){
  R.setAttribute('data-theme',t);
  document.querySelectorAll('[data-theme-pick]').forEach(function(b){
    b.setAttribute('aria-pressed',String(b.getAttribute('data-theme-pick')===t));});
  document.querySelectorAll('iframe').forEach(function(f){
    try{f.contentDocument.documentElement.setAttribute('data-theme',t)}catch(e){}});
}
if(window.self!==window.top)R.classList.add('framed');
apply(get());
document.addEventListener('DOMContentLoaded',function(){
  apply(get());
  document.querySelectorAll('[data-theme-pick]').forEach(function(b){
    b.addEventListener('click',function(){
      var t=b.getAttribute('data-theme-pick');
      try{localStorage.setItem(K,t)}catch(e){}
      apply(t);});});
  document.querySelectorAll('iframe').forEach(function(f){
    f.addEventListener('load',function(){apply(get())});});
});
window.addEventListener('storage',function(e){if(e.key===K)apply(get())});
})();
</script>""" % THEME_KEY

SWITCHER = ('<div class="seg" role="group" aria-label="Theme">'
            '<button type="button" data-theme-pick="light">Light</button>'
            '<button type="button" data-theme-pick="dark">Dark</button>'
            '<button type="button" data-theme-pick="auto">Auto</button></div>')

LOGO = ('<svg width="26" height="26" viewBox="0 0 26 26" fill="none" stroke="currentColor" '
        'stroke-width="2" stroke-linecap="round" aria-hidden="true">'
        '<path d="M7 2c0 7 12 7 12 11S7 17 7 24"/><path d="M19 2c0 7-12 7-12 11s12 4 12 11"/>'
        '<path d="M9.5 7h7M9.5 19h7"/></svg>')
