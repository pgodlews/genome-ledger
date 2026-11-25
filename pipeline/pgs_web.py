"""The score picker at /pgs/ on `run.py serve`: browse the catalog snapshot, tick consent in
the family matrix, and score. A page over pgs_catalog and pgs_chosen; it decides nothing itself.

Reads are open like the rest of the server. **Writes** (recording consent, starting a scoring
job) are accepted only when all of these hold, because any web page the user visits can make
the browser POST to 127.0.0.1:
  - the request comes from a loopback address (the server may be bound to the LAN),
  - it carries this process's random token, which only this page knows,
  - its Origin, when the browser sends one, is this server, and its Host is a loopback name
    (a DNS-rebinding page would present its own host name).
"""

from __future__ import annotations

import json
import secrets
import threading
from html import escape
from urllib.parse import parse_qs

from . import config, pgs_catalog, pgs_chosen, reportpaths
from .util import log

TOKEN = secrets.token_urlsafe(24)
_MAX_BODY = 64 * 1024
_LOOPBACK_HOSTS = {"127.0.0.1", "localhost", "[::1]", "::1"}

_job_lock = threading.Lock()
_job: dict = {"running": False, "log": [], "result": None}


def _json(h, code: int, obj) -> None:
    body = json.dumps(obj).encode()
    h.send_response(code)
    h.send_header("Content-Type", "application/json")
    h.send_header("Cache-Control", "no-store")
    h.send_header("Content-Length", str(len(body)))
    h.end_headers()
    h.wfile.write(body)


def _q(parsed, key: str, default: str = "") -> str:
    return (parse_qs(parsed.query).get(key) or [default])[0]


def _results() -> list[dict]:
    out = []
    status = {(s["callset"], s["pgs_id"]): s for s in pgs_chosen.status()}
    for cs in pgs_chosen.scored_callsets():
        md = reportpaths.chosen_report("prs", cs)
        html = reportpaths.html_for(md)
        try:
            href = "/" + html.relative_to(config.GENOMES_ROOT.resolve()).as_posix()
        except ValueError:
            href = ""
        for r in pgs_chosen.results(cs):
            st = status.get((cs, r["pgs_id"]), {})
            out.append({"callset": cs, "pgs_id": r["pgs_id"], "traits": r["traits"],
                        "pe": r.get("pe"), "raw": r["raw"], "scored": r["scored"],
                        "current": st.get("current", True), "why": st.get("why", ""),
                        "report": href if html.exists() else ""})
    return out


def handle_get(h, parsed) -> None:
    rel = parsed.path.strip("/")
    try:
        if rel == "pgs":
            if not parsed.path.endswith("/"):
                h.send_response(302)
                h.send_header("Location", "/pgs/")
                h.send_header("Content-Length", "0")
                h.end_headers()
                return
            body = _PAGE.replace("__TOKEN__", escape(TOKEN)).encode()
            h.send_response(200)
            h.send_header("Content-Type", "text/html; charset=utf-8")
            h.send_header("Cache-Control", "no-store")
            h.send_header("Content-Security-Policy",
                          "default-src 'self'; style-src 'unsafe-inline'; "
                          "script-src 'unsafe-inline'; frame-ancestors 'none'")
            h.send_header("Content-Length", str(len(body)))
            h.end_headers()
            h.wfile.write(body)
            return
        api = rel[len("pgs/api/"):] if rel.startswith("pgs/api/") else None
        if api == "catalog":
            d = pgs_catalog.current_dir()
            if d is None:
                return _json(h, 404, {"error": "No catalog snapshot. Run `run.py pgs sync` "
                                               "(or `pgs import` a carried one)."})
            cat = pgs_catalog.load()
            return _json(h, 200, {"release": cat.release, "fetched": cat.fetched,
                                  "scores": len(cat.scores), "offline": config.PGS_OFFLINE})
        if api == "categories":
            return _json(h, 200, pgs_catalog.load().categories())
        if api == "traits":
            rows = pgs_catalog.load().traits(_q(parsed, "q"), _q(parsed, "category"))
            return _json(h, 200, rows[:300])
        if api == "scores":
            rows = pgs_catalog.load().scores_for(_q(parsed, "efo"))
            for r in rows:
                r["here"] = pgs_chosen.ready(r["pgs_id"])[0]
                r["performance_text"] = pgs_catalog.describe_performance(r["performance"])
            return _json(h, 200, rows)
        if api == "matrix":
            pids = [p.strip().upper() for p in _q(parsed, "pids").split(",") if p.strip()]
            m = pgs_chosen.matrix(pids)
            for s in m["scores"]:
                s["here"] = pgs_chosen.ready(s["pgs_id"])[0]
            return _json(h, 200, {**m, "warning": pgs_chosen.CHILD_WARNING,
                                  "relatives": pgs_chosen.RELATIVES_NOTE})
        if api == "results":
            return _json(h, 200, _results())
        if api == "job":
            with _job_lock:
                return _json(h, 200, dict(_job, log=list(_job["log"][-200:])))
        return _json(h, 404, {"error": "not found"})
    except SystemExit as e:
        return _json(h, 400, {"error": str(e)})
    except Exception as e:  # noqa: BLE001 — report it to the page, keep serving
        log.warning("pgs picker: %s", e)
        return _json(h, 500, {"error": str(e)})


def _write_allowed(h) -> str | None:
    """Why this write must be refused, or None."""
    if h.client_address[0] not in ("127.0.0.1", "::1"):
        return "writes are accepted only from this machine"
    if not secrets.compare_digest(h.headers.get("X-PGS-Token", ""), TOKEN):
        return "missing or wrong token (reload the page)"
    host = (h.headers.get("Host") or "").rsplit(":", 1)[0]
    if host not in _LOOPBACK_HOSTS:
        return "unexpected Host"
    origin = h.headers.get("Origin")
    if origin and origin.rstrip("/") != f"http://{h.headers.get('Host')}":
        return "cross-origin write refused"
    return None


def _run_job(pids: list[str], people: list[str] | None) -> None:
    def say(msg: str) -> None:
        with _job_lock:
            _job["log"].append(msg)
    try:
        out = pgs_chosen.score(pids, people, progress=say)
        summary = {"scored": len(out["scored"]), "skipped": len(out["skipped"]),
                   "failed": [f"{pid} → {cs}: {err}" for cs, pid, err in out["failed"]]}
        say(f"done: {summary['scored']} scored, {summary['skipped']} skipped (no consent), "
            f"{len(summary['failed'])} failed")
    except (Exception, SystemExit) as e:  # noqa: BLE001
        summary = {"error": str(e)}
        say(f"stopped: {e}")
    with _job_lock:
        _job["running"], _job["result"] = False, summary


def handle_post(h, parsed) -> None:
    why = _write_allowed(h)
    if why:
        return _json(h, 403, {"error": why})
    try:
        n = int(h.headers.get("Content-Length") or 0)
    except ValueError:
        n = 0
    if not 0 < n <= _MAX_BODY:
        return _json(h, 400, {"error": "bad request body"})
    try:
        data = json.loads(h.rfile.read(n))
    except ValueError:
        return _json(h, 400, {"error": "not JSON"})
    try:
        if parsed.path == "/pgs/api/consent":
            person = pgs_chosen.resolve_people([str(data.get("person", ""))])[0]
            pid = str(data.get("pgs_id", "")).upper()
            granted = bool(data.get("granted"))
            pgs_chosen.set_score_consent(person["id"], pid, granted)
            return _json(h, 200, {"person": person["id"], "pgs_id": pid,
                                  "allowed": pgs_chosen.may_score(
                                      person["id"], pgs_chosen.describe(pid)["efo_ids"])})
        if parsed.path == "/pgs/api/score":
            pids = [str(p).upper() for p in data.get("pgs_ids") or []]
            people = [str(p) for p in data.get("people") or []] or None
            if not pids:
                return _json(h, 400, {"error": "no scores selected"})
            with _job_lock:
                if _job["running"]:
                    return _json(h, 409, {"error": "a scoring job is already running"})
                _job.update(running=True, log=[], result=None)
            threading.Thread(target=_run_job, args=(pids, people), daemon=True).start()
            return _json(h, 202, {"started": True})
        return _json(h, 404, {"error": "not found"})
    except SystemExit as e:
        return _json(h, 400, {"error": str(e)})


_PAGE = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Choose polygenic scores</title>
<style>
:root{--bg:#fafaf8;--fg:#1d1d1b;--muted:#6b6b66;--line:#deded8;--card:#fff;--accent:#2f5d8a;
--warn:#8a5a00;--warnbg:#fff6e0;--ok:#2f7a3d}
@media (prefers-color-scheme:dark){:root{--bg:#161615;--fg:#e9e9e4;--muted:#9c9c95;
--line:#34342f;--card:#1f1f1d;--accent:#8db4dc;--warn:#f0c46a;--warnbg:#2c2414;--ok:#7cc48a}}
*{box-sizing:border-box}body{margin:0;font:14px/1.45 system-ui,sans-serif;background:var(--bg);
color:var(--fg)}header{padding:14px 20px;border-bottom:1px solid var(--line)}
h1{font-size:18px;margin:0}h2{font-size:15px;margin:0 0 8px}.muted{color:var(--muted)}
main{display:grid;grid-template-columns:230px 300px 1fr;gap:16px;padding:16px 20px}
@media (max-width:1000px){main{grid-template-columns:1fr}}
.panel{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:12px;
min-width:0}.list{max-height:60vh;overflow:auto}.item{padding:5px 6px;border-radius:5px;
cursor:pointer;display:flex;justify-content:space-between;gap:8px}.item:hover,.item.sel{
background:color-mix(in srgb,var(--accent) 14%,transparent)}input[type=search]{width:100%;
padding:7px;border:1px solid var(--line);border-radius:6px;background:var(--bg);color:var(--fg);
margin-bottom:8px}table{border-collapse:collapse;width:100%}th,td{text-align:left;padding:5px 6px;
border-bottom:1px solid var(--line);vertical-align:top}th{font-weight:600;font-size:12px;
color:var(--muted)}td.num{text-align:right;font-variant-numeric:tabular-nums}
button{font:inherit;padding:5px 10px;border-radius:6px;border:1px solid var(--line);
background:var(--card);color:var(--fg);cursor:pointer}button.primary{background:var(--accent);
color:var(--bg);border-color:var(--accent)}button:disabled{opacity:.5;cursor:default}
.warn{background:var(--warnbg);color:var(--warn);border-radius:6px;padding:8px 10px;margin:8px 0}
.wide{grid-column:1/-1}.tag{font-size:11px;padding:1px 6px;border-radius:9px;border:1px solid
var(--line);color:var(--muted);white-space:nowrap}.here{color:var(--ok);border-color:var(--ok)}
pre{white-space:pre-wrap;max-height:240px;overflow:auto;background:var(--bg);padding:8px;
border-radius:6px;border:1px solid var(--line)}
</style></head><body>
<header><h1>Choose polygenic scores</h1>
<div class="muted" id="cat">Loading the catalog snapshot…</div></header>
<main>
<section class="panel"><h2>Categories</h2><div class="list" id="cats"></div></section>
<section class="panel"><h2>Traits</h2><input type="search" id="q" placeholder="Search traits">
<div class="list" id="traits"></div></section>
<section class="panel"><h2 id="trait-title">Scores</h2><div class="muted" id="hint">Pick a trait.
Published performance is the median over European cohorts without age, sex or clinical
covariates; it says so when only covariate-adjusted figures exist.</div>
<div style="overflow:auto"><table id="scores"></table></div></section>
<section class="panel wide"><h2>Your selection and consent</h2>
<div id="basket" class="muted">Nothing selected yet.</div>
<div id="warn"></div><div style="overflow:auto"><table id="matrix"></table></div>
<p><button class="primary" id="go" disabled>Score where consented</button>
<button id="wish" disabled>Download wishlist (PGS ids)</button>
<span class="muted" id="cli"></span></p>
<pre id="job" hidden></pre></section>
<section class="panel wide"><h2>Computed</h2><div id="results" class="muted">None yet.</div></section>
</main>
<script>
const TOKEN="__TOKEN__";const $=s=>document.querySelector(s);
const el=(t,a={},...kids)=>{const e=document.createElement(t);for(const[k,v]of Object.entries(a)){
if(k==="class")e.className=v;else if(k.startsWith("on"))e.addEventListener(k.slice(2),v);
else e.setAttribute(k,v)}for(const c of kids)e.append(c instanceof Node?c:document.createTextNode(c??""));return e};
const api=async(p,opt)=>{const r=await fetch(p,opt);const j=await r.json();if(!r.ok)throw new Error(j.error||r.status);return j};
let category="",basket=JSON.parse(localStorage.getItem("pgs-basket")||"[]");
const save=()=>localStorage.setItem("pgs-basket",JSON.stringify(basket));
async function init(){try{const c=await api("/pgs/api/catalog");$("#cat").textContent=
`PGS Catalog release ${c.release} (copied ${c.fetched}), ${c.scores.toLocaleString()} scores`+
(c.offline?" — offline: scores not here must be imported as bundles":"")}catch(e){$("#cat").textContent=e.message;return}
const cats=await api("/pgs/api/categories");const box=$("#cats");box.replaceChildren(
el("div",{class:"item sel",onclick:ev=>pickCat("",ev.currentTarget)},el("span",{},"All traits")),
...cats.map(c=>el("div",{class:"item",onclick:ev=>pickCat(c.label,ev.currentTarget)},
el("span",{},c.label),el("span",{class:"muted"},String(c.traits)))));loadTraits();renderBasket();loadResults()}
function pickCat(c,node){category=c;document.querySelectorAll("#cats .item").forEach(n=>n.classList.remove("sel"));
node.classList.add("sel");loadTraits()}
let tq;$("#q").addEventListener("input",()=>{clearTimeout(tq);tq=setTimeout(loadTraits,200)});
async function loadTraits(){const t=await api(`/pgs/api/traits?q=${encodeURIComponent($("#q").value)}&category=${encodeURIComponent(category)}`);
$("#traits").replaceChildren(...t.map(x=>el("div",{class:"item",onclick:ev=>{document.querySelectorAll("#traits .item").forEach(n=>n.classList.remove("sel"));
ev.currentTarget.classList.add("sel");loadScores(x)}},el("span",{},x.label),el("span",{class:"muted"},String(x.scores)))))}
async function loadScores(t){$("#trait-title").textContent=`${t.label} (${t.efo})`;$("#hint").hidden=true;
const rows=await api(`/pgs/api/scores?efo=${encodeURIComponent(t.efo)}`);
$("#scores").replaceChildren(el("tr",{},...["","PGS ID","Variants","Published performance","Training ancestry","Released","Licence"].map(h=>el("th",{},h))),
...rows.map(s=>el("tr",{},el("td",{},s.unsupported?"":el("button",{onclick:()=>add(s.pgs_id),title:"add to selection"},basket.includes(s.pgs_id)?"✓":"Add")),
el("td",{},s.pgs_id," ",s.here?el("span",{class:"tag here"},"here"):"",s.unsupported?el("span",{class:"tag",title:s.unsupported},"cannot score"):""),
el("td",{class:"num"},s.n_variants.toLocaleString()+(s.engine==="plink2"?" (genome-wide)":"")),
el("td",{},s.performance_text),el("td",{},s.training_ancestry||s.gwas_ancestry||"—"),
el("td",{},(s.released||"").slice(0,10)),el("td",{},s.licence)))))}
function add(pid){if(!basket.includes(pid)){basket.push(pid);save();renderBasket()}}
function drop(pid){basket=basket.filter(p=>p!==pid);save();renderBasket()}
async function renderBasket(){$("#go").disabled=$("#wish").disabled=!basket.length;
if(!basket.length){$("#basket").textContent="Nothing selected yet.";$("#matrix").replaceChildren();$("#warn").replaceChildren();$("#cli").textContent="";return}
$("#cli").textContent="CLI: run.py pgs score "+basket.join(" ");
let m;try{m=await api(`/pgs/api/matrix?pids=${basket.join(",")}`)}catch(e){$("#basket").textContent=e.message;return}
$("#basket").replaceChildren(...m.scores.map(s=>el("span",{class:"tag",style:"margin-right:6px"},
`${s.pgs_id} · ${s.label}${s.sensitive?" ⚠":""}${s.here?" · here":""} `,el("a",{href:"#",onclick:ev=>{ev.preventDefault();drop(s.pgs_id)}},"✕"))));
$("#warn").replaceChildren(el("div",{class:"warn"},m.warning," ",m.relatives),
...(m.scores.some(s=>s.sensitive)?[el("div",{class:"warn"},"⚠ marks a trait the curated panel treats as severe with no established intervention.")]:[]));
$("#matrix").replaceChildren(el("tr",{},el("th",{},"Person"),...m.scores.map(s=>el("th",{title:s.label},s.pgs_id))),
...m.people.map(p=>el("tr",{},el("td",{},p.name),...m.scores.map(s=>{const cb=el("input",{type:"checkbox"});cb.checked=m.cells[p.id][s.pgs_id];
cb.addEventListener("change",async()=>{try{await api("/pgs/api/consent",{method:"POST",headers:{"Content-Type":"application/json","X-PGS-Token":TOKEN},
body:JSON.stringify({person:p.id,pgs_id:s.pgs_id,granted:cb.checked})});loadResults()}catch(e){alertBox(e.message);cb.checked=!cb.checked}});
return el("td",{},cb)}))))}
function alertBox(msg){$("#job").hidden=false;$("#job").textContent=msg}
$("#wish").addEventListener("click",()=>{const a=el("a",{href:URL.createObjectURL(new Blob([basket.join("\n")+"\n"],{type:"text/plain"})),download:"pgs-wishlist.txt"});a.click()});
$("#go").addEventListener("click",async()=>{try{await api("/pgs/api/score",{method:"POST",headers:{"Content-Type":"application/json","X-PGS-Token":TOKEN},
body:JSON.stringify({pgs_ids:basket})});$("#go").disabled=true;poll()}catch(e){alertBox(e.message)}});
async function poll(){const j=await api("/pgs/api/job");$("#job").hidden=false;$("#job").textContent=j.log.join("\n");
if(j.running)setTimeout(poll,1500);else{$("#go").disabled=!basket.length;loadResults();renderBasket()}}
async function loadResults(){const r=await api("/pgs/api/results");if(!r.length){$("#results").textContent="None yet.";return}
$("#results").replaceChildren(el("table",{},el("tr",{},...["Callset","PGS ID","Trait","EUR percentile","Scored","State",""].map(h=>el("th",{},h))),
...r.map(x=>el("tr",{},el("td",{},x.callset),el("td",{},x.pgs_id),el("td",{},x.traits.join("; ")),
el("td",{class:"num"},x.pe==null?"—":x.pe.toFixed(1)+"%"),el("td",{},x.scored),el("td",{},x.current?"current":x.why),
el("td",{},x.report?el("a",{href:x.report},"report"):"")))))}
init();
</script></body></html>
"""
