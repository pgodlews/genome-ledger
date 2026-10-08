"""Shared helpers: logging, shell execution, checksums, TSV manifests."""

from __future__ import annotations

import csv
import datetime as _dt
import hashlib
import logging
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Sequence

from . import config

log = logging.getLogger("genome")


def setup_logging(verbose: bool = False) -> None:
    config.LOGS_DIR.mkdir(parents=True, exist_ok=True)
    logfile = config.LOGS_DIR / f"{_dt.date.today():%Y-%m}.log"
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    handlers.append(logging.FileHandler(logfile))
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(message)s",
        datefmt="%H:%M:%S",
        handlers=handlers,
    )


def today() -> str:
    return f"{_dt.date.today():%Y-%m-%d}"


def run(cmd: Sequence[str], **kwargs) -> subprocess.CompletedProcess:
    """Run a command, logging it; raise on non-zero exit."""
    log.debug("$ %s", " ".join(str(c) for c in cmd))
    return subprocess.run([str(c) for c in cmd], check=True, **kwargs)


def require_tools(*tools: str) -> None:
    """Fail fast with a clear message if an external tool is missing."""
    import shutil

    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        raise SystemExit(
            f"Missing required tool(s): {', '.join(missing)}. "
            f"Run `python run.py setup` first."
        )


# Records the caller hard-filtered out are still present in the called VCF (call.py keeps them
# with FILTER carrying the failed filter's name; only `normalize` drops them). A typical GRCh38
# callset carries ~275k of them, ~5.6% of all records — MQ40/QD2/SOR3 failures that are exactly
# the calls we decided not to trust. Reading them as real genotypes contradicts the
# pipeline-wide "absent from a PASS-filtered VCF => hom-ref" convention. Every engine that
# reads genotypes out of a called VCF must apply this. ("PASS,." keeps FILTER='.' records,
# matching normalize.py: a vendor callset can be unfiltered throughout, and a bare PASS would
# score nothing at all for those samples.)
PASS_ONLY = ["-i", 'FILTER="PASS" || FILTER="."']


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_dirs() -> None:
    for d in config.ALL_DIRS:
        d.mkdir(parents=True, exist_ok=True)


# The closing section every report ends with, under one of these headings. Which one says
# how much it holds: plain caveats, caveats plus how the result was derived, or the
# auto-generated block of pinned versions and assumptions (limitations.render).
CLOSING_HEADINGS = ("## Caveats", "## Method & caveats", "## Assumptions & Known Limitations")


def _aligned(lines: Sequence[str]) -> list[str]:
    """The structure every report shares: a title, the findings, one closing section.

    `## Notes` is the same thing as `## Caveats` under another name, and a report whose last
    section is not a closing one gets the minimal closing section appended — so a reader, or
    a retriever splitting on headings, always finds the limits of a result in the same place."""
    out = ["## Caveats" if ln.strip() == "## Notes" else ln for ln in lines]
    h2 = [ln.strip() for ln in out if ln.startswith("## ")]
    if not h2 or h2[-1] not in CLOSING_HEADINGS:
        while out and not out[-1].strip():
            out.pop()
        out += ["", "## Caveats", "- Research-grade, not a clinical diagnosis.", ""]
    return out


def front_matter(path: Path, title: str = "") -> list[str]:
    """YAML front matter for a report in the reports tree: what it is, whose, as of when.
    Derived from where the file is filed (reportpaths.meta), so no writer can get it wrong.
    Empty for a file outside the tree."""
    import json
    from . import config, reportpaths
    m = reportpaths.meta(path)
    if m is None:
        return []
    fields = [("report", m["report"]), ("scope", m["scope"])]
    fields += [(k, m[k]) for k in ("person", "callset", "snapshot") if k in m]
    if "snapshot" in m:
        try:
            from .snapshot import load_manifest
            fields.append(("clinvar_date", load_manifest(m["snapshot"]).get("clinvar_date")))
        except Exception:  # noqa: BLE001 — a snapshot without a readable manifest
            pass
    fields.append(("generated", today()))
    if config.DEMO_NOTICE:
        fields.append(("demo", "true"))
    out = ["---"]
    if title:
        out.append(f"title: {json.dumps(title, ensure_ascii=False)}")
    out += [f"{k}: {v}" for k, v in fields if v not in (None, "")]
    return out + ["---", ""]


def write_report(path: Path, lines: Sequence[str]) -> Path:
    """Write a markdown report (the source of truth), then best-effort render its HTML.

    A report filed in the reports tree (reportpaths) gets the shared shape: front matter
    describing it, its title, its body, and a closing caveats section. HTML failure is
    logged and swallowed — it must never break report generation. Returns the `.md` path."""
    from . import config, reportpaths
    path = Path(path)
    lines = list(lines)
    in_tree = reportpaths.meta(path) is not None
    if in_tree:
        lines = _aligned(lines)
    if config.DEMO_NOTICE:
        at = 1 if lines and lines[0].startswith("# ") else 0
        lines[at:at] = ["", f"> {config.DEMO_NOTICE}"]
    if in_tree:
        title = lines[0][2:].strip() if lines and lines[0].startswith("# ") else ""
        lines = front_matter(path, title) + lines
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))
    try:
        from .render import render_formats  # lazy: keeps util import-light + cycle-free
        render_formats(path)
    except Exception as e:  # noqa: BLE001 — best-effort, markdown already written
        log.warning("Rendering failed for %s: %s (markdown is fine)", path.name, e)
    return path


# --- TSV manifest helpers ----------------------------------------------------
def append_tsv(path: Path, row: dict, fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # Not just exists(): a pre-created/crashed zero-byte file still needs the header,
    # otherwise read_tsv eats the first data row as the header.
    has_content = path.exists() and path.stat().st_size > 0
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        if not has_content:
            w.writeheader()
        w.writerow(row)


def read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def find_row(path: Path, **match) -> dict | None:
    for row in read_tsv(path):
        if all(row.get(k) == v for k, v in match.items()):
            return row
    return None


# --- reference FASTA lookups (shared by prs.py + forcecall.py) --------------
_FAIDX_BATCH = 20_000   # regions per `samtools faidx` invocation — stays well under the
# OS exec() argument-list limit. A ~740k-site marker panel (mostly absent in any one
# sample) once blew past it in one shot: `OSError: [Errno 7] Argument list too long`.


def _fai_contigs() -> set[str]:
    """Contig names present in the reference .fai — cheap, and read once per call."""
    fai = Path(f"{config.REF_FASTA}.fai")
    if not fai.exists():
        return set()
    return {ln.split("\t", 1)[0] for ln in fai.read_text().splitlines() if ln.strip()}


def _faidx_once(regions: list[tuple[str, str]]) -> dict[tuple[str, str], str]:
    """One `samtools faidx` invocation. Returns whatever it managed to emit."""
    args = [f"{c}:{p}-{p}" for c, p in regions]
    proc = subprocess.run(["samtools", "faidx", str(config.REF_FASTA), *args],
                          capture_output=True, text=True)
    out = proc.stdout
    if proc.returncode != 0 and not out.strip():
        err = " | ".join(proc.stderr.strip().splitlines()[-3:])
        raise SystemExit(
            f"samtools faidx failed against {config.REF_FASTA} (exit {proc.returncode}) "
            f"and returned nothing: {err}")
    got: dict[tuple[str, str], str] = {}
    cur = None
    for line in out.splitlines():
        if line.startswith(">"):
            c, rng = line[1:].split(":")
            cur = (c, rng.split("-")[0])
        elif cur:
            got[cur] = line.strip().upper()
            cur = None
    return got


def faidx_bases(positions: Iterable[tuple[str, str]]) -> dict[tuple[str, str], str]:
    """Reference base at each (chrom, pos) via batched `samtools faidx` against REF_FASTA.

    Batched so callers with large position lists (a big PGS, a curated force-call union)
    never hand a single exec() an argument list over the OS limit.

    `samtools faidx` ABORTS at the first region it cannot fetch: everything requested after
    it in the same invocation is never emitted, while the exit code (1) and the partial
    stdout look exactly like "those positions have no reference base". With a 20,000-region
    batch, one position on a contig missing from the FASTA silently cost up to 20,000 real
    lookups. Measured on the force-call site union: 7 positions on alt contigs took 44,217
    canonical ones down with them, so the union built at 29,599 of 73,886 sites and every
    consumer read the difference as "absent ⇒ hom-ref".

    Two defences. Positions on contigs the .fai does not carry are dropped up front, which
    removes the usual cause. And because samtools stops *at* the bad region, any batch that
    comes back short is retried from just past the offender, so a failure costs only itself
    rather than the remainder of the batch.
    """
    positions = list(positions)
    if not positions:
        return {}
    known = _fai_contigs()
    if known:
        unknown = sorted({c for c, _p in positions if c not in known})
        if unknown:
            n_bad = sum(1 for c, _p in positions if c not in known)
            log.warning("faidx: %d position(s) on %d contig(s) absent from %s.fai — skipped "
                        "up front so they cannot truncate a batch: %s", n_bad, len(unknown),
                        config.REF_FASTA.name, ", ".join(unknown[:5]))
        positions = [q for q in positions if q[0] in known]

    bases: dict[tuple[str, str], str] = {}
    unfetchable: list[tuple[str, str]] = []
    for i in range(0, len(positions), _FAIDX_BATCH):
        chunk = positions[i:i + _FAIDX_BATCH]
        while chunk:
            got = _faidx_once(chunk)
            bases.update(got)
            # samtools emits in request order and stops at the failure, so the resolved
            # prefix tells us exactly which region aborted the pass.
            n = 0
            while n < len(chunk) and chunk[n] in got:
                n += 1
            if n >= len(chunk):
                break
            unfetchable.append(chunk[n])
            chunk = chunk[n + 1:]        # strictly shorter each pass, so this terminates
    if unfetchable:
        log.warning("faidx: %d region(s) could not be fetched from %s and are reported as "
                    "absent, not as reference (e.g. %s); the rest of each batch was retried.",
                    len(unfetchable), config.REF_FASTA.name,
                    ", ".join(f"{c}:{p}" for c, p in unfetchable[:3]))
    return bases


# --- derived-artifact cache stamps ------------------------------------------
# Several stages cache an expensive derived artifact (a plink2 pgen, a thinned marker set)
# under a path keyed only by sample name, and skip the rebuild whenever that path exists.
# That silently reuses a stale artifact when its *input* changed — a re-called or re-ingested
# VCF, or a change to the build recipe itself. A stamp file next to the artifact records both,
# so the skip is conditional on the artifact still matching what produced it.


def file_token(path: Path) -> str:
    """Content identity for a pipeline artifact, used as a cache key.

    A hash, not (size, mtime). The old proxy assumed a content change at identical size and
    whole-second mtime was unrealizable — true of an edit in place, false of `rsync -a` and
    `cp -p`, which preserve mtime by design, so a file restored from backup or copied between
    hosts presents the previous token with different bytes. Every caller uses this to decide
    whether an expensive derived artifact may be reused, which is exactly where being wrong
    is expensive. Measured at ~2.6 GB/s here.
    """
    return f"sha256:{sha256(path)}"


def stamp_path(prefix: Path) -> Path:
    return prefix.with_suffix(".cache-stamp")


def stamp_is_current(prefix: Path, token: str) -> bool:
    sp = stamp_path(prefix)
    try:
        return sp.read_text().strip() == token
    except OSError:
        return False      # never stamped (pre-stamp cache) → treat as stale, rebuild once


def write_stamp(prefix: Path, token: str) -> None:
    stamp_path(prefix).write_text(token + "\n")


# --- ClinVar presentation helpers (shared by report + diff) ------------------
def review_stars(revstat: str | None) -> int:
    """ClinVar review status → 0–4 star confidence rating."""
    return config.REVIEW_STATUS_RANK.get((revstat or "").replace(" ", "_"), 0)


def disease_name(clndn: str | None) -> str:
    """First meaningful ClinVar disease name (CLNDN is '|'-separated)."""
    if not clndn:
        return "—"
    skip = {"not_provided", "not_specified", "see_cases", ""}
    parts = [p for p in clndn.split("|") if p.strip().lower() not in skip]
    name = (parts[0] if parts else clndn.split("|")[0]).replace("_", " ")
    return name[:55] + ("…" if len(name) > 55 else "")
