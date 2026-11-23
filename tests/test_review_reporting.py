#!/usr/bin/env python3
"""Regression tests for the DELTA, serving, ancestry, callability, SV, dashboard and shard fixes.

    uv run python tests/test_review_reporting.py

Everything runs in a throwaway GENOMES_ROOT on invented people, genes and coordinates.
"""

from __future__ import annotations

import gzip
import http.client
import json
import os
import pathlib
import re
import subprocess
import sys
import tempfile
import threading
import time
from contextlib import contextmanager

# Before any pipeline import: these tests write a database and reports.
os.environ["GENOMES_ROOT"] = tempfile.mkdtemp(prefix="genome_review_test_")

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from test_units import _failures, check  # noqa: E402

from pipeline import (acmg, ancestry, call, callability, config, dashboard, diff, load,  # noqa: E402
                      render, reportpaths, serve, snapshot, sv, theme)

PREV, CURR = "2031-01-01", "2031-02-01"
GENE = "GENE1"
VUS = {"sig": "Uncertain_significance", "rev": "criteria_provided,_single_submitter"}
FAMILY = ("id\tdisplay_name\tsex\tfather\tmother\tpartner\n"
          "Jan\tJan\tmale\t\t\t\nMaria\tMaria\tfemale\t\t\t\n")


@contextmanager
def patched(obj, **attrs):
    saved = {k: getattr(obj, k) for k in attrs}
    for k, v in attrs.items():
        setattr(obj, k, v)
    try:
        yield
    finally:
        for k, v in saved.items():
            setattr(obj, k, v)


def _row(con, snap, sid, chrom="1", pos=1000, ref="C", alt="T", sig="Pathogenic",
         rev="reviewed_by_expert_panel", csq="stop_gained", zyg="HET", dis="Invented_syndrome"):
    con.execute(
        f"INSERT INTO {load.TABLE} (snapshot_id, sample_id, chrom, pos, ref, alt, gene, "
        "consequence, clinvar_sig, clinvar_revstat, clinvar_alleleid, clinvar_disease, "
        "gnomad_af, sift, polyphen, zygosity) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        [snap, sid, chrom, pos, ref, alt, GENE, csq, sig, rev, "1", dis, None, None, None, zyg])


def _delta(fill) -> tuple[list[str], str]:
    """Run the real DELTA stage on a fresh database: (report rows for Jan, alert)."""
    config.DUCKDB_FILE.unlink(missing_ok=True)
    for s in (PREV, CURR):
        d = snapshot.snapshot_dir(s)
        d.mkdir(parents=True, exist_ok=True)
        (d / "manifest.json").write_text(json.dumps(
            {"snapshot_id": s, "clinvar_date": s, "vep_cache_version": "1"}))
    config.FAMILY_FILE.parent.mkdir(parents=True, exist_ok=True)
    config.FAMILY_FILE.write_text(FAMILY)
    con = load.connect()
    fill(con)
    con.close()
    with patched(acmg, panel_context=lambda: ({GENE}, {GENE: {"AD"}})):
        path, _n, alert = diff.diff(PREV, CURR)
    rows = [ln for ln in pathlib.Path(path).read_text().splitlines()
            if ln.startswith("| ") and "| Jan |" in ln]
    return rows, alert


def _callsets_cell(row: str) -> str:
    return row.rstrip(" |").rsplit("|", 1)[1].strip()


def test_delta_counts_only_compared_callsets() -> None:
    print("DELTA: a callset added in curr is not a disagreement:")

    def added_vendor(con):
        _row(con, PREV, "Jan", **VUS)
        _row(con, CURR, "Jan")
        _row(con, CURR, "Jan_vendor")
    rows, alert = _delta(added_vendor)
    check("a new vendor callset does not turn the finding into 1/2",
          len(rows) == 1 and _callsets_cell(rows[0]) == "-")
    check("the finding still alerts", "ACTIONABLE" in alert)

    def added_bare(con):
        _row(con, PREV, "Jan_t2t", **VUS)
        _row(con, CURR, "Jan_t2t")
        _row(con, CURR, "Jan")
    rows, alert = _delta(added_bare)
    check("a new bare callset carrying it means it is not T2T-only",
          len(rows) == 1 and "T2T only" not in rows[0])
    check("so the alert is kept", "ACTIONABLE" in alert)

    def added_without_it(con):
        _row(con, PREV, "Jan_t2t", **VUS)
        _row(con, CURR, "Jan_t2t")
        _row(con, CURR, "Jan", pos=2000, sig="Benign", csq="synonymous_variant")
    rows, alert = _delta(added_without_it)
    check("a new bare callset without it leaves the T2T cap on",
          len(rows) == 1 and "(T2T only)" in rows[0] and "ACTIONABLE" not in alert)


def test_delta_t2t_cap_holds_for_a_single_callset() -> None:
    print("DELTA: the T2T cap does not depend on how many callsets a person has:")

    def only_t2t(con):
        _row(con, PREV, "Jan_t2t", **VUS)
        _row(con, CURR, "Jan_t2t")
    rows, alert = _delta(only_t2t)
    check("a person's only callset being T2T cannot alert", "ACTIONABLE" not in alert)
    check("and the row says why",
          len(rows) == 1 and _callsets_cell(rows[0]) == "(T2T only)")


def test_delta_mt_is_not_discordant() -> None:
    print("DELTA: lifted chrM callsets do not count against an MT finding:")

    def all_three(con):
        for s in ("Jan", "Jan_vendor", "Jan_t2t"):
            _row(con, PREV, s, chrom="MT", pos=1234, ref="A", alt="G", **VUS)
            _row(con, CURR, s, chrom="MT", pos=1234, ref="A", alt="G")
    rows, _alert = _delta(all_three)
    check("an MT finding is not reported as 1/3",
          len(rows) == 1 and _callsets_cell(rows[0]) == "-")


def test_delta_vanished_variant() -> None:
    print("DELTA: a variant gone from curr keeps its description, not a score:")

    def vanished(con):
        _row(con, PREV, "Jan", zyg="HOM")
        _row(con, PREV, "Jan", pos=2000, sig="Benign", csq="synonymous_variant")
        _row(con, CURR, "Jan", pos=2000, sig="Benign", csq="synonymous_variant")
    rows, alert = _delta(vanished)
    check("zygosity survives", len(rows) == 1 and "| HOM |" in rows[0])
    check("disease survives", len(rows) == 1 and "Invented syndrome" in rows[0])
    check("a vanished stop-gain is not actionable",
          "ACTIONABLE" not in alert and len(rows) == 1 and "🔴" not in rows[0])


def _get(port: int, path: str) -> tuple[int, str]:
    c = http.client.HTTPConnection("127.0.0.1", port, timeout=5)
    c.putrequest("GET", path, skip_host=False)
    c.endheaders()
    r = c.getresponse()
    r.read()
    loc = r.getheader("Location") or ""
    c.close()
    return r.status, loc


def test_serve_percent_encoded_paths() -> None:
    print("serve: percent-encoded page names:")
    html = config.GENOMES_ROOT / "reports" / "html"
    (html / "2031-02-01").mkdir(parents=True, exist_ok=True)
    for name in ("index_Łena.html", "index_Mary Ann.html", "2031-02-01/index_Łena.html"):
        (html / name).write_text("ok")
    (config.GENOMES_ROOT / "private.html").write_text("no")
    httpd = serve.ThreadingHTTPServer(("127.0.0.1", 0), serve._Handler)
    port = httpd.server_address[1]
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    try:
        check("non-ASCII page is served",
              _get(port, "/reports/html/index_%C5%81ena.html")[0] == 200)
        check("page with a space is served",
              _get(port, "/reports/html/index_Mary%20Ann.html")[0] == 200)
        check("NUL byte is a clean 404",
              _get(port, "/reports/html/index%00.html")[0] == 404)
        check("encoded traversal stays blocked",
              _get(port, "/reports/html/..%2F..%2Fprivate.html")[0] == 403)
        check("encoded absolute path stays blocked",
              _get(port, "/%2Fetc%2Fhosts")[0] in (403, 404))
        status, loc = _get(port, "/reports/2031-02-01/index_%C5%81ena.html")
        check("old bookmark redirects with an encoded Location",
              status == 302 and loc == "/reports/html/2031-02-01/index_%C5%81ena.html")
    finally:
        httpd.shutdown()
        httpd.server_close()


def test_ancestry_subset_is_pass_only() -> None:
    print("ancestry: hard-filtered MT/Y calls do not reach the haplogroup tools:")
    d = pathlib.Path(tempfile.mkdtemp(dir=config.GENOMES_ROOT))
    vcf = d / "adam.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"p\">\n"
        "##FILTER=<ID=QD2,Description=\"q\">\n##contig=<ID=MT,length=16000>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAdam\n"
        "MT\t100\t.\tA\tG\t50\tPASS\t.\tGT\t1\n"
        "MT\t200\t.\tC\tT\t50\tQD2\t.\tGT\t1\n"
        "MT\t300\t.\tG\tA\t50\t.\t.\tGT\t1\n")
    gz = d / "adam.vcf.gz"
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(gz), str(vcf)], check=True)
    subprocess.run(["tabix", "-p", "vcf", str(gz)], check=True)
    out = ancestry._subset(gz, "MT", d / "adam.MT.vcf.gz")
    kept = subprocess.run(["bcftools", "query", "-f", "%POS\n", str(out)], check=True,
                          capture_output=True, text=True).stdout.split()
    check("PASS and unfiltered kept, QD2 dropped", kept == ["100", "300"])


def _thresholds(path: str, depth: int) -> None:
    with gzip.open(path, "wt") as fh:
        fh.write(f"#chrom\tstart\tend\tregion\t1X\t{depth}X\t20X\n")
        fh.write(f"1\t0\t1000\t{GENE}\t1000\t990\t900\n")


def test_callability_stale_cache() -> None:
    print("callability: a cached mosdepth run without the depth column:")
    sample = "Lena"
    cram = call.aligned_cram(sample)
    cram.parent.mkdir(parents=True, exist_ok=True)
    cram.write_bytes(b"invented")
    config.LOGS_DIR.mkdir(parents=True, exist_ok=True)
    out = config.CALLABLE_DIR / sample
    out.mkdir(parents=True, exist_ok=True)
    prefix = out / sample
    pathlib.Path(f"{prefix}.mosdepth.summary.txt").write_text(
        "chrom\tlength\tbases\tmean\tmin\tmax\ntotal\t100\t3000\t30\t0\t60\n")
    pathlib.Path(f"{prefix}.mosdepth.global.dist.txt").write_text(
        "total\t15\t0.95\ntotal\t10\t0.97\ntotal\t0\t1.0\n")
    with gzip.open(f"{prefix}.regions.bed.gz", "wt") as fh:
        fh.write(f"1\t0\t1000\t{GENE}\t30.0\n")
    _thresholds(f"{prefix}.thresholds.bed.gz", 10)
    panel = config.GENOMES_ROOT / "sv_genes_invented.tsv"
    panel.write_text(f"gene\tchrom\tstart\tend\tdisease\n{GENE}\t1\t1\t1000\tInvented\n")

    ran = []

    def fake_mosdepth(cmd, **_kw):
        ran.append(cmd)
        _thresholds(f"{prefix}.thresholds.bed.gz", 15)
        return subprocess.CompletedProcess(cmd, 0)

    with patched(config, CALLABLE_MIN_DEPTH=15, SV_GENES_FILE=panel), \
            patched(callability, _mosdepth_bin=lambda: "mosdepth"), \
            patched(callability.subprocess, run=fake_mosdepth):
        md = callability.callability_report(sample)
    check("mosdepth is re-run for the new depth", len(ran) == 1)
    check("the gene is not reported as 0.0%",
          f"| {GENE} |" in md.read_text() and "| 0.0% |" not in md.read_text())

    bad = config.GENOMES_ROOT / "no_column.thresholds.bed.gz"
    _thresholds(str(bad), 10)
    try:
        callability._gene_callability(bad, pathlib.Path(f"{prefix}.regions.bed.gz"), 15)
        loud = False
    except SystemExit:
        loud = True
    check("a missing column is an error, never 0%", loud)


def test_sv_partial_no_calls() -> None:
    print("sv: only genotypes with an alt allele are carried:")
    gts = {"0/0": False, "./.": False, ".": False, "0|0": False, "0/.": False, "./0": False,
           "0|.": False, ".|0": False, "0": False, "./1": True, "1/.": True, "0/1": True,
           "1|1": True, "1": True}
    body = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAdam\n"
    starts = {}
    for i, gt in enumerate(gts):
        pos = 1000 + i * 100
        starts[gt] = pos
        body += (f"1\t{pos}\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END={pos + 50};SVLEN=-50\t"
                 f"GT\t{gt}\n")
    path = config.GENOMES_ROOT / "invented_sv.vcf.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(body)
    kept = {s["start"] for s in sv._parse_sv_vcf(path)[0]}
    wrong = [gt for gt, carried in gts.items() if (starts[gt] in kept) != carried]
    check(f"carried status is right for every genotype (wrong: {wrong})", not wrong)


def test_viewing_as_card_is_highlighted() -> None:
    print("render: the selected person's card has a style:")
    cls = re.search(r"\[data-person\]'\)\.forEach\(el =>\s*el\.classList\.toggle\('(\w+)'",
                    render._PERSPECTIVE_SCRIPT).group(1)
    css = theme.STYLE + theme.SHELL_STYLE
    check("the class the script toggles is styled, for linked and unsequenced cards",
          f".people .card.{cls} " in css)


def test_genome_map_shows_y_and_mt() -> None:
    print("dashboard: the genome map covers Y and MT:")

    def finding(chrom, pos, gene):
        return {"chrom": chrom, "pos": pos, "gene": gene, "support": "", "low_ab": None,
                "copies": {"n": "1", "text": "one copy", "label": "One copy"}}
    html = dashboard.genome_map([finding("1", 1000, "GENEA"), finding("Y", 2_000_000, "GENEY"),
                                 finding("MT", 1234, "GENEM")])
    check("autosome pin still drawn", ">GENEA<" in html)
    check("Y pin drawn", ">GENEY<" in html)
    check("MT pin drawn", ">GENEM<" in html)
    check("MT is placed on its own length",
          'left:7.4%">GENEM<' in html)


def test_failed_shard_cancels_the_queue() -> None:
    print("call: a failed shard stops the queued windows:")
    total = 30
    sample = "Adam"
    config.GATK_JAR.parent.mkdir(parents=True, exist_ok=True)
    config.GATK_JAR.write_bytes(b"invented")
    cram = call.aligned_cram(sample)
    cram.parent.mkdir(parents=True, exist_ok=True)
    cram.write_bytes(b"invented")
    config.LOGS_DIR.mkdir(parents=True, exist_ok=True)
    tasks = [(f"1:{i}-{i}", 2, f"1_{i:09d}") for i in range(total)]
    ran = []

    def fake_run(cmd, **kw):
        kw["stdout"].close()
        ran.append(cmd)
        if len(ran) == 1:
            raise subprocess.CalledProcessError(1, "gatk")
        time.sleep(0.03)

    with patched(call, _scatter_tasks=lambda: tasks, _shard_costs=lambda labels: {},
                 _java_env=lambda: {}, run=fake_run), patched(config, HC_WORKERS=1):
        try:
            call.call(sample)
            raised = False
        except subprocess.CalledProcessError:
            raised = True
    check("the failure is raised", raised)
    check(f"the rest of the queue is not run ({len(ran)} of {total} started)", len(ran) < 5)


if __name__ == "__main__":
    test_delta_counts_only_compared_callsets()
    test_delta_t2t_cap_holds_for_a_single_callset()
    test_delta_mt_is_not_discordant()
    test_delta_vanished_variant()
    test_serve_percent_encoded_paths()
    test_ancestry_subset_is_pass_only()
    test_callability_stale_cache()
    test_sv_partial_no_calls()
    test_viewing_as_card_is_highlighted()
    test_genome_map_shows_y_and_mt()
    test_failed_shard_cancels_the_queue()
    print()
    if _failures:
        print(f"{len(_failures)} test(s) FAILED: {_failures}")
        sys.exit(1)
    print("ALL TESTS PASSED")
