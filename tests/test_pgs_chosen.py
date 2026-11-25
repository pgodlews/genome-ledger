#!/usr/bin/env python3
"""Tests for scores chosen from the PGS Catalog (pgs_catalog, pgs_chosen, pgs_web).
Invented data throughout: a three-person family, a 400-base reference, a three-variant score
for an invented trait, and catalog CSVs written here.

    uv run python tests/test_pgs_chosen.py
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import pathlib
import subprocess
import sys
import tempfile
import traceback

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import test_units  # noqa: E402
from test_units import check  # noqa: E402

from pipeline import (config, consent, pgs_catalog, pgs_chosen, prs,  # noqa: E402
                      reportpaths)

ROOT = pathlib.Path(tempfile.mkdtemp(prefix="genome_pgs_chosen_"))
TRAIT_A, TRAIT_B = "EFO_9999901", "EFO_9999902"
PID_A = "PGS999001"          # one trait, three variants, files present locally
PID_AB = "PGS999002"         # two traits, known only to the catalog snapshot
PID_B = "PGS999003"          # trait B only, files present locally
HEADER = ("##fileformat=VCFv4.2\n##contig=<ID=1,length=400>\n"
          '##FILTER=<ID=PASS,Description="ok">\n'
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n')


# --- fixtures ------------------------------------------------------------------------------
def _reference() -> None:
    seq = list(("GGCCTTGGCC" * 40)[:400])
    seq[99], seq[149], seq[249] = "A", "C", "A"          # positions 100, 150, 250
    config.REF_FASTA.parent.mkdir(parents=True, exist_ok=True)
    config.REF_FASTA.write_text(">1\n" + "".join(seq) + "\n")
    subprocess.run(["samtools", "faidx", str(config.REF_FASTA)], check=True)


def _vcf(sample: str, rows: list[tuple[str, str, str, str]]) -> pathlib.Path:
    from pipeline.call import called_vcf
    out = called_vcf(sample)
    out.parent.mkdir(parents=True, exist_ok=True)
    raw = out.with_name(f"{sample}.raw.vcf")
    raw.write_text(HEADER + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
                   + "".join(f"1\t{p}\t.\t{r}\t{a}\t50\tPASS\t.\tGT\t{g}\n" for p, r, a, g in rows))
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(out), str(raw)], check=True,
                   capture_output=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(out)], check=True)
    return out


def _scoring(pid: str, efos: str, traits: str, rows: list[tuple]) -> pathlib.Path:
    p = pgs_chosen.scoring_path(pid)
    p.parent.mkdir(parents=True, exist_ok=True)
    head = (f"###PGS CATALOG SCORING FILE\n#format_version=2.0\n#pgs_id={pid}\n"
            f"#trait_mapped={traits}\n#trait_efo={efos}\n#variants_number={len(rows)}\n"
            "#HmPOS_build=GRCh38\n#HmPOS_date=2026-01-01\n")
    body = "effect_allele\tother_allele\teffect_weight\thm_chr\thm_pos\n" + "".join(
        f"{ea}\t{oa}\t{w}\t{c}\t{pos}\n" for c, pos, ea, oa, w in rows)
    with gzip.open(p, "wt") as fh:
        fh.write(head + body)
    return p


def _af(pid: str, rows: list[tuple]) -> None:
    prs._kg_cache(pid).write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))


def _catalog() -> pathlib.Path:
    d = config.PGS_CATALOG_DIR / "2026-01-02"
    d.mkdir(parents=True, exist_ok=True)

    def write(name, header, rows):
        with (d / name).open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            w.writerows(rows)

    sc_h = ["Polygenic Score (PGS) ID", "PGS Name", "Reported Trait", "Mapped Trait(s) (EFO label)",
            "Mapped Trait(s) (EFO ID)", "PGS Development Method",
            "PGS Development Details/Relevant Parameters", "Original Genome Build",
            "Number of Variants", "Number of Interaction Terms", "Type of Variant Weight",
            "PGS Publication (PGP) ID", "Publication (PMID)", "Publication (doi)",
            "Score and results match the original publication",
            "Ancestry Distribution (%) - Source of Variant Associations (GWAS)",
            "Ancestry Distribution (%) - Score Development/Training",
            "Ancestry Distribution (%) - PGS Evaluation", "FTP link", "Release Date",
            "License/Terms of Use"]
    lic = "PGS obtained from the Catalog should be cited appropriately"

    def score(pid, label, efo, n, released):
        return [pid, f"{pid}_name", label, label, efo, "C+T", "", "GRCh38", str(n), "0", "beta",
                "PGP999001", "", "", "", "European:100", "European:100", "", "", released, lic]
    write(pgs_catalog.FILES["scores"], sc_h, [
        score(PID_A, "invented trait alpha", TRAIT_A, 3, "2024-01-01"),
        score(PID_AB, "invented trait alpha|invented trait beta", f"{TRAIT_A}|{TRAIT_B}", 3,
              "2025-01-01"),
        score(PID_B, "invented trait beta", TRAIT_B, 2, "2024-06-01")])
    write(pgs_catalog.FILES["traits"],
          ["Ontology Trait ID", "Ontology Trait Label", "Ontology Trait Description",
           "Ontology URL"],
          [[TRAIT_A, "invented trait alpha", "", ""], [TRAIT_B, "invented trait beta", "", ""]])
    write(pgs_catalog.FILES["publications"],
          ["PGS Publication/Study (PGP) ID", "First Author", "Title", "Journal Name",
           "Publication Date", "Release Date", "Authors", "digital object identifier (doi)",
           "PubMed ID (PMID)"],
          [["PGP999001", "Example A", "t", "J Invented", "2023-05-05", "", "", "", ""]])
    ss_h = ["PGS Sample Set (PSS)", "Polygenic Score (PGS) ID", "Number of Individuals",
            "Number of Cases", "Number of Controls", "Percent of Participants Who are Male",
            "Sample Age", "Broad Ancestry Category", "Ancestry (e.g. French, Chinese)",
            "Country of Recruitment", "Additional Ancestry Description",
            "Phenotype Definitions and Methods", "Followup Time", "GWAS Catalog Study ID (GCST...)",
            "Source PubMed ID (PMID)", "Source DOI", "Cohort(s)",
            "Additional Sample/Cohort Information"]

    def sset(pss, n, cases, anc="European"):
        return [pss, "", str(n), str(cases) if cases is not None else "", "", "", "", anc,
                "", "", "", "", "", "", "", "", "C1", ""]
    write(pgs_catalog.FILES["sample_sets"], ss_h, [
        sset("PSS1", 5000, 500), sset("PSS2", 6000, 600), sset("PSS3", 7000, 700),
        sset("PSS4", 8000, 8000),                       # impossible: as many cases as people
        sset("PSS5", 9000, 50),                         # too few cases
        sset("PSS6", 4000, 400, "African")])
    pm_h = ["PGS Performance Metric (PPM) ID", "Evaluated Score", "PGS Sample Set (PSS)",
            "PGS Publication (PGP) ID", "Reported Trait", "Covariates Included in the Model",
            "PGS Performance: Other Relevant Information", "Publication (PMID)",
            "Publication (doi)", "Hazard Ratio (HR)", "Odds Ratio (OR)", "Beta",
            "Area Under the Receiver-Operating Characteristic Curve (AUROC)",
            "Concordance Statistic (C-index)", "Other Metric(s)"]

    def perf(pid, pss, cov, auroc):
        return ["PPM", pid, pss, "PGP999001", "", cov, "", "", "", "", "", "", auroc, "", ""]
    write(pgs_catalog.FILES["performance"], pm_h, [
        perf(PID_A, "PSS1", "", "0.60"), perf(PID_A, "PSS2", "PCs 1-10", "0.70 [0.6,0.8]"),
        perf(PID_A, "PSS3", "", "0.80"), perf(PID_A, "PSS1", "age, sex", "0.95"),
        perf(PID_A, "PSS4", "", "0.99"), perf(PID_A, "PSS5", "", "0.98"),
        perf(PID_A, "PSS6", "", "0.97"),
        perf(PID_B, "PSS1", "0", "0.66"),
        perf(PID_AB, "PSS2", "age, sex, PCs 1-10", "0.75")])
    (d / pgs_catalog.CATEGORIES_FILE).write_text(json.dumps(
        [{"label": "Invented category", "efotraits": [{"id": TRAIT_A, "label": "alpha"},
                                                      {"id": TRAIT_B, "label": "beta"}]}]))
    pgs_catalog._write_manifest(d, "2026-01-02", 3)
    pgs_catalog._set_current("2026-01-02")
    pgs_catalog._loaded.clear()
    return d


def _family() -> None:
    config.FAMILY_FILE.parent.mkdir(parents=True, exist_ok=True)
    config.FAMILY_FILE.write_text(
        "id\tdisplay_name\tsex\tfather\tmother\tpartner\n"
        "Adam\tAdam\tmale\t\t\tMaria\nMaria\tMaria\tfemale\t\t\tAdam\n"
        "Lena\tLena\tfemale\tAdam\tMaria\t\n")
    config.SAMPLES_TSV.parent.mkdir(parents=True, exist_ok=True)
    config.SAMPLES_TSV.write_text("sample_id\nAdam\nAdam_vendor\nMaria\nLena\n")


def setup() -> None:
    config.set_genomes_root(ROOT)
    config.PGS_OFFLINE = False
    _reference()
    _family()
    _catalog()
    # Score A: 1:100 effect A is the REFERENCE (absent = 2 copies), 1:150 C>T, 1:250 A>G.
    _scoring(PID_A, TRAIT_A, "invented trait alpha",
             [("1", "100", "A", "G", 0.5), ("1", "150", "T", "C", 0.2), ("1", "250", "G", "A", 0.3)])
    _af(PID_A, [("1", "100", "A", "G", "0.2", "0.2"), ("1", "150", "C", "T", "0.3", "0.3"),
                ("1", "250", "A", "G", "0.4", "0.4")])
    _scoring(PID_B, TRAIT_B, "invented trait beta",
             [("1", "150", "T", "C", 1.0), ("1", "250", "G", "A", 1.0)])
    _af(PID_B, [("1", "150", "C", "T", "0.3", "0.3"), ("1", "250", "A", "G", "0.4", "0.4")])
    _vcf("Adam", [("150", "C", "T", "0/1"), ("250", "A", "G", "1/1")])
    _vcf("Adam_vendor", [("150", "C", "T", "0/1")])
    _vcf("Maria", [])
    _vcf("Lena", [("250", "A", "G", "0/1")])


# --- tests --------------------------------------------------------------------------------
def test_covariates() -> None:
    print("covariate classes:")
    for text, want in [("", "none"), ("PCs 1-10", "pcs"), ("PC1-PC10", "pcs"),
                       ("first 10 principal components", "pcs"), ("10 genetic PCs", "pcs"),
                       ("0", "unclear"), ("-", "unclear"), ("NA", "unclear"),
                       ("age, sex", "adjusted"), ("PCs 1-10, sex", "adjusted"),
                       ("age, sex, 10 PCs", "adjusted")]:
        check(f"{text!r} is {want}", pgs_catalog.covariate_class(text) == want)


def test_catalog_summary() -> None:
    print("catalog snapshot: search and the performance summary:")
    cat = pgs_catalog.load()
    check("search finds both invented traits", {t["efo"] for t in cat.traits("invented")}
          == {TRAIT_A, TRAIT_B})
    check("a category lists its traits", cat.categories()[0]["traits"] == 2)
    p = cat.performance(PID_A)
    check(f"median of covariate-free European cohorts only (got {p})",
          p.get("median") == 0.7 and p.get("cohorts") == 3 and p.get("basis") == "covariate-free")
    check("an evaluation with age and sex does not count when covariate-free ones exist",
          p.get("high") == 0.8)
    check("a field recorded as '0' is 'covariates not stated', not covariate-free",
          cat.performance(PID_B).get("basis") == "covariates not stated")
    check("only an adjusted model → labelled as such",
          cat.performance(PID_AB).get("basis") == "with covariates"
          and "only with age/sex" in pgs_catalog.describe_performance(cat.performance(PID_AB)))
    check("a newer score for the same trait is found", cat.newer(PID_A) == [PID_AB])
    check("the snapshot verifies", pgs_catalog.verify(pgs_catalog.current_dir()) == [])
    (pgs_catalog.current_dir() / pgs_catalog.FILES["traits"]).write_text("tampered\n")
    check("a changed file fails verification", pgs_catalog.verify(pgs_catalog.current_dir()) != [])
    _catalog()


def test_gate_before_download() -> None:
    print("no consent, nothing fetched or computed:")
    calls = []
    saved = prs._scoring_file
    prs._scoring_file = lambda pid: calls.append(pid) or (_ for _ in ()).throw(
        AssertionError("downloaded"))
    try:
        out = pgs_chosen.score([PID_AB, PID_A])
    finally:
        prs._scoring_file = saved
    check("nothing scored", out["scored"] == [] and out["failed"] == [])
    check("nothing downloaded", calls == [])
    check("every person × score is reported as skipped", len(out["skipped"]) == 6)
    check("no result directory was created", not any(config.PRS_CHOSEN_DIR.glob("*/*.json")))


def test_consent_and_score() -> None:
    print("consent → score → report:")
    pgs_chosen.set_score_consent("Adam", PID_A, True)
    check("a grant is recorded per person and trait",
          consent.has_consented("Adam", "pgs:" + TRAIT_A)
          and consent.has_consented("Adam_vendor", "pgs:" + TRAIT_A))
    m = pgs_chosen.matrix([PID_A, PID_AB])
    check("the matrix ticks Adam for A only", m["cells"]["Adam"] == {PID_A: True, PID_AB: False})
    check("no birth year or age anywhere in the consent ledger",
          "birth" not in config.CONSENT_FILE.read_text().lower())
    out = pgs_chosen.score([PID_A])
    check(f"only the consenting person is scored (got {out['scored']})",
          out["scored"] == [("Adam", PID_A)])
    r = json.loads(pgs_chosen.result_path("Adam", PID_A).read_text())
    # 1:100 absent, effect = REF → 2×0.5; 1:150 het → 0.2; 1:250 hom-alt → 2×0.3
    check(f"the raw score is the panel engine's (got {r['raw']})", abs(r["raw"] - 1.8) < 1e-9)
    check("a percentile came from the cached frequencies", r["pe"] is not None)
    check("the result is stamped", r["stamp"]["hm_date"] == "2026-01-01"
          and r["stamp"]["scoring_file"].startswith("sha256:"))
    md = reportpaths.chosen_report("prs", "Adam")
    check("the report is filed under genome/<person>/chosen/",
          md == config.REPORTS_DIR / "md" / "genome" / "Adam" / "chosen" / "prs_Adam.md"
          and md.exists())
    check("its front matter names it a chosen report",
          "report: chosen_prs" in md.read_text() and "callset: Adam" in md.read_text())
    check("…and it renders to the HTML tree at the same path",
          reportpaths.html_for(md).exists()
          and reportpaths.html_for(md).parent.name == "chosen")
    check("the curated PRS report is not touched",
          not reportpaths.genome_report("prs", "Adam").exists())
    check("the chosen folder is not mistaken for a snapshot",
          reportpaths.md_snapshots() == [])
    check("the report lists the covariate-free median",
          "0.70 median of 3 cohorts" in md.read_text())

    pgs_chosen.set_score_consent("Adam", PID_AB, True)
    check("a two-trait score needs both traits: consenting via the score grants both",
          consent.has_consented("Adam", "pgs:" + TRAIT_B))
    consent.set_consent("Adam", "pgs:" + TRAIT_B, False)
    check("withdrawing one of its traits withdraws the score",
          not pgs_chosen.may_score("Adam", [TRAIT_A, TRAIT_B])
          and pgs_chosen.may_score("Adam", [TRAIT_A]))

    pgs_chosen.score([PID_A], ["Adam"], all_callsets=True)
    check("--all-callsets scores the person's other callsets under the same grant",
          pgs_chosen.result_path("Adam_vendor", PID_A).exists())


def test_withdrawal_deletes() -> None:
    print("withdrawal deletes the data and rebuilds the report:")
    pgs_chosen.set_score_consent("Adam", PID_B, True)
    pgs_chosen.score([PID_B], ["Adam"])
    md = reportpaths.chosen_report("prs", "Adam")
    check("two traits listed", PID_A in md.read_text() and PID_B in md.read_text())
    consent.set_consent("Adam", "pgs:" + TRAIT_B, False)
    check("the withdrawn trait's result is deleted",
          not pgs_chosen.result_path("Adam", PID_B).exists())
    check("the report is rebuilt without it",
          md.exists() and PID_B not in md.read_text() and PID_A in md.read_text())
    consent.set_consent("Adam", "pgs:" + TRAIT_A, False)
    check("the last withdrawal removes the results of every callset of the person",
          not pgs_chosen.result_path("Adam", PID_A).exists()
          and not pgs_chosen.result_path("Adam_vendor", PID_A).exists())
    check("…and the report in every format",
          not any(f.exists() for f in pgs_chosen.report_outputs("Adam")))


def test_stale() -> None:
    print("stamps and re-scoring:")
    pgs_chosen.set_score_consent("Lena", PID_A, True)
    pgs_chosen.score([PID_A], ["Lena"])
    check("fresh result is current", all(s["current"] for s in pgs_chosen.status()))
    _vcf("Lena", [("250", "A", "G", "1/1")])
    st = [s for s in pgs_chosen.status() if s["callset"] == "Lena"][0]
    check(f"new genotypes make it stale (got {st})", not st["current"] and "vcf" in st["why"])
    out = pgs_chosen.rescore_stale()
    r = json.loads(pgs_chosen.result_path("Lena", PID_A).read_text())
    check("--stale recomputes it", out["scored"] == [("Lena", PID_A)]
          and abs(r["raw"] - 1.6) < 1e-9 and all(s["current"] for s in pgs_chosen.status()))


def test_bundle_roundtrip() -> None:
    print("bundles carry a score to a machine that never fetched it:")
    carry = ROOT / "carry"
    pgs_chosen.write_bundle(PID_A, carry)
    b, problems = pgs_chosen.verify_bundle(carry / "bundles" / PID_A)
    check("the bundle verifies", problems == [] and b["efo_ids"] == [TRAIT_A])
    for f in [pgs_chosen.scoring_path(PID_A), prs._kg_cache(PID_A)]:
        f.unlink()
    check("gone here", not pgs_chosen.ready(PID_A)[0])
    out = pgs_chosen.import_dir(carry)
    check("imported and ready again", out["bundles"] == [PID_A] and pgs_chosen.ready(PID_A)[0])
    target = carry / "bundles" / PID_A / f"{PID_A}_1kg_af.v2.tsv"
    target.write_text(target.read_text().replace("0.2", "0.9"))
    out = pgs_chosen.import_dir(carry)
    check("a tampered bundle is refused whole", out["bundles"] == [] and len(out["refused"]) == 1)
    bad = carry / "bundles" / "PGS999999"
    bad.mkdir(parents=True)
    (bad / "bundle.json").write_text(json.dumps({"format": 1, "pgs_id": "PGS999999",
                                                 "files": {"../escape.txt": "0"}}))
    out = pgs_chosen.import_dir(carry)
    check("a path outside the bundle is refused", any(p == "PGS999999" for p, _ in out["refused"]))


def test_offline() -> None:
    print("offline: a missing score is reported, not fetched:")
    config.PGS_OFFLINE = True
    try:
        pgs_chosen.set_score_consent("Maria", PID_AB, True)
        out = pgs_chosen.score([PID_AB], ["Maria"])
        check("failed with an offline reason", out["failed"]
              and "offline" in out["failed"][0][2])
        try:
            pgs_chosen.fetch([PID_AB])
            check("fetch refuses offline", False)
        except SystemExit:
            check("fetch refuses offline", True)
    finally:
        config.PGS_OFFLINE = False


def test_ref_table() -> None:
    print("the compact reference-allele table stands in for the pvars:")
    wd = config.PGS_DIR / "PGS999100_1kg"
    wd.mkdir(parents=True, exist_ok=True)
    (wd / "chr1.pvar").write_text("#CHROM\tPOS\tID\tREF\tALT\n1\t100\t.\tA\tG\n1\t150\t.\tc\tT\n")
    variants = [("1", "100", "A", "G", 1.0), ("1", "150", "T", "C", 1.0), ("1", "300", "G", "A", 1.0)]
    from_pvar = prs._ref_alleles(variants, wd)
    path = prs.export_ref_table("PGS999100", variants[:2])
    check("the table is written from complete pvars", path is not None and path.exists())
    (wd / "chr1.pvar").unlink()
    check("…and read back identically without them",
          prs._ref_alleles(variants, wd) == from_pvar == {("1", "100"): "A", ("1", "150"): "C"})
    check("no table from incomplete pvars", prs.export_ref_table("PGS999101", variants) is None)


def test_remote_header() -> None:
    print("upstream header read from the first chunk only:")
    import urllib.request
    body = gzip.compress(("#pgs_id=PGS999001\n#HmPOS_date=2026-03-03\nchr\tpos\n"
                          + "1\t1\n" * 100000).encode())

    class Resp(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False
    saved = urllib.request.urlopen
    urllib.request.urlopen = lambda req, timeout=0: Resp(body)
    try:
        h = pgs_chosen.remote_header(PID_A)
    finally:
        urllib.request.urlopen = saved
    check("HmPOS_date parsed", h.get("HmPOS_date") == "2026-03-03")


def test_web() -> None:
    print("picker API: reads open, writes need the token, a loopback client and our origin:")
    import threading
    import urllib.error
    import urllib.request
    from http.server import ThreadingHTTPServer

    from pipeline import pgs_web, serve
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), serve._Handler)
    port = httpd.server_address[1]
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{port}"

    def get(path):
        with urllib.request.urlopen(base + path, timeout=10) as r:
            return r.status, r.read().decode()

    def post(path, payload, token=None, origin=None):
        req = urllib.request.Request(base + path, data=json.dumps(payload).encode(),
                                     method="POST", headers={"Content-Type": "application/json"})
        if token:
            req.add_header("X-PGS-Token", token)
        if origin:
            req.add_header("Origin", origin)
        try:
            with urllib.request.urlopen(req, timeout=10) as r:
                return r.status, r.read().decode()
        except urllib.error.HTTPError as e:
            return e.code, e.read().decode()
    try:
        code, page = get("/pgs/")
        check("the picker page is served", code == 200 and "PGS" in page)
        check("the page carries the session token", pgs_web.TOKEN in page)
        import re
        import shutil
        node = shutil.which("node")
        if node:
            js = ROOT / "picker.js"
            js.write_text(re.search(r"<script>(.*)</script>", page, re.S).group(1))
            ok = subprocess.run([node, "--check", str(js)], capture_output=True)
            check("the page's script parses (node --check)", ok.returncode == 0)
        code, body = get("/pgs/api/traits?q=invented")
        check("traits search", code == 200 and TRAIT_A in body)
        code, body = get(f"/pgs/api/scores?efo={TRAIT_A}")
        check("a trait's scores", code == 200 and PID_A in body)
        code, body = get(f"/pgs/api/matrix?pids={PID_B}")
        data = json.loads(body)
        check("the matrix lists the family and the warning",
              code == 200 and {p["id"] for p in data["people"]} == {"Adam", "Maria", "Lena"}
              and "child" in data["warning"])
        grant = {"person": "Lena", "pgs_id": PID_B, "granted": True}
        check("Lena starts without consent for trait B",
              not consent.has_consented("Lena", "pgs:" + TRAIT_B))
        code, _ = post("/pgs/api/consent", grant)
        check("a write without the token is refused", code == 403
              and not consent.has_consented("Lena", "pgs:" + TRAIT_B))
        code, _ = post("/pgs/api/consent", grant, pgs_web.TOKEN, origin="http://evil.example")
        check("a write from another origin is refused", code == 403)
        code, _ = post("/pgs/api/consent", grant, pgs_web.TOKEN, origin=base)
        check("a write with the token from this page is accepted", code == 200
              and consent.has_consented("Lena", "pgs:" + TRAIT_B))
        code, body = post("/pgs/api/score", {"pgs_ids": [PID_B], "people": ["Lena"]},
                          pgs_web.TOKEN, origin=base)
        check("a scoring job starts", code == 202)
        import time
        for _ in range(100):
            if not json.loads(get("/pgs/api/job")[1])["running"]:
                break
            time.sleep(0.1)
        code, body = get("/pgs/api/results")
        check("…and its result is listed with a link to the report",
              any(r["callset"] == "Lena" and r["pgs_id"] == PID_B and r["report"]
                  for r in json.loads(body)))
        check("the reports tree is still the only file area served",
              serve._ALLOWED == ("reports/html",))
    finally:
        httpd.shutdown()


def test_unsupported() -> None:
    print("scores the additive engine cannot compute are refused, not mis-scored:")
    check("plain odds ratios", pgs_catalog.unsupported("Odds Ratio over expected risk")
          is not None)
    check("log odds are fine", pgs_catalog.unsupported("log(OR)") is None
          and pgs_catalog.unsupported("beta") is None and pgs_catalog.unsupported("NR") is None
          and pgs_catalog.unsupported("ln(OR)") is None)
    check("interaction terms", pgs_catalog.unsupported("beta", 2) is not None)
    pid = "PGS999004"
    p = pgs_chosen.scoring_path(pid)
    with gzip.open(p, "wt") as fh:
        fh.write(f"#pgs_id={pid}\n#trait_efo={TRAIT_A}\n#trait_mapped=invented trait alpha\n"
                 "#weight_type=beta\neffect_allele\tother_allele\teffect_weight\thm_chr\t"
                 "hm_pos\tdosage_0_weight\tdosage_1_weight\tdosage_2_weight\n"
                 "A\tG\t0.1\t1\t100\t0\t0.1\t0.3\n")
    check("per-dosage weight columns are found in the file",
          "dosage_0_weight" in (pgs_chosen.engine_problem(pid) or ""))
    p2 = pgs_chosen.scoring_path("PGS999005")
    with gzip.open(p2, "wt") as fh:
        fh.write("#pgs_id=PGS999005\n#trait_efo=" + TRAIT_A + "\n#weight_type=beta\n"
                 "effect_allele\tother_allele\teffect_weight\thm_chr\thm_pos\tis_interaction"
                 "\tis_haplotype\tdosage_0_weight\n"
                 "A\tG\t0.1\t1\t100\tFALSE\tFALSE\t\nT\tC\t0.2\t1\t150\tFALSE\t\t\n")
    check("the columns present but every row FALSE or blank: an ordinary additive score",
          pgs_chosen.engine_problem("PGS999005") is None)
    p3 = pgs_chosen.scoring_path("PGS999006")
    with gzip.open(p3, "wt") as fh:
        fh.write("#pgs_id=PGS999006\n#trait_efo=" + TRAIT_A + "\n#weight_type=beta\n"
                 "effect_allele\tother_allele\teffect_weight\thm_chr\thm_pos\tis_haplotype\n"
                 "A\tG\t0.1\t1\t100\tFALSE\nT\tC\t0.2\t1\t150\tTRUE\n")
    check("one TRUE haplotype row is refused", "is_haplotype" in (pgs_chosen.engine_problem("PGS999006") or ""))
    pgs_chosen.set_score_consent("Maria", pid, True)
    out = pgs_chosen.score([pid], ["Maria"])
    check("…and refused at scoring with that reason", out["scored"] == []
          and "does not model" in out["failed"][0][2])
    check("nothing was written", not pgs_chosen.result_path("Maria", pid).exists())


def test_select_and_resume() -> None:
    print("selecting ids for a long fetch, and resuming one:")
    cat = pgs_catalog.load()
    check("by category", set(pgs_catalog.select(cat, ["Invented category"])) == {PID_A, PID_AB, PID_B})
    check("european evidence keeps all three (each has a usable European evaluation)",
          set(pgs_catalog.select(cat, None, "european")) == {PID_A, PID_AB, PID_B})
    check("covariate-free keeps only the score with such an evaluation",
          pgs_catalog.select(cat, None, "covariate-free") == [PID_A])
    try:
        pgs_catalog.select(cat, ["No such category"])
        check("an unknown category is refused", False)
    except SystemExit:
        check("an unknown category is refused", True)
    carry = ROOT / "carry2"
    pgs_chosen.write_bundle(PID_A, carry)
    calls = []
    saved = prs._scoring_file
    prs._scoring_file = lambda pid: calls.append(pid) or saved(pid)
    try:
        done = pgs_chosen.fetch([PID_A], to=carry)
    finally:
        prs._scoring_file = saved
    check("an intact bundle already written is skipped, not rebuilt",
          done == [PID_A] and calls == [])
    order = []
    saved = (pgs_chosen._reference_batches, pgs_chosen.write_bundle)
    pgs_chosen._reference_batches = lambda part, to, refresh, size: order.append(("ref", list(part)))
    pgs_chosen.write_bundle = lambda pid, to, link=False: order.append(("bundle", pid)) or saved[1](pid, to, link)
    try:
        many = [PID_A, PID_B] * 5                       # 10 ids, batch 2 → windows of 8
        pgs_chosen.fetch(many, ROOT / "carry5", batch=2)
    finally:
        pgs_chosen._reference_batches, pgs_chosen.write_bundle = saved
    firsts = [i for i, (kind, _x) in enumerate(order) if kind == "bundle"]
    refs = [i for i, (kind, _x) in enumerate(order) if kind == "ref"]
    check(f"batched fetch writes bundles window by window (order {[k for k, _ in order]})",
          len(refs) == 2 and firsts and firsts[0] < refs[1])
    linked = pgs_chosen.write_bundle(PID_A, ROOT / "carry3", link=True)
    check("--link hard-links bundle files to the cache",
          (linked / f"{PID_A}_hmPOS_GRCh38.txt.gz").stat().st_ino
          == pgs_chosen.scoring_path(PID_A).stat().st_ino
          and not pgs_chosen.verify_bundle(linked)[1])


def test_batch_equals_per_score() -> None:
    print("one pass for a batch gives each score exactly its per-score reference data:")
    import random
    import shutil
    if not (list(config.PLINK2_DIR.glob("plink2")) or shutil.which("plink2")):
        print("  SKIP  plink2 not available here")
        return
    rnd = random.Random(11)
    kg = config.KG_VCF_LOCAL_DIR
    kg.mkdir(parents=True, exist_ok=True)
    samples = [f"HG{i:05d}" for i in range(40)]
    config.KG_PANEL_FILE.parent.mkdir(parents=True, exist_ok=True)
    config.KG_PANEL_FILE.write_text("SampleID\tSex\tBiosample\tPopulation\tPopName\t"
                                    "Superpopulation\n" + "".join(
                                        f"{s}\t1\tx\tP\tp\t{'EUR' if i % 2 else 'AFR'}\n"
                                        for i, s in enumerate(samples)))
    sites: dict[str, list] = {}
    for chrom in ("1", "2"):
        rows, pos = [], 1000
        for i in range(300):
            pos += rnd.randint(5, 50)
            ref, alt = rnd.sample("ACGT", 2)
            if i == 50:                       # a multi-allelic record: --max-alleles 2 drops it
                alt = alt + "," + next(b for b in "ACGT" if b not in (ref, alt))
            gts = ["|".join(rnd.choice("01") for _ in range(2)) for _ in samples]
            rows.append((pos, ref, alt, gts))
            if i == 80:                       # two records at one position: --rm-dup drops both
                rows.append((pos, ref, next(b for b in "ACGT" if b not in (ref, alt)), gts))
        sites[chrom] = rows
        raw = kg / f"chr{chrom}.raw.vcf"
        raw.write_text("##fileformat=VCFv4.2\n" + f"##contig=<ID=chr{chrom}>\n"
                       + '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
                       + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                       + "\t".join(samples) + "\n"
                       + "".join(f"chr{chrom}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t"
                                 + "\t".join(g) + "\n" for p, r, a, g in rows))
        out = kg / f"chr{chrom}.vcf.gz"
        subprocess.run(["bcftools", "view", "-Oz", "-o", str(out), str(raw)], check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", str(out)], check=True)

    def score(n, offset):
        out = []
        for chrom, rows in sites.items():
            for p, ref, alt, _g in rows[offset:offset + n]:
                a1 = alt.split(",")[0]
                out.append((chrom, str(p), rnd.choice([ref, a1]), a1, round(rnd.uniform(-1, 1), 4)))
        return out
    scores = {"PGS999201": score(150, 0), "PGS999202": score(120, 40),   # overlaps the first
              "PGS999203": score(60, 220),
              # chr2 holds a single variant at a position 1000G does not have: plink2 finds
              # nothing to score there, which must count as 0, not as a failed chromosome.
              "PGS999204": [v for v in score(50, 0) if v[0] == "1"]
                           + [("2", "999999", "A", "G", 0.5)]}
    per = {}
    for pid, vs in scores.items():
        shutil.rmtree(config.PGS_DIR / f"{pid}_1kg", ignore_errors=True)
        ref = prs._reference_scores(pid, vs)
        table = prs.export_ref_table(pid, vs)
        per[pid] = ({iid: sc for iid, sc in
                     ((l.split("\t")[0], float(l.split("\t")[2])) for l in
                      (config.PGS_DIR / f"{pid}_1kg" / "scores_all.tsv").read_text().splitlines())},
                    sorted(gzip.open(table, "rt").read().splitlines()))
        check(f"{pid}: per-score reference built", ref is not None and len(per[pid][0]) == 40)
        shutil.rmtree(config.PGS_DIR / f"{pid}_1kg")
    done = prs.reference_scores_batch(scores, config.PGS_DIR / ".batch-test")
    check("the batch completes every score", sorted(done) == sorted(scores))
    for pid in scores:
        d = config.PGS_DIR / f"{pid}_1kg"
        got = {l.split("\t")[0]: float(l.split("\t")[2])
               for l in (d / "scores_all.tsv").read_text().splitlines()}
        same = got.keys() == per[pid][0].keys() and all(
            abs(got[k] - per[pid][0][k]) < 1e-9 for k in got)
        check(f"{pid}: reference scores equal the per-score build", same)
        check(f"{pid}: REF table equals the per-score build",
              sorted(gzip.open(prs.ref_table(d), "rt").read().splitlines()) == per[pid][1])


def test_batching_is_machine_aware() -> None:
    print("batches close at the score, row and distinct-position limits; limits fit the machine:")
    seen = []
    saved = (prs.reference_scores_batch, config.PRS_BATCH_MAX_ROWS, config.PRS_BATCH_MAX_UNION)
    prs.reference_scores_batch = lambda batch, wd: seen.append(sorted(batch)) or list(batch)
    try:
        config.PRS_BATCH_MAX_ROWS, config.PRS_BATCH_MAX_UNION = 10**9, 10**9
        site = lambda c, p: (c, str(p), "A", "G", 1.0)
        shared = [site("1", p) for p in range(1, 1001)]
        items = [(f"S{i}", shared) for i in range(5)]           # five scores on one site set
        built = prs.prebuild_references(iter(items), size=3)
        check(f"closed at 3 scores (got {seen})", seen == [["S0", "S1", "S2"], ["S3", "S4"]]
              and len(built) == 5)
        seen.clear()
        config.PRS_BATCH_MAX_UNION = 1_000_000                   # the floor batch_limits keeps
        big = lambda k: [site("2", k * 10**6 + p) for p in range(600_000)]
        prs.prebuild_references(iter([("B1", big(1)), ("B2", big(2)), ("B3", big(3))]), size=40)
        check(f"unrelated large scores close the batch at the union cap (got {seen})",
              seen == [["B1"], ["B2"], ["B3"]])
    finally:
        prs.reference_scores_batch, config.PRS_BATCH_MAX_ROWS, config.PRS_BATCH_MAX_UNION = saved
    rows, union = prs.batch_limits(config.PGS_DIR)
    check(f"limits are positive and no larger than the configured caps ({rows:,}, {union:,})",
          1_000_000 <= rows <= config.PRS_BATCH_MAX_ROWS
          and 1_000_000 <= union <= config.PRS_BATCH_MAX_UNION)
    share = prs._plink2_share()
    threads = int(share[share.index("--threads") + 1])
    check(f"each plink2 gets its share of the cores ({share})",
          1 <= threads <= max(1, (__import__("os").cpu_count() or 2) // max(1, config.PRS_REF_WORKERS)))

    print("the panel stage prebuilds only genome-wide scores that lack reference data:")
    got = []
    saved = (prs.prebuild_references, prs._scoring_file, prs._load_variants)
    prs.prebuild_references = lambda items, size=None: got.extend(p for p, _v in items) or []
    prs._scoring_file = lambda pid: pid
    prs._load_variants = lambda pid: [("1", str(i), "A", "G", 1.0)
                                      for i in range(config.PRS_PLINK_MIN + (1 if pid != "SMALL" else -10))]
    try:
        (config.PGS_DIR / "HAVE_1kg").mkdir(parents=True, exist_ok=True)
        (config.PGS_DIR / "HAVE_1kg" / "scores_all.tsv").write_text("x\ty\t1\n")
        prs._prebuild_panel([{"pgs_id": "NEED", "n_variants": ""},
                             {"pgs_id": "HAVE", "n_variants": "900000"},
                             {"pgs_id": "SMALL", "n_variants": "40"}])
    finally:
        prs.prebuild_references, prs._scoring_file, prs._load_variants = saved
    check(f"only the genome-wide one without reference data (got {got})", got == ["NEED"])


def test_local_zero_rows_are_an_answer() -> None:
    print("a chromosome with none of the score's positions in a LOCAL 1000G copy is cached:")
    kg = config.KG_VCF_LOCAL_DIR
    kg.mkdir(parents=True, exist_ok=True)
    for chrom, rows in (("1", "chr1\t100\t.\tA\tG\t.\tPASS\tAF=0.3;AF_EUR=0.2\n"),
                        ("2", "")):
        raw = kg / f"af{chrom}.vcf"
        raw.write_text("##fileformat=VCFv4.2\n" + f"##contig=<ID=chr{chrom}>\n"
                       '##INFO=<ID=AF,Number=A,Type=Float,Description="a">\n'
                       '##INFO=<ID=AF_EUR,Number=A,Type=Float,Description="e">\n'
                       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + rows)
        out = kg / f"chr{chrom}.vcf.gz"
        subprocess.run(["bcftools", "view", "-Oz", "-o", str(out), str(raw)], check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", str(out)], check=True)
    pid = "PGS999301"
    prs._kg_cache(pid).unlink(missing_ok=True)
    freqs = prs._reference_freqs(pid, [("1", "100", "G", "A", 1.0), ("2", "500", "T", "C", 1.0)])
    check("the position that exists is found", ("1", "100") in freqs)
    check("…and the cache is written although chr2 had none of the positions",
          prs._kg_cache(pid).exists() and prs._kg_cache(pid).stat().st_size > 0)


def test_update() -> None:
    print("pgs update: nothing unless the release changed; builds what the store lacks:")
    from pipeline import pgs_update
    store = ROOT / "store"
    (store / "2026-01-02" / "bundles").mkdir(parents=True, exist_ok=True)
    saved = (pgs_update.remote_release, pgs_catalog.sync, pgs_chosen.fetch, pgs_update.notify.notify,
             pgs_chosen.remote_header)
    pgs_chosen.remote_header = lambda pid: {}            # no network in tests
    calls, notes = [], []
    pgs_update.notify.notify = lambda m: notes.append(m) or True
    try:
        pgs_update.remote_release = lambda: {"date": "2026-01-02"}
        out = pgs_update.update(store)
        check("an unchanged release does nothing and notifies nothing",
              not out["changed"] and notes == [])
        pgs_update.remote_release = lambda: {"date": "2026-02-01"}
        pgs_catalog.sync = lambda *a, **k: pgs_catalog.current_dir()

        def fake_fetch(pids, to=None, refresh=False, prune=False, batch=1, link=False):
            calls.append(list(pids))
            done = []
            for p in pids:
                if p == PID_B:
                    pgs_chosen.fetch_errors[p] = "cannot be scored here: has is_haplotype rows, which this engine does not model"
                    continue
                (Path(to) / "bundles" / p).mkdir(parents=True, exist_ok=True)
                (Path(to) / "bundles" / p / "bundle.json").write_text("{}")
                done.append(p)
            return done
        from pathlib import Path
        pgs_chosen.fetch = fake_fetch
        out = pgs_update.update(store)
        check(f"a new release builds every computable score the store lacks (got {calls})",
              out["changed"] and sorted(calls[0]) == sorted([PID_A, PID_AB, PID_B]))
        check("a haplotype refusal is remembered in refused.tsv",
              PID_B in pgs_update.refused(store))
        check("one notification, naming the release and the counts",
              len(notes) == 1 and "2026-02-01" in notes[0] and "2 new score bundle(s)" in notes[0]
              and "1 refused" in notes[0])
        calls.clear()
        out = pgs_update.update(store, force=True)
        check("next time neither built nor refused ones are retried", calls == [])
        for p in (PID_A, PID_AB):
            import shutil
            shutil.rmtree(store / "2026-02-01" / "bundles" / p)
        calls.clear(); notes.clear()
        out = pgs_update.update(store, max_new=1, force=True)
        check("more new scores than --max-new: nothing built, and it says so",
              out["too_many"] and calls == [] and "too many" in notes[0])
    finally:
        (pgs_update.remote_release, pgs_catalog.sync, pgs_chosen.fetch, pgs_update.notify.notify,
         pgs_chosen.remote_header) = saved


TESTS = [test_covariates, test_catalog_summary, test_gate_before_download, test_consent_and_score,
         test_withdrawal_deletes, test_stale, test_bundle_roundtrip, test_offline,
         test_ref_table, test_remote_header, test_unsupported, test_select_and_resume, test_batch_equals_per_score, test_batching_is_machine_aware, test_local_zero_rows_are_an_answer, test_update, test_web]


def main() -> int:
    setup()
    failed = []
    for test in TESTS:
        before = len(test_units._failures)
        try:
            test()
        except Exception:  # noqa: BLE001 — one crash must not hide the other results
            traceback.print_exc()
            test_units._failures.append(test.__name__)
        if len(test_units._failures) > before:
            failed.append(test.__name__)
    print()
    if failed:
        print(f"{len(failed)} of {len(TESTS)} test(s) FAILED: {failed}")
        return 1
    print(f"ALL {len(TESTS)} TESTS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
