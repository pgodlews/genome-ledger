#!/usr/bin/env python3
"""Regression tests for the PRS, traits, HLA and PGx review fixes. Invented data throughout.

    uv run python tests/test_review_prs_hla_pgx.py
"""

from __future__ import annotations

import email.message
import io
import json
import pathlib
import subprocess
import sys
import tempfile
import traceback
import urllib.request
from urllib.response import addinfourl

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import test_units  # noqa: E402  (also points the reports tree at a temp directory)
from test_units import check  # noqa: E402

from pipeline import config, forcecall, hla, hla_typing, pgx, prs, traits  # noqa: E402

TMP = pathlib.Path(tempfile.mkdtemp(prefix="genome_review_"))
HEADER = ("##fileformat=VCFv4.2\n##contig=<ID=1,length=400>\n##contig=<ID=2,length=400>\n"
          '##FILTER=<ID=PASS,Description="ok">\n'
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAdam\n")


def _reference() -> pathlib.Path:
    """Position 100-101 reads AT (A/AT is ambiguous); 200-201 reads AC (A/AT is an insertion)."""
    seq = list(("GGCCTTGGCC" * 40)[:400])
    seq[99], seq[100], seq[199], seq[200] = "A", "T", "A", "C"
    fa = TMP / "ref.fa"
    fa.write_text("".join(f">{c}\n{''.join(seq)}\n" for c in ("1", "2")))
    subprocess.run(["samtools", "faidx", str(fa)], check=True)
    return fa


def _vcf(name: str, rows: list[tuple[str, str, str]], pos: str = "100", norm: bool = False):
    raw, out = TMP / f"{name}.vcf", TMP / f"{name}.vcf.gz"
    raw.write_text(HEADER + "".join(f"1\t{pos}\t.\t{r}\t{a}\t50\tPASS\t.\tGT\t{g}\n"
                                    for r, a, g in rows))
    cmd = (["bcftools", "norm", "-f", str(config.REF_FASTA), "-m", "-both"] if norm
           else ["bcftools", "view"])
    subprocess.run([*cmd, "-Oz", "-o", str(out), str(raw)], check=True, capture_output=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(out)], check=True)
    return out


def _pvar(wd: pathlib.Path, chrom: str, rows: str) -> None:
    wd.mkdir(parents=True, exist_ok=True)
    (wd / f"chr{chrom}.pvar").write_text("#CHROM\tPOS\tID\tREF\tALT\n" + rows)


def test_trans_dq25_is_recognised() -> None:
    print("hla_typing._risk_calls (DQ2.5 in trans):")
    calls = hla_typing._risk_calls({"DQA1": ["DQA1*02:01:01", "DQA1*05:05:01"],
                                    "DQB1": ["DQB1*02:02:01", "DQB1*03:01:01"]})
    check("DQA1*05:05 + DQB1*02:02 is DQ2.5", ": DQ2.5 " in calls[1])
    alone = hla_typing._risk_calls({"DQA1": ["DQA1*02:01"], "DQB1": ["DQB1*02:02"]})
    check("DQ2.2 alone is not called DQ2.5", "neither" in alone[1])


def test_tag_negative_does_not_reassure() -> None:
    print("hla.hla_report (a tag negative is not an exclusion):")
    empty = _vcf("hla_empty", [])
    rules = TMP / "hla_tags.tsv"
    rules.write_text("trait\trsid\tchrom\tpos\tref\teffect\tinterp0\tinterp1\tinterp2\n"
                     "Celiac DQ2.5 haplotype\trs0000001\t1\t150\tG\tT\tabsent\thet\thom\n"
                     "Celiac DQ8 haplotype\trs0000002\t1\t160\tG\tT\tabsent\thet\thom\n")
    saved = config.HLA_DIR, config.HLA_FILE, hla._called_vcf, forcecall.read_genotypes
    config.HLA_DIR, config.HLA_FILE = TMP / "hla", rules
    hla._called_vcf = lambda s: empty
    hla.forcecall.read_genotypes = lambda s: {}
    try:
        text = pathlib.Path(hla.hla_report("Maria", force=True)).read_text()
    finally:
        config.HLA_DIR, config.HLA_FILE, hla._called_vcf, forcecall.read_genotypes = saved
    check("no 'very unlikely' claim", "very unlikely" not in text)
    check("no '>99%' claim", ">99%" not in text)
    check("says the negative does not rule it out", "not ruled out" in text)
    check("names what the tags cannot see", "trans" in text and "DQ2.2" in text)


def test_cram_file_uri() -> None:
    print("hla_typing._cram_ref_dir (file: URIs in @SQ UR):")
    saved_run, saved_ref = subprocess.run, config.REF_FASTA
    config.REF_FASTA = pathlib.Path("/refs/lena/ref.fa")

    def header(ur):
        out = f"@HD\tVN:1.6\n@SQ\tSN:1\tLN:400\tUR:{ur}\n"
        subprocess.run = lambda *a, **k: type("R", (), {"stdout": out})()
        return hla_typing._cram_ref_dir(pathlib.Path("/data/Jan.cram"))

    try:
        check("file:/// keeps only the directory", header("file:///refs/jan/ref.fa") == "/refs/jan")
        check("percent-encoding is decoded",
              header("file:///refs/jan%20old/ref.fa") == "/refs/jan old")
        check("file: URI of our own reference needs no mount",
              header("file:///refs/lena/ref.fa") is None)
        check("a bare path still works", header("/refs/jan/ref.fa") == "/refs/jan")
        check("a remote URL has nothing to mount", header("https://example.org/r/ref.fa") is None)
    finally:
        subprocess.run, config.REF_FASTA = saved_run, saved_ref


def test_pgx_flags_by_exact_value() -> None:
    print("pgx._write_summary (exact non-actionable values):")
    cases = {"GENE_A": ["Normal Metabolizer"], "GENE_B": ["Favorable response allele"],
             "GENE_C": ["ivacaftor non-responsive in CF patients"], "GENE_D": ["Indeterminate"],
             "GENE_E": ["Unfavorable response allele"], "GENE_F": ["Poor Metabolizer"],
             "GENE_G": ["Normal Metabolizer", "Intermediate Metabolizer"],
             "GENE_H": ["Abnormally Named New Phenotype"]}
    genes = {g: {"recommendationDiplotypes": [{"label": "*1/*2", "phenotypes": p}]}
             for g, p in cases.items()}
    genes["GENE_I"] = {"recommendationDiplotypes": []}
    rep, out = TMP / "pgx.json", TMP / "pgx.md"
    rep.write_text(json.dumps({"genes": genes}))
    pgx._write_summary("Lena", rep, out)
    table = out.read_text().split("## Actionable")[1].split("## All gene calls")[0]
    for g in ("GENE_E", "GENE_F", "GENE_G", "GENE_H"):
        check(f"{'; '.join(cases[g])!r} is flagged", f"| {g} |" in table)
    for g in ("GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_I"):
        check(f"{g} is not flagged", f"| {g} |" not in table)


def test_absent_indel_rows() -> None:
    print("prs._score_one (indel rows at sites the VCF does not carry):")
    other = _vcf("indel_other", [("G", "T", "0/1")], pos="300")
    for label, row in (("effect=A", ("1", "100", "A", "AT", 1.0)),
                       ("effect=AT", ("1", "100", "AT", "A", 1.0))):
        score, used, _c, contrib = prs._score_one([other], [row])
        check(f"ambiguous A/AT ({label}) is skipped, not guessed",
              score == 0.0 and used == 0 and not contrib)
    carrier = _vcf("indel_carrier", [("AT", "A", "0/1")])
    _s, used, _c, _x = prs._score_one([carrier], [("1", "100", "A", "AT", 1.0)])
    check("...and skipped for a carrier too (no selection on genotype)", used == 0)
    check("resolvable insertion: hom-ref has two copies of REF",
          prs._score_one([other], [("1", "200", "A", "AT", 1.0)])[0] == 2.0)
    check("resolvable insertion: hom-ref has none of the ALT",
          prs._score_one([other], [("1", "200", "AT", "A", 1.0)])[:2] == (0.0, 1))

    wd = TMP / "indel_no_pvars"
    check("FASTA supplement does not credit an ambiguous row",
          prs._ref_effect_supplement([("1", "100", "A", "AT", 1.0)], set(), wd) == (0.0, 0))
    check("FASTA supplement credits a resolved insertion REF",
          prs._ref_effect_supplement([("1", "200", "A", "AT", 1.0)], set(), wd) == (2.0, 1))

    print("prs._score_one (1000G states the REF of an ambiguous indel row):")
    deletion = {("1", "100"): ("AT", "A", "0.1", "0.1")}
    insertion = {("1", "100"): ("A", "AT", "0.1", "0.1")}
    check("panel deletion: hom-ref has none of the shorter allele",
          prs._score_one([other], [("1", "100", "A", "AT", 1.0)], deletion)[:2] == (0.0, 1))
    check("panel deletion: hom-ref has two copies of REF",
          prs._score_one([other], [("1", "100", "AT", "A", 1.0)], deletion)[:2] == (2.0, 1))
    check("panel insertion: hom-ref has two copies of the shorter allele",
          prs._score_one([other], [("1", "100", "A", "AT", 1.0)], insertion)[:2] == (2.0, 1))
    check("panel deletion: a het carrier scores one copy",
          prs._score_one([carrier], [("1", "100", "A", "AT", 1.0)], deletion)[:2] == (1.0, 1))
    merged = {("1", "100"): ("ATG", "AG,GTG", "0.1,0.2", "0.1,0.2")}
    check("panel record with a longer merged REF still resolves the row",
          prs._score_one([other], [("1", "100", "AT", "A", 1.0)], merged)[:2] == (2.0, 1))
    unrelated = {("1", "100"): ("A", "G", "0.1", "0.1")}
    check("a panel record for a different variant leaves the row skipped",
          prs._score_one([other], [("1", "100", "A", "AT", 1.0)], unrelated)[1] == 0)
    check("a panel record stating both orientations resolves to nothing",
          prs._panel_row_ref(("ATT", "AT,ATTT", "0.1,0.1", "0.1,0.1"), "A", ["AT"]) is None)


def test_colocated_genotypes() -> None:
    print("prs / traits (several records at one POS):")
    snv_first = _vcf("coloc_snv_first", [("AT", "GT,A", "1/2")], norm=True)
    check("PRS keeps the SNV allele when the indel record is last",
          prs._score_one([snv_first], [("1", "100", "G", "A", 1.0)])[0] == 1.0)
    indel_first = _vcf("coloc_indel_first", [("AT", "A,GT", "1/2")], norm=True)
    pair = traits._alleles_at(indel_first, "1", "100", "A")
    check("traits sees the SNV allele when the indel record is first", "G" in pair)
    check("...with the deletion's anchor base as the other copy", sorted(pair) == ["A", "G"])

    deletion = _vcf("coloc_del", [("AT", "A", "0/1")])
    check("an anchored deletion leaves the marker base on both copies (traits)",
          traits._alleles_at(deletion, "1", "100", "A") == ("A", "A"))
    check("...and in PRS",
          prs._score_one([deletion], [("1", "100", "A", "G", 1.0)])[0] == 2.0)
    snv = _vcf("coloc_snv", [("A", "G", "0/1")])
    check("a SNV is not read as a co-located insertion's REF allele",
          prs._score_one([snv], [("1", "200", "A", "AT", 1.0)])[1] == 1)
    three = traits._merge_records([("A", "G", "1/1"), ("A", "T", "0/1")])
    check("more ALT copies than chromosomes is unresolved", three[2] == "./.")
    check("a repeated record counts once",
          traits._merge_records([("A", "G", "0/1"), ("A", "G", "0/1")]) == ("A", "G", "0/1"))


def test_colocated_frequencies() -> None:
    print("prs._reference_freqs / _effect_freq (co-located 1000G records):")
    saved = subprocess.run, config.PGS_DIR, config.KG_VCF_LOCAL_DIR
    config.PGS_DIR = TMP / "pgs_coloc"
    config.PGS_DIR.mkdir()
    config.KG_VCF_LOCAL_DIR = config.PGS_DIR / "none"
    (config.PGS_DIR / "PGSADAM_1kg_af.tsv").write_text("1\t100\tAT\tA\t0.01\t0.01\n")
    subprocess.run = lambda *a, **k: type("R", (), {
        "stdout": "chr1\t100\tA\tG\t0.30\t0.20\nchr1\t100\tAT\tA\t0.01\t0.02\n",
        "returncode": 0, "stderr": ""})()
    try:
        rec = prs._reference_freqs("PGSADAM", [("1", "100", "G", "A", 1.0)])[("1", "100")]
        cached = prs._reference_freqs("PGSADAM", [("1", "100", "G", "A", 1.0)])[("1", "100")]
    finally:
        subprocess.run, config.PGS_DIR, config.KG_VCF_LOCAL_DIR = saved
    check("a cache from before the merge is not read back", "G" in rec[1])
    check("the merged record round-trips through the cache", cached == rec)
    check("SNV ALT frequency", prs._effect_freq("G", "A", rec) == (0.30, 0.20))
    fe, fa = prs._effect_freq("A", "G", rec)
    check("SNV REF frequency ignores the indel", abs(fe - 0.70) < 1e-9 and abs(fa - 0.80) < 1e-9)
    check("the deletion keeps its own frequency", prs._effect_freq("A", "AT", rec) == (0.01, 0.02))
    check("an indel-only record says nothing about a SNV row",
          prs._effect_freq("A", "G", ("AT", "A", "0.01", "0.01")) is None)


def test_missing_frequency() -> None:
    print("prs._effect_freq (missing AF):")
    check("REF effect allele with AF '.' is unresolved, not 1.0",
          prs._effect_freq("A", "G", ("A", "G", ".", ".")) is None)
    check("one missing ALT of two is unresolved too",
          prs._effect_freq("A", "G", ("A", "G,T", "0.1,.", "0.1,.")) is None)


def test_partial_pvars() -> None:
    print("prs._ref_alleles (partial 1000G pvars):")
    wd = TMP / "partial"
    _pvar(wd, "1", "1\t100\t1:100\tA\tG\n")
    rows = [("1", "100", "A", "G", 1.0), ("2", "100", "A", "G", 10.0)]
    check("both chromosomes are supplemented, not just the one with a pvar",
          prs._ref_effect_supplement(rows, set(), wd) == (22.0, 2))
    _pvar(wd, "2", "2\t100\t2:100\tC\tG\n")
    check("a complete pvar set stays authoritative",
          prs._ref_effect_supplement(rows, set(), wd) == (2.0, 1))


def test_population_map_download() -> None:
    print("prs._kg_population_map (atomic, bounded download):")
    body = b"sample\tpop\tsuper_pop\n" + b"".join(b"HG%03d\tGBR\tEUR\n" % i for i in range(20))
    state = {"short": True, "timeout": None}

    class Fake(urllib.request.BaseHandler):
        def fake_open(self, req):
            state["timeout"] = req.timeout
            head = email.message.Message()
            head["Content-Length"] = str(len(body))
            sent = body[:len(body) // 2] if state["short"] else body
            return addinfourl(io.BytesIO(sent), head, req.full_url, 200)

    saved = config.KG_PANEL_FILE, config.KG_PANEL_URL, config.PGS_DIR
    config.PGS_DIR = TMP / "pgs_map"
    config.KG_PANEL_FILE, config.KG_PANEL_URL = config.PGS_DIR / "pop.txt", "fake://x/pop"
    urllib.request.install_opener(urllib.request.build_opener(Fake))
    try:
        check("a short body gives no map", prs._kg_population_map() == {})
        check("...and leaves nothing under the final name", not config.KG_PANEL_FILE.exists())
        check("...nor a temp file", not list(config.PGS_DIR.glob("*.part")))
        check("the download has a timeout", state["timeout"] == config.PGS_DOWNLOAD_TIMEOUT)
        state["short"] = False
        check("the next run retries and gets every sample",
              len(prs._kg_population_map()) == 20)
    finally:
        urllib.request.install_opener(None)
        config.KG_PANEL_FILE, config.KG_PANEL_URL, config.PGS_DIR = saved


def test_duplicate_coordinate_not_credited() -> None:
    print("prs._ref_effect_supplement (coordinate plink2 drops on both sides):")
    wd = TMP / "dups"
    _pvar(wd, "1", "1\t100\t1:100\tA\tG\n1\t150\t1:150\tC\tT\n")
    rows = [("1", "100", "A", "G", 0.5), ("1", "100", "G", "A", 0.3),
            ("1", "150", "C", "T", 0.7)]
    supp, n = prs._ref_effect_supplement(rows, set(), wd)
    check("only the unique coordinate is credited", abs(supp - 1.4) < 1e-9 and n == 1)


def test_multi_record_site_is_not_hom_ref() -> None:
    print("prs._score_entry (site removed by --rm-dup exclude-all):")
    vcf = _vcf("multi", [("A", "G,C", "1/2")], norm=True)
    rows = [("1", "100", "A", "G", 1.0)]
    stub = TMP / "stub_pgen"
    stub.with_suffix(".pvar").write_text("#CHROM\tPOS\tID\tREF\tALT\n")
    saved = (prs._load_variants, prs._scoring_file, prs._score_plink2, prs._pgen,
             prs._reference_scores, prs._called_vcf, config.PGS_DIR, config.PRS_PLINK_MIN)
    config.PGS_DIR = TMP / "pgs_multi"
    _pvar(config.PGS_DIR / "PGSJAN_1kg", "1", "1\t100\t1:100\tA\tG\n")
    prs._load_variants = lambda p: rows
    prs._scoring_file = lambda pid: TMP / "unused"
    prs._score_plink2 = lambda s, pid, v, o: (0.0, 0, len(v))
    prs._pgen = lambda s: stub
    prs._reference_scores = lambda pid, v: None
    prs._called_vcf = lambda s: vcf
    config.PRS_PLINK_MIN = 0
    try:
        got = prs._score_entry("Jan", {"pgs_id": "PGSJAN", "trait": "t"}, [], TMP)
        check("a G/C sample has no REF copy to credit", got["raw"] == 0.0)
        rows[0] = ("1", "100", "G", "A", 1.0)
        got = prs._score_entry("Jan", {"pgs_id": "PGSJAN", "trait": "t"}, [], TMP)
        check("...and its real ALT copy is scored", got["raw"] == 1.0)
        check("the multi-record site is found",
              set(prs._multi_record_sites(vcf)) == {("1", "100")})
    finally:
        (prs._load_variants, prs._scoring_file, prs._score_plink2, prs._pgen,
         prs._reference_scores, prs._called_vcf, config.PGS_DIR, config.PRS_PLINK_MIN) = saved


TESTS = [test_trans_dq25_is_recognised, test_tag_negative_does_not_reassure, test_cram_file_uri,
         test_pgx_flags_by_exact_value, test_absent_indel_rows, test_colocated_genotypes,
         test_colocated_frequencies, test_missing_frequency, test_partial_pvars,
         test_population_map_download, test_duplicate_coordinate_not_credited,
         test_multi_record_site_is_not_hom_ref]


def main() -> int:
    config.REF_FASTA = _reference()
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
