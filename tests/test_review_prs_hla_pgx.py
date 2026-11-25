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

    print("hla_typing._risk_calls (an untyped locus is not a negative):")
    untyped = hla_typing._risk_calls({"DQA1": ["DQA1*01:01"], "DQB1": ["DQB1*05:01"]})
    check("no HLA-B call is not reported as B27 absent",
          "not determined" in untyped[0] and "absent" not in untyped[0])
    check("typed HLA-B without B*27 is still absent",
          "absent" in hla_typing._risk_calls({"B": ["B*07:02"]})[0])
    check("no DQB1 call leaves celiac undetermined",
          "not determined" in hla_typing._risk_calls({"DQA1": ["DQA1*05:01"]})[1])
    check("a DQ2 beta chain without DQA1 leaves celiac undetermined",
          "not determined" in hla_typing._risk_calls({"DQB1": ["DQB1*02:01"]})[1])
    check("DQB1 with neither risk chain excludes both without DQA1",
          "neither" in hla_typing._risk_calls({"DQB1": ["DQB1*05:01"]})[1])
    check("DQ8 needs only DQB1",
          ": DQ8 " in hla_typing._risk_calls({"DQB1": ["DQB1*03:02"]})[1])

    print("hla_typing._risk_calls (typing and tag SNPs disagree):")
    other = {"B": ["B*07:02"], "DQA1": ["DQA1*01:01"], "DQB1": ["DQB1*05:01"]}
    near = {"DQA1": ["DQA1*05:01"], "DQB1": ["DQB1*02:03"]}
    trans = {"DQA1": ["DQA1*05:05"], "DQB1": ["DQB1*02:02"]}
    b27 = {"B": ["B*27:05"]}
    dq8 = {"DQB1": ["DQB1*03:02"]}
    near_line = hla_typing._risk_calls(near)[1]
    check("a near-miss DQ2 allele is not an exclusion",
          "not determined" in near_line and "neither" not in near_line)
    missed = hla_typing._risk_calls(near, {"DQ2.5": True})[1]
    check("a tag-positive DQ2.5 that typing misses is unresolved",
          "unresolved" in missed and "neither" not in missed)
    check("a tag-positive DQ2.5 with unrelated alleles is unresolved",
          "unresolved" in hla_typing._risk_calls(other, {"DQ2.5": True})[1])
    check("a tag-positive DQ2.5 without DQ typing is unresolved",
          "unresolved" in hla_typing._risk_calls({}, {"DQ2.5": True})[1])
    check("a typed DQ2.5 with a negative tag stands (trans)",
          ": DQ2.5 " in hla_typing._risk_calls(trans, {"DQ2.5": False})[1])
    check("a typed DQ8 with a negative tag is unresolved",
          "unresolved" in hla_typing._risk_calls(dq8, {"DQ8": False})[1])
    check("agreeing negatives stay neither",
          "neither" in hla_typing._risk_calls(other, {"DQ2.5": False, "DQ8": False})[1])
    check("a tag no-call is not a disagreement",
          "neither" in hla_typing._risk_calls(other, {"DQ2.5": None, "DQ8": None})[1])
    check("B27 tag positive with no typed B*27 is unresolved",
          "unresolved" in hla_typing._risk_calls(other, {"B27": True})[0])
    check("typed B*27 with a negative tag is unresolved",
          "unresolved" in hla_typing._risk_calls(b27, {"B27": False})[0])
    check("agreeing B27 negatives stay absent",
          "absent" in hla_typing._risk_calls(other, {"B27": False})[0])


def test_tag_status() -> None:
    print("hla.tag_status:")
    vcf = _vcf("hla_status", [])
    rules = TMP / "hla_status.tsv"
    rules.write_text("trait\trsid\tchrom\tpos\tref\teffect\tinterp0\tinterp1\tinterp2\n"
                     "HLA-B27 (ankylosing spondylitis)\trs0000003\t1\t140\tG\tA\tno\thet\thom\n"
                     "Celiac DQ2.5 haplotype\trs0000001\t1\t150\tG\tT\tno\thet\thom\n"
                     "Celiac DQ8 haplotype\trs0000002\t1\t160\tG\tT\tno\thet\thom\n")
    forced = {("1", "140"): {"nocall": False, "alleles": ("G", "G")},
              ("1", "150"): {"nocall": False, "alleles": ("G", "T")},
              ("1", "160"): {"nocall": True, "alleles": (".", ".")}}
    saved = config.HLA_FILE, hla._called_vcf, forcecall.read_genotypes
    config.HLA_FILE, hla._called_vcf = rules, lambda s: vcf
    hla.forcecall.read_genotypes = lambda s: forced
    try:
        status = hla.tag_status("Maria")
        hla._called_vcf = lambda s: TMP / "no_such.vcf.gz"
        missing = hla.tag_status("Maria")
    finally:
        config.HLA_FILE, hla._called_vcf, forcecall.read_genotypes = saved
    check("reads carrier, non-carrier and no-call",
          status == {"B27": False, "DQ2.5": True, "DQ8": None})
    check("an unreadable VCF means no cross-check", missing == {})


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


def test_t1k_parsing() -> None:
    print("hla_typing._parse_t1k (quality, ties and rejected alleles):")
    tsv = TMP / "t1k_genotype.tsv"
    tsv.write_text("\n".join([
        "HLA-A\t2\tHLA-A*01:01:01\t50.0\t30\tHLA-A*02:01:01\t48.0\t28\t.",
        "HLA-B\t2\tHLA-B*07:02:01\t40.0\t25\tHLA-B*08:01:01\t12.0\t5\t.",
        "HLA-C\t2\tHLA-C*04:01:01,HLA-C*04:09\t20.0\t30\tHLA-C*07:01:01\t20.0\t30\t.",
        "HLA-DPA1\t1\tHLA-DPA1*01:03:01\t90.0\t60\t.\t0\t-1\t.",
        "HLA-DQA1\t2\tHLA-DQA1*01:01:01\t50.0\t30\tHLA-DQA1*05:01:01\t2.0\t0\t.",
        "HLA-DRB3\t1\tHLA-DRB3*01:01:02\t30.0\t30\t.\t0\t-1\t.",
        "HLA-DQB1\t0\t.\t0\t-1\t.\t0\t-1\t.",
    ]) + "\n")
    calls = hla_typing._parse_t1k(tsv)
    check("a clean call is not low confidence", calls["A"]["low"] is False)
    check("names are cut to two fields", calls["A"]["alleles"][0]["names"] == ["A*01:01"])
    check("a quality under the threshold is low", calls["B"]["low"] is True)
    check("ties across 2-field alleles are low",
          calls["C"]["low"] is True and calls["C"]["alleles"][0]["names"] == ["C*04:01", "C*04:09"])
    check("a homozygous call keeps one allele",
          len(calls["DPA1"]["alleles"]) == 1 and calls["DPA1"]["low"] is False)
    check("a quality-0 allele is dropped and flags the locus",
          len(calls["DQA1"]["alleles"]) == 1 and calls["DQA1"]["low"] is True)
    check("loci outside the panel are ignored", "DRB3" not in calls)
    check("an empty locus is a no call", "DQB1" not in calls)
    check("the report row shows low confidence",
          hla_typing._locus_row("B", calls["B"]).endswith("| 25, 5 | low |"))
    check("the report row shows a no call",
          "no call" in hla_typing._locus_row("DQB1", calls.get("DQB1")))


def test_low_confidence_risk_calls() -> None:
    print("hla_typing._risk_calls (low-confidence loci):")
    other = {"B": ["B*07:02", "B*08:01"], "DQA1": ["DQA1*01:01"], "DQB1": ["DQB1*05:01"]}
    check("low B without a tag is not a negative",
          "not determined" in hla_typing._risk_calls(other, {}, {"B"})[0])
    check("low B with an agreeing negative tag is absent",
          ": absent (low-confidence" in hla_typing._risk_calls(other, {"B27": False}, {"B"})[0])
    b27 = {"B": ["B*27:05", "B*08:01"]}
    check("low B*27 with an agreeing tag is present",
          ": present (low-confidence" in hla_typing._risk_calls(b27, {"B27": True}, {"B"})[0])
    check("low B*27 without a tag is not determined",
          "not determined" in hla_typing._risk_calls(b27, {}, {"B"})[0])
    check("low B against the tag stays unresolved",
          "unresolved" in hla_typing._risk_calls(other, {"B27": True}, {"B"})[0])
    subtype = {"B": ["B*27:09", "B*08:01"]}
    check("a low non-associated subtype is not called",
          "not determined" in hla_typing._risk_calls(subtype, {"B27": True}, {"B"})[0])
    check("low DQ cannot read as neither",
          "not determined" in hla_typing._risk_calls(other, {"DQ2.5": False}, {"DQB1"})[1])
    dq25 = {"DQA1": ["DQA1*05:01"], "DQB1": ["DQB1*02:01"]}
    check("low DQ2.5 with an agreeing tag stands",
          ": DQ2.5 (low-confidence" in hla_typing._risk_calls(dq25, {"DQ2.5": True}, {"DQA1"})[1])
    check("low DQ2.5 without a tag is not determined",
          "not determined" in hla_typing._risk_calls(dq25, {}, {"DQA1"})[1])
    check("a low locus elsewhere leaves the risk calls alone",
          ": absent (" in hla_typing._risk_calls(other, {"B27": False}, {"C"})[0])


def test_t1k_read_selection() -> None:
    print("hla_typing._mhc_contig / _read_filter:")
    saved = subprocess.run

    def contig(sq):
        out = f"@HD\tVN:1.6\n@SQ\tSN:1\tLN:400\n{sq}\n"
        subprocess.run = lambda *a, **k: type("R", (), {"stdout": out})()
        try:
            return hla_typing._mhc_contig(pathlib.Path("/data/Jan.cram"))
        except SystemExit as e:
            return f"error: {e}"
        finally:
            subprocess.run = saved

    length = config.HLA_CHR6_LENGTH
    check("Ensembl naming", contig(f"@SQ\tSN:6\tLN:{length}") == "6")
    check("UCSC naming", contig(f"@SQ\tSN:chr6\tLN:{length}") == "chr6")
    check("another build stops", contig("@SQ\tSN:6\tLN:171115067").startswith("error"))
    check("a missing chromosome stops", contig("@SQ\tSN:7\tLN:400").startswith("error"))
    expr = hla_typing._read_filter("chr6")
    check("keeps unmapped reads and their mates", "flag.unmap" in expr and "flag.munmap" in expr)
    check("keeps mates of MHC reads", 'rnext=="chr6"' in expr and "pnext>=28510120" in expr)


def test_t1k_index_checksum() -> None:
    print("hla_typing.ensure_t1k_index (pinned reference):")
    saved = config.T1K_DIR, config.T1K_INDEX_DIR, hla_typing.require_tools
    root = TMP / "t1k"
    root.mkdir()
    config.T1K_DIR, config.T1K_INDEX_DIR = root, root / "imgt-test"
    hla_typing.require_tools = lambda *a: None
    try:
        (root / "hla.dat.zip").write_bytes(b"not the pinned file")
        try:
            hla_typing.ensure_t1k_index()
            refused = False
        except SystemExit:
            refused = True
        check("a reference with the wrong checksum is refused", refused)
        config.T1K_INDEX_DIR.mkdir()
        built = config.T1K_INDEX_DIR / "hlaidx_dna_seq.fa"
        built.write_text(">x\nA\n")
        check("an existing index is reused", hla_typing.ensure_t1k_index() == built)
    finally:
        config.T1K_DIR, config.T1K_INDEX_DIR, hla_typing.require_tools = saved


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


TESTS = [test_trans_dq25_is_recognised, test_tag_status, test_tag_negative_does_not_reassure,
         test_t1k_parsing,
         test_low_confidence_risk_calls,
         test_t1k_read_selection,
         test_t1k_index_checksum,
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
