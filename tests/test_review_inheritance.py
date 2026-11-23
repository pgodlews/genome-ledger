#!/usr/bin/env python3
"""Regression tests for phasing, de-novo, segregation and family-inheritance verdicts.

    uv run python tests/test_review_inheritance.py

Runs against a throwaway GENOMES_ROOT and invented data. bcftools, samtools and tabix do the
real normalizing and piling up, because several of these bugs lived in what those tools
actually write rather than in what the code assumed they write. Exits non-zero on failure.
"""

from __future__ import annotations

import contextlib
import os
import pathlib
import random
import re
import subprocess
import sys
import tempfile
import traceback

_ROOT = pathlib.Path(tempfile.mkdtemp(prefix="genome_review_test_"))
os.environ["GENOMES_ROOT"] = str(_ROOT)      # before pipeline.config reads it
REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import duckdb  # noqa: E402

from pipeline import alphagenome_client as ag  # noqa: E402
from pipeline import config, denovo as dn, enformer, inheritance, pedigree, phasing  # noqa: E402
from pipeline import secondhit  # noqa: E402
from pipeline import segregation as sg  # noqa: E402
from pipeline.load import connect  # noqa: E402
from pipeline.normalize import normalized_path  # noqa: E402

_failures: list[str] = []


def check(name: str, cond: bool) -> None:
    print(f"  {'PASS' if cond else 'FAIL'}  {name}")
    if not cond:
        _failures.append(name)


@contextlib.contextmanager
def patched(obj, **attrs):
    saved = {k: getattr(obj, k) for k in attrs}
    for k, v in attrs.items():
        setattr(obj, k, v)
    try:
        yield
    finally:
        for k, v in saved.items():
            setattr(obj, k, v)


# --- fixtures ----------------------------------------------------------------
_rng = random.Random(11)
_seq = [_rng.choice("ACGT") for _ in range(3000)]
for _pos, _bases in ((998, "CAATGCA"), (1498, "CGATTCG"), (1998, "TAGCATC"), (2198, "TCAGTCA")):
    _seq[_pos - 1:_pos - 1 + len(_bases)] = list(_bases)
REF = "".join(_seq)
WORK = _ROOT / "fixtures"
WORK.mkdir()
FASTA = WORK / "ref.fa"
FASTA.write_text(">1\n" + "\n".join(REF[i:i + 60] for i in range(0, len(REF), 60)) + "\n")
subprocess.run(["samtools", "faidx", str(FASTA)], check=True)
config.REF_FASTA = FASTA
config.LOGS_DIR.mkdir(parents=True, exist_ok=True)

VCF_HEADER = ("##fileformat=VCFv4.2\n##contig=<ID=1,length=3000>\n"
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
              '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="p">\n'
              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="d">\n'
              '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="a">\n'
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")


def write_vcf(path: pathlib.Path, records, fmt="GT:DP:AD") -> pathlib.Path:
    """records: (pos, ref, alt, sample-column). Bgzipped and indexed when the name says so."""
    body = "".join(f"1\t{p}\t.\t{r}\t{a}\t50\tPASS\t.\t{fmt}\t{s}\n" for p, r, a, s in records)
    if path.suffix != ".gz":
        path.write_text(VCF_HEADER + body)
        return path
    raw = path.with_suffix("")
    raw.write_text(VCF_HEADER + body)
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(path), str(raw)], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(path)], check=True)
    return path


def read(name: str, start: int, sub: dict | None = None, ins: tuple | None = None) -> str:
    """One 150 bp SAM record, with base substitutions and at most one insertion after a POS."""
    bases = [(sub or {}).get(p, REF[p - 1]) for p in range(start, start + 150)]
    if ins:
        cut = ins[0] - start + 1
        left, right = bases[:cut], bases[cut:len(bases) - len(ins[1])]
        seq, cigar = "".join(left) + ins[1] + "".join(right), f"{cut}M{len(ins[1])}I{len(right)}M"
    else:
        seq, cigar = "".join(bases), "150M"
    return f"{name}\t0\t1\t{start}\t60\t{cigar}\t*\t0\t0\t{seq}\t{'I' * len(seq)}\tRG:Z:r"


def bam(name: str, reads: list[str]) -> pathlib.Path:
    sam = "@HD\tVN:1.6\n@SQ\tSN:1\tLN:3000\n@RG\tID:r\tSM:S\n" + "\n".join(reads) + "\n"
    out = WORK / f"{name}.bam"
    subprocess.run(["samtools", "sort", "-o", str(out), "-"], input=sam.encode(), check=True)
    subprocess.run(["samtools", "index", str(out)], check=True)
    return out


def set_family(samples: list[str]) -> None:
    """Adam + Maria -> Lena, with the given callsets ingested."""
    config.FAMILY_FILE.parent.mkdir(parents=True, exist_ok=True)
    config.SAMPLES_TSV.parent.mkdir(parents=True, exist_ok=True)
    config.FAMILY_FILE.write_text(
        "id\tdisplay_name\tsex\tfather\tmother\tpartner\n"
        "Adam\tAdam\tmale\t\t\t\nMaria\tMaria\tfemale\t\t\t\n"
        "Lena\tLena\tfemale\tAdam\tMaria\t\n")
    config.SAMPLES_TSV.write_text("sample_id\n" + "".join(f"{s}\n" for s in samples))


def add_row(con, snap, sid, chrom, pos, ref, alt, sig="Pathogenic", gene="GENE1"):
    con.execute("INSERT INTO variants (snapshot_id, sample_id, chrom, pos, ref, alt, gene, "
                "consequence, clinvar_sig, gnomad_af, zygosity) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, 'intron_variant', ?, 0.001, 'HET')",
                [snap, sid, chrom, pos, ref, alt, gene, sig])


def verdict(phased, *variants) -> dict:
    cg = {"gene": "GENE1", "chrom": "1",
          "variants": [{"pos": p, "ref": r, "alt": a, "csq": "frameshift", "sig": None}
                       for p, r, a in variants]}
    return phasing._verdict(cg, phasing._parse_phased(phased))


# --- phasing -----------------------------------------------------------------
def test_split_multiallelic_is_never_in_cis() -> None:
    print("phasing: the two alleles of a split 1/2 indel")
    # The shape WhatsHap writes: first record phased, duplicate POS left unphased.
    whatshap = write_vcf(WORK / "split.vcf", [(1000, "AT", "A", "1|0:1000"),
                                              (1000, "ATG", "A", "0/1:."),
                                              (1040, "C", "T", "0|1:1000")], fmt="GT:PS")
    v = verdict(whatshap, (1001, "T", "-"), (1001, "TG", "-"))
    check("WhatsHap's output is not called in-cis", "in-cis" not in v["verdict"])
    both = write_vcf(WORK / "split2.vcf", [(1000, "AT", "A", "1|0:1000"),
                                           (1000, "ATG", "A", "0|1:1000")], fmt="GT:PS")
    v = verdict(both, (1001, "T", "-"), (1001, "TG", "-"))
    check("two phased records at one POS are not called in-cis", "in-cis" not in v["verdict"])


def test_indel_beside_snv_keeps_its_own_phase() -> None:
    print("phasing: an indel and a SNV sharing a VEP coordinate")
    dele = write_vcf(WORK / "adj1.vcf", [(1000, "AT", "A", "0|1:1000"),
                                         (1001, "T", "C", "1|0:1000")], fmt="GT:PS")
    v = verdict(dele, (1001, "T", "-"), (1001, "T", "C"))
    check("deletion and SNV at the deleted base read in-trans", v["in_trans"])
    ins = write_vcf(WORK / "adj2.vcf", [(1000, "A", "AG", "1|0:1000"),
                                        (1001, "T", "C", "0|1:1000")], fmt="GT:PS")
    v = verdict(ins, (1001, "-", "G"), (1001, "T", "C"))
    check("insertion and SNV at the next base read in-trans", v["in_trans"])
    cis = write_vcf(WORK / "adj3.vcf", [(1000, "A", "AG", "1|0:1000"),
                                        (1001, "T", "C", "1|0:1000")], fmt="GT:PS")
    check("the same pair on one copy is still in-cis",
          verdict(cis, (1001, "-", "G"), (1001, "T", "C"))["verdict"].startswith("in-cis"))


# --- parental evidence -------------------------------------------------------
def test_partial_no_call_is_not_hom_ref() -> None:
    print("denovo.parent_verdict: partial no-calls")
    for gt in ("0/.", "0|."):
        check(f"{gt} is indeterminate", dn.parent_verdict(gt, 30, 0)[0] == dn.INDETERMINATE)
    check("0/0 is still evidenced hom-ref", dn.parent_verdict("0/0", 30, 0)[0] == dn.CONFIRMED)
    check("1/. still carries", dn.parent_verdict("1/.", 30, 0)[0] == dn.INHERITED)


def test_indel_record_does_not_hide_snv_reads() -> None:
    print("denovo: mpileup's SNP and INDEL records at one POS")
    reads = [read(f"r{i}", 2120 + i) for i in range(11)]
    reads += [read(f"s{i}", 2120 + i, sub={2200: "G"}) for i in range(11)]
    reads += [read(f"x{i}", 2125 + i, ins=(2200, "TT")) for i in range(2)]
    parent_bam = bam("adam_snv", reads)
    homref = write_vcf(WORK / "homref.vcf.gz", [(2200, "A", "G", "0/0:24:24,0")])
    wd = dn._work_dir("Lena")
    wd.mkdir(parents=True, exist_ok=True)
    cand = {"chrom": "1", "pos": 2200, "ref": "A", "alt": "G", "dp": 40, "gq": 50, "vaf": 0.5,
            "kind": "SNV"}
    with patched(dn, aligned_cram=lambda p: parent_bam,
                 _force_call_parent=lambda *a, **k: homref):
        res = dn.confirm("Lena", [cand], ["Adam"], force=True)
    check("11 alt reads are not hidden by 2 indel reads", res[0]["status"] == dn.INHERITED)


def _evidence(reads: list[str], force_called: list, child_alt: str, name: str):
    """Run segregation._read_evidence for one parent with a real pileup."""
    parent_bam = bam(name, reads)
    v = {"chrom": "1", "pos": 2200, "ref": "A", "alt": child_alt, "vcf_pos": 2200,
         "vcf_ref": "A", "vcf_alt": child_alt, "key": ("1", 2200, "A", child_alt)}

    def fake_call(parent, cram, vcf_sites, intervals, out, key, tag="x"):
        write_vcf(pathlib.Path(out), force_called)
        return out

    with patched(sg, aligned_cram=lambda p: parent_bam, call_at_sites=fake_call), \
         patched(dn, aligned_cram=lambda p: parent_bam):
        ev = sg._read_evidence(f"Child_{name}", ["Adam"], [v], {"Adam": set()}, force=True)
    return ev[("Adam", v["key"])][0]


def test_parent_evidence_is_per_allele() -> None:
    print("segregation: a different ALT at the position is not the child's allele")
    other = [read(f"r{i}", 2120 + i) for i in range(12)]
    other += [read(f"c{i}", 2120 + i, sub={2200: "C"}) for i in range(12)]
    check("parent with a different ALT is unknown, not a carrier",
          _evidence(other, [(2200, "A", "G,C", "0/2:24:12,0,12")], "G", "adam_other") == "unknown")
    same = [read(f"r{i}", 2120 + i) for i in range(12)]
    same += [read(f"g{i}", 2120 + i, sub={2200: "G"}) for i in range(12)]
    check("parent with the child's own ALT still carries it",
          _evidence(same, [(2200, "A", "G", "0/1:24:12,12")], "G", "adam_same") == "carries")
    clean = [read(f"r{i}", 2120 + i) for i in range(24)]
    check("evidenced hom-ref is still absent",
          _evidence(clean, [(2200, "A", "G", "0/0:24:24,0")], "G", "adam_clean") == "absent")


def test_joint_record_alleles_resolve() -> None:
    print("denovo: a joint force-call record resolves to the candidate's allele")
    joint = write_vcf(WORK / "joint.vcf.gz", [(1000, "ATG", "AG,A", "0/1:30:15,15,0")])
    wd = dn._work_dir("Lena")
    empty = wd / "empty_regions.txt"
    empty.write_text("1\t1000\n")
    no_reads = bam("adam_joint", [read(f"r{i}", 930 + i) for i in range(30)])
    for alt_ref, label, want in (("AT", "the allele they carry", dn.INHERITED),
                                 ("ATG", "the allele they do not", dn.INDETERMINATE)):
        cand = {"chrom": "1", "pos": 1000, "ref": alt_ref, "alt": "A", "dp": 40, "gq": 50,
                "vaf": 0.5, "kind": "indel"}
        with patched(dn, aligned_cram=lambda p: no_reads,
                     _force_call_parent=lambda *a, **k: joint):
            res = dn.confirm("Lena", [cand], ["Adam"], force=True)
        check(f"{alt_ref}>A, {label}: {want}", res[0]["status"] == want)


# --- segregation -------------------------------------------------------------
def test_each_parent_means_in_trans() -> None:
    print("segregation.resolve: a variant from each parent")
    for third in (sg.UNRESOLVED, sg.DE_NOVO, sg.BOTH):
        tr = [{"transmission": t} for t in (sg.PATERNAL, sg.MATERNAL, third)]
        check(f"paternal + maternal + {third!r} is in-trans", sg.resolve(tr)[0] == sg.IN_TRANS)
    check("paternal + unresolved stays unresolved",
          sg.resolve([{"transmission": sg.PATERNAL},
                      {"transmission": sg.UNRESOLVED}])[0] == sg.UNRESOLVED)


def test_delins_matches_normalized_vcf() -> None:
    print("segregation._vcf_alleles: delins as bcftools norm writes them")
    raw = write_vcf(WORK / "lena.raw.vcf", [(1000, "AT", "A", "0/1:40:20,20"),
                                            (1500, "ATT", "AG", "0/1:40:20,20"),
                                            (2000, "GCA", "T", "0/1:40:20,20")])
    dest = normalized_path("Lena")
    dest.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["bcftools", "norm", "-f", str(FASTA), "-m", "-both", "-Oz", "-o", str(dest),
                    str(raw)], check=True, capture_output=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(dest)], check=True)
    cands = [{"chrom": "1", "pos": 1001, "ref": "T", "alt": "-"},
             {"chrom": "1", "pos": 1501, "ref": "TT", "alt": "G"},
             {"chrom": "1", "pos": 2000, "ref": "GCA", "alt": "T"},
             {"chrom": "1", "pos": 1502, "ref": "TT", "alt": "G"}]   # right alleles, wrong POS
    sg._vcf_alleles("Lena", cands)
    check("deletion resolves", cands[0]["vcf_pos"] == 1000)
    check("trimmed delins resolves", cands[1]["vcf_pos"] == 1501)
    check("unanchored delins resolves", cands[2]["vcf_pos"] == 2000)
    check("same alleles one base away do not", cands[3]["vcf_pos"] is None)
    check("a different event sharing the tail does not match",
          not sg._same_event({"ref": "TT", "alt": "G"}, "CTT", "GG"))


# --- de-novo coordinates -----------------------------------------------------
def test_denovo_uses_vep_form() -> None:
    print("denovo: score lookups and confirmed.tsv use VEP's allele form")
    config.ENFORMER_SCORES_FILE.parent.mkdir(parents=True, exist_ok=True)
    config.ENFORMER_SCORES_FILE.write_text(
        "chrom\tpos\tref\talt\tdelta_max\tl2_center\ttop_track\ttop_desc\n"
        "3\t101\tT\t-\t7.5\t1.0\t5\tCAGE:x\n3\t201\t-\tT\t7.5\t1.0\t5\tCAGE:x\n")
    enformer._table = None
    try:
        for label, pos, ref, alt in (("deletion", 100, "AT", "A"), ("insertion", 200, "A", "AT")):
            row = dn._row("Lena", [], {"chrom": "3", "pos": pos, "ref": ref, "alt": alt,
                                       "kind": "indel", "dp": 30, "vaf": 0.5, "gene": "GENE1",
                                       "csq": "intron_variant"})
            check(f"{label} finds its regulatory score", "Δ7.50" in row)
            check(f"{label} is still displayed in VCF form", f"{ref}>{alt}" in row)
    finally:
        enformer._table = None

    set_family(["Lena", "Adam", "Maria"])
    cands = [{"chrom": "3", "pos": 100, "ref": "AT", "alt": "A", "status": dn.CONFIRMED},
             {"chrom": "3", "pos": 300, "ref": "C", "alt": "G", "status": dn.CONFIRMED}]
    with patched(dn, candidates=lambda *a, **k: (cands, {}), confirm=lambda c, cs, *a, **k: cs,
                 _annotate=lambda *a: None, _report=lambda *a, **k: WORK / "report.md"):
        dn.denovo("Lena", force=True)
    lines = (dn._work_dir("Lena") / "confirmed.tsv").read_text().split()
    check("confirmed.tsv holds the deletion as DuckDB does", "3\t101\tT\t-" in
          (dn._work_dir("Lena") / "confirmed.tsv").read_text())
    check("confirmed.tsv keeps the SNV", "300" in lines)


# --- inheritance page --------------------------------------------------------
def test_old_snapshot_ignores_later_parent() -> None:
    print("inheritance: a parent ingested after the snapshot was built")
    set_family(["Lena", "Maria", "Adam"])
    snap = "2031-01-01"
    con = connect()
    add_row(con, snap, "Lena", "1", 100, "C", "T")
    add_row(con, snap, "Maria", "1", 100, "C", "T")
    add_row(con, snap, "Lena", "1", 200, "G", "A")
    con.close()
    html = inheritance.build_page(snap, out_dir=_ROOT / "html_old").read_text()
    origins = re.findall(r'"origin": \{"Lena": "(\w+)"', html)
    check("neither variant gets a parental origin or de-novo label",
          origins == [inheritance.UNRESOLVED, inheritance.UNRESOLVED])
    unknown = re.search(r'"unknown": \[([^\]]*)\]', html)
    check("the later parent shows as not sequenced", bool(unknown) and "Adam" in unknown.group(1))


def test_mt_is_not_partial_support() -> None:
    print("inheritance: MT reads only the primary callset")
    set_family(["Lena", "Lena_vendor", "Maria", "Adam"])
    snap = "2031-02-01"
    con = connect()
    for sid in ("Lena", "Lena_vendor", "Maria", "Adam"):
        add_row(con, snap, sid, "1", 5000, "A", "G", sig="Benign")
    for sid in ("Lena", "Lena_vendor", "Maria"):
        add_row(con, snap, sid, "MT", 7000, "A", "G", gene="GENE2")
    add_row(con, snap, "Lena", "1", 300, "C", "T")
    con.close()
    fam = pedigree.load(sequenced=pedigree.sequenced_ids())
    by_chrom = {v["key"][0]: v for v in inheritance.collect(snap, fam)}
    check("an MT call is not discordant", inheritance.discordant(by_chrom["MT"], fam) == [])
    check("...and has no n/m marker",
          inheritance.callset_support(by_chrom["MT"]["carriers"]["Lena"], "Lena", fam) == "")
    check("a real single-callset autosomal call is still flagged",
          inheritance.discordant(by_chrom["1"], fam) == ["Lena"])


# --- second hit --------------------------------------------------------------
_A = {"sample": "Lena", "gene": "GENE1", "chrom": "2", "pos": 1000, "ref": "C", "alt": "T",
      "consequence": "missense_variant", "sig": "Pathogenic", "af": 0.001}
_B = {**_A, "pos": 900000, "ref": "G", "alt": "A", "sig": None}


def test_unsequenced_parent_gives_unknown_phase() -> None:
    print("second-hit: parental origin needs both parents")
    set_family(["Lena", "Adam"])
    snap = "2031-03-01"
    con = connect()
    for sid in ("Lena", "Adam"):
        for x in (_A, _B):
            add_row(con, snap, sid, "2", x["pos"], x["ref"], x["alt"], sig=x["sig"])
    con.close()
    check("one sequenced parent carrying both is unknown, not cis",
          secondhit.parental_phase(snap, {**_B, "known": _A}) == "unknown")
    set_family(["Lena", "Adam", "Maria"])
    con = connect()
    add_row(con, snap, "Maria", "1", 1, "A", "C", sig=None)
    con.close()
    check("both parents present, one carrying both, is still cis",
          secondhit.parental_phase(snap, {**_B, "known": _A}).startswith("cis"))


def test_unknown_read_phase_keeps_origin() -> None:
    print("second-hit: --phase with no read evidence")
    set_family(["Lena", "Adam", "Maria"])
    snap = "2031-04-01"
    con = connect()
    for sid, x in (("Lena", _A), ("Lena", _B), ("Adam", _A), ("Maria", _B)):
        add_row(con, snap, sid, "2", x["pos"], x["ref"], x["alt"], sig=x["sig"])
    con.close()

    def phase_cell(read_backed: str) -> str:
        with patched(config, INTERPRET_MODEL="alphagenome2"), \
             patched(secondhit, select=lambda s: [{**_B, "known": _A}]), \
             patched(ag, score_selected=lambda *a, **k: ({}, 0), evidence=lambda *a, **k: [],
                     summary=lambda *a, **k: [], splice_label=lambda *a, **k: "-",
                     conditions_for=lambda g: [{"label": "t", "condition": "c"}]), \
             patched(phasing, phase_pair=lambda *a, **k: read_backed):
            rep = secondhit.second_hit(snap, phase=True)
        row = next(l for l in rep.read_text().splitlines() if l.startswith("| Lena"))
        return row.split("|")[5].strip()

    check("'unknown' does not erase a parental-origin call",
          phase_cell("unknown") == "trans (parental-origin candidate)")
    check("a real read-backed call still replaces it",
          phase_cell("cis (read-backed)") == "cis (read-backed)")


# --- PED ---------------------------------------------------------------------
def test_ped_has_no_dangling_parent() -> None:
    print("pedigree: a parent ingested only under a suffix")
    p = pedigree.Person
    fam = pedigree.Family(
        {"Jan": p("Jan", "Jan", "male", None, None, None),
         "Maria": p("Maria", "Maria", "female", None, None, None),
         "Lena": p("Lena", "Lena", "female", "Jan", "Maria", None)},
        sequenced={"Jan_vendor", "Maria", "Maria_t2t", "Lena", "Lena_t2t"})
    out = pedigree.generate_pedigree_ped(fam, WORK / "family.ped")
    rows = [l.split("\t") for l in out.read_text().splitlines()]
    ids = {r[1] for r in rows}
    check("every parent column names a PED row",
          all(par == "0" or par in ids for r in rows for par in r[2:4]))
    by_id = {r[1]: r for r in rows}
    check("the suffixed parent is the one used", by_id["Lena"][2] == "Jan_vendor")
    check("a same-suffix parent is still preferred", by_id["Lena_t2t"][3] == "Maria_t2t")


# --- query -------------------------------------------------------------------
def test_segregation_query() -> None:
    print("queries/segregation.sql")
    sql = (REPO / "queries" / "segregation.sql").read_text().replace("${snapshot}", "S")
    con = duckdb.connect()
    con.execute("CREATE TABLE variants (snapshot_id VARCHAR, sample_id VARCHAR, gene VARCHAR, "
                "chrom VARCHAR, pos BIGINT, ref VARCHAR, alt VARCHAR, clinvar_sig VARCHAR, "
                "zygosity VARCHAR)")
    wanted = ["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic",
              "Pathogenic,_low_penetrance", "Pathogenic|risk_factor",
              "Pathogenic/Pathogenic,_low_penetrance|risk_factor"]
    unwanted = ["Conflicting_classifications_of_pathogenicity", "Benign"]
    for i, s in enumerate(wanted + unwanted):
        con.execute("INSERT INTO variants VALUES ('S','Lena','G','1',?,'A','C',?,'HET')",
                    [100 + i, s])
    got = {r[5] for r in con.execute(sql).fetchall()}
    check("every (likely) pathogenic term is kept", set(wanted) <= got)
    check("non-pathogenic terms stay out", not got & set(unwanted))
    con.execute("DELETE FROM variants")
    con.execute("INSERT INTO variants VALUES ('S','Lena','G','1',500,'A','C','Pathogenic','HET')")
    con.execute("INSERT INTO variants VALUES ('S','Adam','G','1',500,'A','C','Pathogenic',NULL)")
    row = con.execute(sql).fetchone()
    check("a NULL-zygosity carrier stays in the carrier list",
          "Adam" in row[6] and "Lena" in row[6] and row[7] == 2)


if __name__ == "__main__":
    tests = [v for k, v in list(globals().items()) if k.startswith("test_") and callable(v)]
    for test in tests:
        try:
            test()
        except Exception:  # noqa: BLE001 — a crash is a failure of that test, not of the run
            traceback.print_exc()
            _failures.append(f"{test.__name__} raised")
    print()
    if _failures:
        print(f"{len(_failures)} test(s) FAILED: {_failures}")
        sys.exit(1)
    print("ALL TESTS PASSED")
