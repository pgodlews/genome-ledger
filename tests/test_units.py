#!/usr/bin/env python3
"""Dependency-free unit tests for the new pipeline logic (no pytest needed).

    uv run python tests/test_units.py

Covers the inheritance-mode interpretation (panel._status) and the sub-chromosomal
scatter (call._scatter_tasks) — the two pieces with non-trivial branching that the
end-to-end run doesn't exercise in isolation.
"""

from __future__ import annotations

import gzip
import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import call, config, panel  # noqa: E402

# Reports are filed under config.REPORTS_DIR whichever stage writes them, so a test that
# points one stage directory at a tempdir would still write its report into the real
# $GENOMES_ROOT/reports. Point the reports tree at a throwaway directory for the whole run;
# tests that care about its contents set their own and restore this one.
import tempfile as _tempfile  # noqa: E402
config.REPORTS_DIR = pathlib.Path(_tempfile.mkdtemp(prefix="genome_test_reports_")) / "reports"
config.REPORT_DIRS = (config.REPORTS_DIR,)

_failures: list[str] = []


def render_formats_count() -> dict:
    """The registered formats — a revoked report must go from every one of their trees."""
    from pipeline import render
    return render.FORMATS


def check(name: str, cond: bool) -> None:
    print(f"  {'PASS' if cond else 'FAIL'}  {name}")
    if not cond:
        _failures.append(name)


def _v(gene: str, key, zyg: str) -> dict:
    return {"gene": gene, "key": key, "zyg": zyg, "sig": "Pathogenic",
            "rev": "", "disease": "x", "label": "x"}


def test_panel_status() -> None:
    print("panel._status (inheritance-aware carrier/affected):")
    s = panel._status
    check("AD het            -> at-risk",
          s("AD", "female", [_v("BRCA1", ("17", 1, "A", "T"), "HET")]) == "at-risk")
    check("AR het            -> carrier",
          s("AR", "female", [_v("CFTR", ("7", 1, "A", "T"), "HET")]) == "carrier")
    check("AR hom            -> affected",
          s("AR", "female", [_v("CFTR", ("7", 1, "A", "T"), "HOM")]) == "affected")
    check("AR 2 diff het     -> compound",
          s("AR", "female", [_v("CFTR", ("7", 1, "A", "T"), "HET"),
                             _v("CFTR", ("7", 2, "C", "G"), "HET")]) == "compound")
    check("AR same het twice -> carrier (not compound)",
          s("AR", "female", [_v("CFTR", ("7", 1, "A", "T"), "HET"),
                             _v("CFTR", ("7", 1, "A", "T"), "HET")]) == "carrier")
    check("XL male het       -> affected (hemizygous)",
          s("XL", "male", [_v("F8", ("X", 1, "A", "T"), "HET")]) == "affected")
    check("XL female het     -> carrier",
          s("XL", "female", [_v("F8", ("X", 1, "A", "T"), "HET")]) == "carrier")
    check("XL female hom     -> affected",
          s("XL", "female", [_v("F8", ("X", 1, "A", "T"), "HOM")]) == "affected")
    check("no variants       -> None",
          s("AD", "female", []) is None)


def test_family_risk_rows() -> None:
    print("panel.family_risk_rows / _pairing_risk_text (AR carrier-index -> family report):")
    dz = panel._Disease("Cystic fibrosis", ["CFTR"], "AR", "")

    # Item 5: a hom/compound (affected) member paired with a carrier must still produce a
    # family-report row — and at the correct (higher) risk, not the carrier x carrier 25%.
    rows = panel.family_risk_rows(
        [dz], {"CFTR": {"SAMPLE002": panel._CARRIER, "SAMPLE001": panel._AFFECTED}})
    check("gene appears in family report", any(r[1] == "CFTR" for r in rows))
    check("affected x carrier -> 50% risk wording",
          [r for r in rows if r[1] == "CFTR"][0][3] == "50%")

    check("carrier x carrier -> 25%",
          panel._pairing_risk_text({"a": panel._CARRIER, "b": panel._CARRIER}) == "25%")
    check("affected x affected -> 100%",
          panel._pairing_risk_text({"a": panel._AFFECTED, "b": panel._AFFECTED}) == "100%")

    # A "possible compound het" is two different damaging alleles with the phase UNKNOWN.
    # In-trans it transmits like an affected homozygote (1.0), in-cis like a carrier (0.5),
    # so the per-child risk is bounded, not determined. Scoring it as a flat 1.0 printed a
    # definitive 50%/100% for a question nobody had answered yet.
    compound_x_carrier = panel._pairing_risk_text(
        {"a": panel._COMPOUND, "b": panel._CARRIER})
    check("compound x carrier -> a range, not a definitive 50%",
          compound_x_carrier.startswith("25-50%".replace("-", "\u2013")))
    check("...and says what would resolve it", "phase" in compound_x_carrier)
    compound_x_compound = panel._pairing_risk_text(
        {"a": panel._COMPOUND, "b": panel._COMPOUND})
    check("compound x compound -> 25-100%, not a definitive 100%",
          compound_x_compound.startswith("25\u2013100%"))
    check("compound x affected -> 50-100%",
          panel._pairing_risk_text(
              {"a": panel._COMPOUND, "b": panel._AFFECTED}).startswith("50\u2013100%"))
    # resolved statuses must still collapse to a single figure
    check("carrier x affected stays a flat 50% (nothing unresolved)",
          panel._pairing_risk_text({"a": panel._CARRIER, "b": panel._AFFECTED}) == "50%")
    # both widening reasons at once are both named
    three = panel._pairing_risk_text(
        {"a": panel._COMPOUND, "b": panel._CARRIER, "c": panel._AFFECTED})
    check("several members AND an unresolved phase name both reasons",
          "which two are partners" in three and "phase unresolved" in three)
    check("the compound-het label says 'possible'",
          panel._SHORT_LABEL[panel._COMPOUND] == "possible compound het")

    # A single AR hit with no second family member shares nothing to report on.
    check("solo carrier (no shared member) -> no row",
          panel.family_risk_rows([dz], {"CFTR": {"Maria": panel._CARRIER}}) == [])

    # Item 7: gene lists interpolated into SQL must handle the empty set and embedded quotes.
    check("sql_in_list([]) -> NULL (valid SQL, matches nothing; not the syntax-error '()')",
          panel.sql_in_list([]) == "NULL")
    check("sql_in_list escapes embedded quotes",
          panel.sql_in_list(["O'Brien"]) == "'O''Brien'")


def test_ar_status_index_severity() -> None:
    """A gene can sit in two AR panel diseases, so one sample can be affected for one and a
    carrier for the other. The index must keep the MOST SEVERE status — last-wins would let
    a later 'carrier' overwrite an earlier 'affected' and understate the child risk."""
    print("panel.index_ar_status (most-severe wins across AR diseases sharing a gene):")
    from pipeline.panel import _AFFECTED, _CARRIER, _COMPOUND, index_ar_status

    idx: dict[str, dict[str, str]] = {}
    index_ar_status(idx, "G", "kid", _AFFECTED)
    index_ar_status(idx, "G", "kid", _CARRIER)      # later, weaker — must not overwrite
    check("affected survives a later carrier", idx["G"]["kid"] == _AFFECTED)

    idx2: dict[str, dict[str, str]] = {}
    index_ar_status(idx2, "G", "kid", _CARRIER)
    index_ar_status(idx2, "G", "kid", _COMPOUND)    # later, stronger — must upgrade
    check("carrier upgrades to compound het", idx2["G"]["kid"] == _COMPOUND)
    check("distinct samples stay independent",
          (index_ar_status(idx2, "G", "other", _CARRIER),
           idx2["G"] == {"kid": _COMPOUND, "other": _CARRIER})[1])

    # The risk text must follow the severity that survived: affected x carrier is 50%,
    # not the 25% the report used to hard-code.
    from pipeline.panel import _pairing_risk_text
    check("surviving severity drives the risk text",
          _pairing_risk_text({"a": _AFFECTED, "b": _CARRIER}) == "50%")


def test_read_platform_detection() -> None:
    """The @RG PL tag came from a READ_PLATFORM default that named one sequencer, so reads from
    any other were mislabelled unless the user knew to override it. It is now read off the
    first read name, with READ_PLATFORM as an override and ILLUMINA (warned) as the fallback."""
    print("call.read_platform (from the FASTQ read names):")
    import gzip
    import tempfile
    from pipeline import call

    saved = config.READ_PLATFORM
    tmp = pathlib.Path(tempfile.mkdtemp())

    def fq(name, header, gz=True):
        path = tmp / name
        body = f"{header}\nACGT\n+\nIIII\n"
        if gz:
            with gzip.open(path, "wt") as fh:
                fh.write(body)
        else:
            path.write_text(body)
        return path
    try:
        config.READ_PLATFORM = None
        cases = {
            "MGI/DNBSEQ (G400)": ("@V300012345L1C001R0010000001/1", "DNBSEQ"),
            "MGI/DNBSEQ (T7)": ("@E100012345L2C012R03400000123/2", "DNBSEQ"),
            "Illumina CASAVA 1.8+": ("@A00123:8:HABCDXX:1:1101:1000:2000 1:N:0:ATCACG",
                                     "ILLUMINA"),
            "Illumina, no comment": ("@NB501234:12:HXXXXBGXY:1:11101:10000:1000", "ILLUMINA"),
            "Illumina pre-1.8": ("@HWUSI-EAS100R:6:73:941:1973#0/1", "ILLUMINA"),
        }
        for label, (header, want) in cases.items():
            check(f"{label} -> {want}", call.read_platform(fq("r1.fq.gz", header)) == want)
        check("an uncompressed FASTQ is read too",
              call.read_platform(fq("r1.fq", cases["MGI/DNBSEQ (G400)"][0], gz=False))
              == "DNBSEQ")
        check("an unrecognised read name (SRA re-export) falls back to ILLUMINA",
              call.read_platform(fq("r1.fq.gz", "@SRR1234567.1 1 length=150")) == "ILLUMINA")
        check("an unreadable file falls back too", call.read_platform(tmp / "missing.fq.gz")
              == "ILLUMINA")
        config.READ_PLATFORM = "ELEMENT"
        check("READ_PLATFORM overrides detection",
              call.read_platform(fq("r1.fq.gz", cases["MGI/DNBSEQ (G400)"][0])) == "ELEMENT")
    finally:
        config.READ_PLATFORM = saved


def test_fastq_naming_and_lanes() -> None:
    """`align` accepted one provider's naming only (`<id>_1.fq.gz`), so Illumina-named or
    lane-split reads stopped it with "FASTQ not found"."""
    print("call._fastq_files / _lane_stream (FASTQ naming and lanes):")
    import gzip
    import shutil
    import tempfile
    from pipeline import call

    saved = config.INCOMING_DIR
    tmp = pathlib.Path(tempfile.mkdtemp())
    try:
        config.INCOMING_DIR = tmp

        def touch(*names):
            for n in names:
                (tmp / n).write_bytes(b"")

        touch("ACME1_1.fq.gz", "ACME1_2.fq.gz", "ACME12_1.fq.gz", "ACME12_2.fq.gz")
        r1, r2 = call._fastq_files("ACME1")
        check("provider naming: one pair, and ACME1 does not pick up ACME12's files",
              [p.name for p in r1 + r2] == ["ACME1_1.fq.gz", "ACME1_2.fq.gz"])
        touch("KID_R1.fastq.gz", "KID_R2.fastq.gz")
        check("_R1/_R2 naming", [p.name for p in sum(call._fastq_files("KID"), [])]
              == ["KID_R1.fastq.gz", "KID_R2.fastq.gz"])
        touch(*[f"LANE_S1_L00{n}_R{r}_001.fastq.gz" for n in (2, 1) for r in (1, 2)],
              "LANE_S1_L001_I1_001.fastq.gz")
        r1, r2 = call._fastq_files("LANE")
        check("Illumina lane split: lanes in order, index reads ignored",
              [p.name for p in r1] == ["LANE_S1_L001_R1_001.fastq.gz",
                                       "LANE_S1_L002_R1_001.fastq.gz"]
              and [p.name for p in r2] == ["LANE_S1_L001_R2_001.fastq.gz",
                                           "LANE_S1_L002_R2_001.fastq.gz"])
        touch("HALF_L01_1.fq.gz", "HALF_L01_2.fq.gz", "HALF_L02_1.fq.gz")
        try:
            call._fastq_files("HALF")
            refused = False
        except SystemExit as e:
            refused = "HALF_L02_1.fq.gz" in str(e)
        check("a lane without its mate is refused, by name", refused)
        try:
            call._fastq_files("NOBODY")
            missing = False
        except SystemExit as e:
            missing = "_R1_001" in str(e)
        check("no files: the error lists the accepted names", missing)

        # Lanes are streamed, decompressed and in order, into one pipe.
        lanes = []
        for i in (1, 2):
            path = tmp / f"S_L00{i}.fq.gz"
            with gzip.open(path, "wt") as fh:
                fh.write(f"@r{i}\nACGT\n+\nIIII\n")
            lanes.append(path)
        plain = tmp / "S_L003.fq"
        plain.write_text("@r3\nTTTT\n+\nIIII\n")
        fifo = tmp / "pipe.fastq"
        writer = call._lane_stream(lanes + [plain], fifo)
        with open(fifo) as fh:
            got = fh.read()
        check("lane stream: gzip and plain lanes concatenated in order",
              writer.wait() == 0 and [ln for ln in got.splitlines() if ln.startswith("@")]
              == ["@r1", "@r2", "@r3"])
    finally:
        config.INCOMING_DIR = saved
        shutil.rmtree(tmp, ignore_errors=True)


def test_share_rebuilds_headers() -> None:
    """A VCF or CRAM handed to someone else carried the working file's history: command
    lines with home-folder paths, the reference path, run dates, internal and vendor ids."""
    print("share (VCF / CRAM headers rebuilt from an allowlist):")
    import gzip
    import shutil
    import subprocess
    import tempfile
    from pipeline import share

    if not all(shutil.which(t) for t in ("bcftools", "tabix", "samtools", "bgzip")):
        print("  SKIP  (bcftools/tabix/samtools/bgzip not on PATH)")
        return
    names = ("SHARE_DIR", "NORMALIZED_DIR", "CALLED_DIR", "ALIGNED_DIR", "FAMILY_FILE")
    saved = {n: getattr(config, n) for n in names}
    tmp = pathlib.Path(tempfile.mkdtemp())
    try:
        for n, sub in (("SHARE_DIR", "share"), ("NORMALIZED_DIR", "norm"),
                       ("CALLED_DIR", "called"), ("ALIGNED_DIR", "aligned")):
            setattr(config, n, tmp / sub)
            (tmp / sub).mkdir()
        config.FAMILY_FILE = tmp / "family.tsv"
        config.FAMILY_FILE.write_text("id\tdisplay_name\tsex\tfather\tmother\tpartner\n"
                                      "Adam\tAdam\tmale\t\t\t\n")
        private = "/Users/someone/genomes"
        vcf_text = (
            "##fileformat=VCFv4.2\n"
            f'##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller --output {private}/x">\n'
            f"##bcftools_viewCommand=view {private}/raw/x.vcf.gz; Date=Sun Jun 28 14:17:54 2026\n"
            f"##reference=file://{private}/refs/ref.fa\n"
            "##source=HaplotypeCaller\n"
            '##FILTER=<ID=PASS,Description="All filters passed">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "##contig=<ID=1,length=1000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tKIT123456\n"
            "1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n")
        src = config.NORMALIZED_DIR / "Adam_vendor.norm.vcf.gz"
        subprocess.run(["bgzip", "-c"], input=vcf_text.encode(), stdout=open(src, "wb"),
                       check=True)
        out = share.share("Adam_vendor")
        text = gzip.open(out, "rt").read()
        check("VCF: written under the person's name", out.name == f"Adam.{config.BUILD}.vcf.gz"
              and out.parent.name == "Adam")
        check("VCF: no command lines, paths or dates left in the header",
              "/Users/" not in text and "Command" not in text and "Date=" not in text
              and "##source" not in text)
        check("VCF: sample column is the person, not the kit id",
              "\tKIT123456" not in text and text.splitlines()[-2].endswith("\tAdam"))
        check("VCF: records untouched, build named", text.splitlines()[-1]
              == "1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1" and f"##reference={config.BUILD}" in text)

        ref = tmp / "ref.fa"
        ref.write_text(">1\n" + "ACGT" * 250 + "\n")
        subprocess.run(["samtools", "faidx", str(ref)], check=True)
        sam = tmp / "in.sam"
        sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            f"@SQ\tSN:1\tLN:1000\tUR:file://{private}/refs/ref.fa\n"
            "@RG\tID:father\tSM:father\tPL:DNBSEQ\tLB:father\tPU:FLOWCELL1.L1\tDT:2026-06-01\n"
            f"@PG\tID:bwa\tPN:bwa\tCL:bwa mem {private}/refs/ref.fa\n"
            "@CO\tsequenced at home\n"
            "r1\t0\t1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:father\n")
        cram = config.ALIGNED_DIR / "Adam.cram"
        subprocess.run(["samtools", "view", "--no-PG", "-C", "-T", str(ref), "-o", str(cram),
                        str(sam)], check=True)
        out = share.share("Adam", "cram")
        hdr = subprocess.run(["samtools", "view", "-H", "--no-PG", str(out)], check=True,
                             capture_output=True, text=True).stdout
        reads = subprocess.run(["samtools", "view", "--no-PG", "-T", str(ref), str(out)],
                               check=True, capture_output=True, text=True).stdout
        check("CRAM: no @PG/@CO, no paths, no flowcell or date in the read group",
              "@PG" not in hdr and "@CO" not in hdr and "/Users/" not in hdr
              and "PU:" not in hdr and "DT:" not in hdr)
        check("CRAM: read group names the person, keeps the ID the reads point at",
              "\tSM:Adam" in hdr and "\tLB:Adam" in hdr and "ID:father" in hdr)
        check("CRAM: reads intact and indexed", reads.startswith("r1\t0\t1\t1\t60\t4M")
              and pathlib.Path(f"{out}.crai").exists())
    finally:
        for n, v in saved.items():
            setattr(config, n, v)
        shutil.rmtree(tmp, ignore_errors=True)


def test_scatter_tasks() -> None:
    print("call._scatter_tasks (sub-chromosomal windows):")
    # Stub contig lengths so the test needs no reference .fai on disk.
    orig = call._fai_lengths
    call._fai_lengths = lambda: {**{c: 100_000_000 for c in config.CALL_CONTIGS},
                                 config.MT_CONTIG: 16569}
    try:
        tasks = call._scatter_tasks()
    finally:
        call._fai_lengths = orig

    # 100Mb / 30Mb window -> 4 windows per diploid contig, plus one MT shard.
    check("window count", len(tasks) == len(config.CALL_CONTIGS) * 4 + 1)
    check("MT last, whole-contig, haploid",
          tasks[-1] == (config.MT_CONTIG, config.MT_PLOIDY, config.MT_CONTIG))
    check("diploid windows are ploidy 2", all(t[1] == 2 for t in tasks[:-1]))

    c0 = config.CALL_CONTIGS[0]
    wins = [t[0] for t in tasks if t[0].startswith(c0 + ":")]
    prev_end, contiguous = None, True
    for region in wins:
        start, end = map(int, region.split(":")[1].split("-"))
        if prev_end is not None and start != prev_end + 1:
            contiguous = False
        prev_end = end
    check("windows contiguous & non-overlapping within a contig", contiguous)
    labels = [t[2] for t in tasks if t[0].startswith(c0 + ":")]
    check("labels sort in genomic order (zero-padded)", labels == sorted(labels))


def test_diff_classify() -> None:
    """The DELTA triage engine: which annotation changes are real, and how the multi-factor
    score sorts them into the four tiers (actionable_now … ignore)."""
    print("triage.change_kind / score_delta (transitions + tiers):")
    from pipeline import triage
    REV0 = "criteria_provided,_single_submitter"
    REVX = "reviewed_by_expert_panel"

    def kind(*a):
        r = triage.change_kind(*a)
        return r[0] if r else None

    # change detection — (prev_sig, curr_sig, prev_rev, curr_rev, prev_af, curr_af)
    check("VUS -> Pathogenic = new_pathogenic",
          kind("Uncertain_significance", "Pathogenic", REV0, REV0, None, None) == "new_pathogenic")
    check("absent -> Pathogenic = new_pathogenic",
          kind(None, "Pathogenic", None, REV0, None, None) == "new_pathogenic")
    check("Benign -> Pathogenic = upgraded",
          kind("Benign", "Pathogenic", REV0, REV0, None, None) == "upgraded")
    check("Pathogenic -> Benign = downgraded",
          kind("Pathogenic", "Benign", REV0, REV0, None, None) == "downgraded")
    check("Benign -> Likely_benign = reclassified",
          kind("Benign", "Likely_benign", REV0, REV0, None, None) == "reclassified")
    check("same sig, review ↑ on actionable = reclassified",
          kind("Pathogenic", "Pathogenic", REV0, REVX, None, None) == "reclassified")
    check("rare -> common AF crossing = af_shift",
          kind("Benign", "Benign", REV0, REV0, 0.0005, 0.02) == "af_shift")
    check("no change at all -> None",
          triage.change_kind("Benign", "Benign", REV0, REV0, 0.01, 0.01) is None)
    check("same actionable, no review change -> None",
          triage.change_kind("Pathogenic", "Pathogenic", REVX, REVX, None, None) is None)

    # tier scoring — a new pathogenic, rare, on-panel coding variant is actionable; a common
    # benign frequency drift is ignored; an unchanged variant yields no verdict.
    null_call = {"tier": "Likely pathogenic", "codes": ["PVS1"], "div": "", "csq": "stop_gained"}
    act = triage.score_delta(prev_sig="Uncertain_significance", curr_sig="Pathogenic",
                             prev_rev=REV0, curr_rev=REVX, prev_af=None, curr_af=None,
                             consequence="stop_gained", ref="C", alt="T",
                             acmg_call=null_call, in_panel=True, constrained=True)
    check("new pathogenic, rare, on-panel null -> actionable_now",
          act["tier"] == "actionable_now")
    drift = triage.score_delta(prev_sig="Benign", curr_sig="Benign", prev_rev=REV0,
                               curr_rev=REV0, prev_af=0.0005, curr_af=0.02,
                               consequence="intron_variant", ref="A", alt="G", in_panel=False)
    check("common benign AF drift -> ignore", drift["tier"] == "ignore")
    check("unchanged variant -> no verdict",
          triage.score_delta(prev_sig="Pathogenic", curr_sig="Pathogenic", prev_rev=REVX,
                             curr_rev=REVX, prev_af=None, curr_af=None) is None)

    # stability factor over a variant's cross-snapshot signature trail (Phase 2)
    check("monotone strengthening -> +0.5",
          triage.stability_factor(["Uncertain_significance", "Likely_pathogenic",
                                   "Pathogenic"]) == 0.5)
    check("oscillating call -> -0.5",
          triage.stability_factor(["Pathogenic", "Benign", "Pathogenic", "Benign",
                                   "Pathogenic"]) == -0.5)
    check("single data point -> neutral",
          triage.stability_factor(["Pathogenic"]) == 0.0)


def test_render() -> None:
    """The markdown→HTML converter: tables, escaping, emoji, self-containment."""
    print("render.to_html (markdown → self-contained HTML):")
    from pipeline import render
    md = "\n".join([
        "# Report — dad",
        "",
        "| Gene | Variant | ClinVar | ! |",
        "|------|---------|:-------:|---|",
        "| BRCA1 | 17:43000000 T>G | **Pathogenic** | 🔴 |",
        "",
        "- a `code` bullet",
        "_Research-grade, not a clinical diagnosis._",
    ])
    h = render.to_html(md, "Report — dad")
    check("doctype + table present", "<!DOCTYPE html>" in h and "<table>" in h)
    check("escapes '>' in T>G (load-bearing)", "T&gt;G" in h and "T>G" not in h)
    check("**bold** rendered", "<strong>Pathogenic</strong>" in h)
    check("`code` rendered", "<code>code</code>" in h)
    check("emoji passes through", "🔴" in h)
    check("h1 heading", "<h1>" in h)
    check("full-line _italic_ → note", 'class="note"' in h)
    check("self-contained (no external refs)",
          not any(x in h for x in ("http://", "https://", "src=", "<link")))
    # star-allele cells must not be mangled by the bold regex
    h2 = render.to_html("| Gene | Call |\n|---|---|\n| SLCO1B1 | *1/*15 |", "x")
    check("single-* star alleles survive", "*1/*15" in h2 and "<strong>" not in h2)


def test_forcecall_decode() -> None:
    """forcecall._decode: the load-bearing no-call-vs-evidenced-hom-ref distinction."""
    print("forcecall._decode (force-called genotype → alleles / no-call):")
    from pipeline import forcecall as fc
    d = fc._decode
    check("0/0 → evidenced hom-ref",
          d("A", "C", "0/0", "40") == {"ref": "A", "alt": "C", "gt": "0/0",
                                       "alleles": ["A", "A"], "dp": 40, "nocall": False})
    check("0/1 het → ref+alt", d("A", "C", "0/1", "30")["alleles"] == ["A", "C"])
    check("1/1 hom-alt → alt+alt", d("A", "C", "1/1", "30")["alleles"] == ["C", "C"])
    check("./. → no-call (NOT hom-ref)",
          d("A", "C", "./.", "0")["nocall"] and d("A", "C", "./.", "0")["alleles"] is None)
    check("phased 0|1 handled", d("A", "C", "0|1", "30")["alleles"] == ["A", "C"])
    check("multiallelic 1/2 → both alts", d("A", "C,T", "1/2", "30")["alleles"] == ["C", "T"])
    check("non-variant ALT '.' 0/0 → hom-ref", d("A", ".", "0/0", "30")["alleles"] == ["A", "A"])
    check("haploid '1' → homozygous", d("A", "C", "1", "30")["alleles"] == ["C", "C"])
    check("missing DP → 0", d("A", "C", "0/1", ".")["dp"] == 0)
    # Regression (bug review #12): a half-called ./1 was doubled into a confident hom-alt.
    check("partial no-call './1' → no-call (NOT confident hom-alt)",
          d("A", "C", "./1", "30")["nocall"] and d("A", "C", "./1", "30")["alleles"] is None)


def test_forcecall_sites() -> None:
    """forcecall._sites_records: REF from the reference, ALT = the marker's other allele."""
    print("forcecall._sites_records (build force-call records):")
    from pipeline import forcecall as fc
    orig_ref, orig_fai = fc._ref_bases, fc._fai_lengths
    # Stub the reference so the test needs no FASTA/.fai on disk.
    fc._fai_lengths = lambda: {"1": 1_000, "2": 1_000, "X": 1_000}
    fc._ref_bases = lambda pos: {
        ("1", "100"): "A", ("1", "200"): "G", ("2", "50"): "C",
        ("X", "10"): "T", ("1", "300"): "A", ("1", "400"): "A"}
    try:
        recs = fc._sites_records({
            ("1", "200"): {"G", "A"},      # ref G, effect A → ALT A
            ("1", "100"): {"A", "C"},      # ref A, effect C → ALT C
            ("2", "50"): {"C", "T"},       # ref C → ALT T
            ("X", "10"): {"T", "T"},       # both alleles == ref → dropped (no ALT)
            ("1", "300"): {"AT", "A"},     # indel allele → dropped (SNV-only)
            ("1", "400"): {"A", "N"},      # non-ACGT ALT → dropped
        })
    finally:
        fc._ref_bases, fc._fai_lengths = orig_ref, orig_fai
    check("drops ref-only / indel / non-ACGT sites", len(recs) == 3)
    check("REF is the reference base, ALT the other allele",
          recs == [("1", 100, "A", "C"), ("1", 200, "G", "A"), ("2", 50, "C", "T")])


def test_acmg_constraint_pvs1() -> None:
    """ACMG PVS1 grading by gnomAD constraint + the PVS1_moderate combine rules."""
    print("acmg._pvs1_code / _combine (constraint-graded PVS1):")
    from pipeline import acmg
    p, c = acmg._pvs1_code, acmg._combine
    check("recessive gene keeps full PVS1 even when LOF-tolerant",
          p(dominant=False, constrained=False) == "PVS1")
    check("dominant + LOF-constrained keeps full PVS1",
          p(dominant=True, constrained=True) == "PVS1")
    check("dominant + LOF-tolerant downgrades to PVS1_moderate",
          p(dominant=True, constrained=False) == "PVS1_moderate")
    check("no constraint data does not penalize (full PVS1)",
          p(dominant=True, constrained=None) == "PVS1")
    check("lone full PVS1 -> Likely pathogenic", c(["PVS1"]) == "Likely pathogenic")
    check("PVS1 + 2 supporting -> Pathogenic", c(["PVS1", "PM2", "PP3"]) == "Pathogenic")
    check("lone PVS1_moderate -> VUS", c(["PVS1_moderate"]) == "VUS")
    check("PVS1_moderate + supporting -> Likely pathogenic",
          c(["PVS1_moderate", "PM2"]) == "Likely pathogenic")
    check("PVS1_moderate + benign signal -> VUS", c(["PVS1_moderate", "BP4"]) == "VUS")
    check("PVS1_moderate can't reach Pathogenic",
          c(["PVS1_moderate", "PM2", "PP3"]) == "Likely pathogenic")


def test_gnomad_constraint() -> None:
    """gnomAD constraint: is_constrained thresholds + per-gene table reduction."""
    print("gnomad_constraint.is_constrained / setup._reduce_constraint:")
    from pipeline import gnomad_constraint as gc
    orig = gc._TABLE
    gc._TABLE = {"HI": (0.99, 0.21), "TOL": (0.0, 1.30), "LO": (0.0, 0.45),
                 "NA": (None, None)}
    try:
        check("high pLI -> constrained", gc.is_constrained("HI") is True)
        check("LOF-tolerant -> not constrained", gc.is_constrained("TOL") is False)
        check("low LOEUF (<0.6) -> constrained", gc.is_constrained("LO") is True)
        check("no metric -> None (distinct from False)", gc.is_constrained("NA") is None)
        check("gene absent from table -> None", gc.is_constrained("ZZZ") is None)
        check("label formats pLI + LOEUF", gc.label("HI") == "pLI 0.99 / LOEUF 0.21")
    finally:
        gc._TABLE = orig

    from pipeline import setup
    lines = [
        "gene\tcanonical\tmane_select\tlof.pLI\tlof.oe_ci.upper",
        "G1\tfalse\tfalse\t0.10\t1.20",      # lowest priority for G1
        "G1\ttrue\ttrue\t0.90\t0.30",        # MANE -> wins
        "G2\ttrue\tfalse\t0.50\t0.70",       # canonical -> wins (no MANE row)
        "G3\tfalse\tfalse\t0.05\t1.50",      # first-seen wins over equal-priority later
        "G3\tfalse\tfalse\t0.06\t1.40",
    ]
    out = setup._reduce_constraint(lines)
    check("one row per gene, sorted, best transcript kept",
          out == [("G1", "0.90", "0.30"), ("G2", "0.50", "0.70"), ("G3", "0.05", "1.50")])


def test_spliceai() -> None:
    """SpliceAI: INFO parsing, max-effect selection, bands/predicates, ACMG splice codes."""
    print("spliceai.lookup / bands / acmg splice PP3-BP4:")
    from pipeline import spliceai as sa

    # lookup(): parse the SpliceAI INFO field, pick the strongest effect across genes, match
    # ref/alt. Stub _query/_file_for so no real VCF or tabix is needed.
    line = ("7\t117559590\t.\tA\tG\t.\t.\t"
            "SpliceAI=G|CFTR|0.02|0.93|0.10|0.00|-5|2|-12|7,"
            "G|OTHER|0.01|0.10|0.00|0.00|1|2|3|4")
    orig_q, orig_f, orig_a = sa._query, sa._file_for, sa.available
    sa._query = lambda path, chrom, lo, hi: [line]
    sa._file_for = lambda ref, alt: "stub"
    sa.available = lambda: True
    try:
        hit = sa.lookup("7", 117559590, "A", "G")
        check("picks the maximal delta (acceptor loss 0.93)",
              hit is not None and abs(hit[0] - 0.93) < 1e-9 and hit[1] == "acceptor loss")
        check("returns the gene symbol of the top effect", hit[2] == "CFTR")
        check("ref/alt mismatch -> None", sa.lookup("7", 117559590, "A", "T") is None)

        # Regression (bug review #2): VEP-style '-' indel alleles must match the padded
        # VCF representation — previously every indel lookup silently returned None.
        ins_line = ("1\t100\t.\tA\tAT\t.\t.\t"
                    "SpliceAI=T|GENE|0.10|0.20|0.30|0.88|1|2|3|4")
        sa._query = lambda path, chrom, lo, hi: [ins_line] \
            if (lo, hi) == (100, 102) else []
        hit_ins = sa.lookup("1", 101, "-", "T")          # VEP insertion at 101
        check("VEP '-' insertion matches padded VCF record",
              hit_ins is not None and abs(hit_ins[0] - 0.88) < 1e-9)
        del_line = ("1\t100\t.\tATC\tA\t.\t.\t"
                    "SpliceAI=A|GENE|0.70|0.10|0.10|0.10|1|2|3|4")
        sa._query = lambda path, chrom, lo, hi: [del_line] \
            if (lo, hi) == (100, 102) else []
        hit_del = sa.lookup("1", 101, "TC", "-")        # VEP deletion at 101
        check("VEP '-' deletion matches padded VCF record",
              hit_del is not None and abs(hit_del[0] - 0.70) < 1e-9)
        check("_norm: padded ≡ VEP minimal (insertion)",
              sa._norm(100, "A", "AT") == sa._norm(101, "-", "T"))
        check("_norm: padded ≡ VEP minimal (deletion)",
              sa._norm(100, "ATC", "A") == sa._norm(101, "TC", "-"))
    finally:
        sa._query, sa._file_for, sa.available = orig_q, orig_f, orig_a

    # Bands + predicates keyed off the config cut-offs (0.2 / 0.5 / 0.8).
    hi, mid, lo = (0.93, "acceptor loss", "CFTR"), (0.6, "donor gain", "X"), (0.05, "donor loss", "X")
    check("≥0.8 -> high band", sa.label(hi).endswith("(high)"))
    check("≥0.5 -> likely band", sa.label(mid).endswith("(likely)"))
    check("<0.2 -> no-effect band", sa.label(lo).endswith("(no effect)"))
    check("None -> em dash", sa.label(None) == "—")
    check("is_pathogenic at recommended cut-off", sa.is_pathogenic(mid) and sa.is_pathogenic(hi))
    check("is_pathogenic False below cut-off", not sa.is_pathogenic(lo))
    check("is_benign only when confidently no effect", sa.is_benign(lo) and not sa.is_benign(mid))
    check("predicates on None are False", not sa.is_pathogenic(None) and not sa.is_benign(None))

    # ACMG integration: SpliceAI drives PP3/BP4 for near-splice variants and boosts missense.
    from pipeline import acmg
    splice_path = acmg._codes("splice_region_variant", None, "", "", None, hi, True)
    check("splice-altering near-splice variant earns PP3", "PP3" in splice_path)
    splice_ben = acmg._codes("splice_region_variant", None, "", "", None, lo, True)
    check("no-effect near-splice variant earns BP4", "BP4" in splice_ben)
    splice_none = acmg._codes("splice_region_variant", None, "", "", None, None, True)
    check("near-splice with no SpliceAI data earns neither", "PP3" not in splice_none
          and "BP4" not in splice_none)
    # A positive splice signal vetoes a missense benign in-silico call.
    veto = acmg._codes("missense_variant", None, "tolerated", "benign",
                       ("ambiguous-am", "likely_benign"), mid, True)
    check("splice signal vetoes missense BP4 -> PP3", "PP3" in veto and "BP4" not in veto)


def test_prs_genomewide_percentile() -> None:
    """Genome-wide PRS percentile: empirical ranking, population-map + .pvar parsing."""
    print("prs._empirical_percentile / _kg_population_map / _sample_covered:")
    import pathlib
    import tempfile
    from pipeline import config, prs

    ref = [("EUR", 1.0), ("EUR", 2.0), ("EUR", 3.0), ("EUR", 4.0),  # 4 EUR
           ("AFR", 5.0), ("EAS", 0.0)]                              # 2 non-EUR
    eur, allp = prs._empirical_percentile(2.5, ref)
    check("EUR percentile = share of EUR below (2/4 = 50%)", eur == 50.0)
    check("overall percentile counts all pops (3/6 = 50%)", allp == 50.0)
    top_eur, _ = prs._empirical_percentile(99.0, ref)
    check("score above everyone -> 100th EUR percentile", top_eur == 100.0)
    bot_eur, _ = prs._empirical_percentile(-1.0, ref)
    check("score below everyone -> 0th percentile", bot_eur == 0.0)
    check("empty reference -> None", prs._empirical_percentile(1.0, []) == (None, None))

    d = pathlib.Path(tempfile.mkdtemp())
    panel = d / "pop.txt"
    panel.write_text("SampleID\tFatherID\tMotherID\tSex\tPopulation\tSuperpopulation\n"
                     "HG001\t0\t0\t1\tGBR\tEUR\nNA002\t0\t0\t2\tYRI\tAFR\n")
    orig = config.KG_PANEL_FILE
    config.KG_PANEL_FILE = panel
    try:
        pm = prs._kg_population_map()
        check("population map keys sample -> superpop",
              pm == {"HG001": "EUR", "NA002": "AFR"})
    finally:
        config.KG_PANEL_FILE = orig

    pvar = d / "s"
    (d / "s.pvar").write_text("##fileformat\n#CHROM\tPOS\tID\tREF\tALT\n"
                              "1\t100\t1:100\tA\tG\n7\t250\t7:250\tC\tT\n")
    cov = prs._sample_covered(pvar)
    check("covered set parsed from .pvar (skips headers)",
          cov == {("1", "100"), ("7", "250")})

    # Regression (bug review #1): the plink2 --score file must use BARE 'chrom:pos' IDs —
    # plink2 strips any 'chr' prefix on pgen output, so a 'chr' prefix
    # here matched zero variants and silently killed the empirical percentile.
    sf = d / "score.txt"
    prs._write_score_file(sf, [("1", "100", "G", "A", 0.5), ("7", "250", "T", "C", 1.0)])
    ids = [ln.split("\t")[0] for ln in sf.read_text().splitlines()[1:]]
    check("score-file IDs are bare 'chrom:pos' (same namespace as .pvar)",
          ids == ["1:100", "7:250"])


def test_mito_heteroplasmy() -> None:
    """Mito: heteroplasmy classification, curated pathogenic match, Mutect2 VCF parse."""
    print("mito._classify / _curated / _parse_calls:")
    from pipeline import config, mito

    check("AF ≥ 0.95 -> homoplasmic", mito._classify(0.98) == "homoplasmic")
    check("AF < 0.95 -> heteroplasmic", mito._classify(0.42) == "heteroplasmic")

    cur = mito._curated()
    check("curated keyed by (pos,ref,alt)", cur.get(("3243", "A", "G"))["locus"] == "MT-TL1")
    check("curated distinguishes alt allele (8993 T>G vs T>C)",
          cur[("8993", "T", "G")]["mutation"] == "m.8993T>G"
          and cur[("8993", "T", "C")]["mutation"] == "m.8993T>C")
    check("non-pathogenic position absent from table", ("750", "A", "G") not in cur)

    # _parse_calls shells out to bcftools query — stub it with canned tab output.
    canned = ("3243\tA\tG\tPASS\t0.350\t4200\n"      # heteroplasmic pathogenic
              "750\tA\tG\tPASS\t0.998\t5100\n"        # homoplasmic polymorphism
              "8270\tC\tT\tweak_evidence\t0.04\t3800\n"  # filtered
              "malformed line without enough fields\n")

    class _R:
        stdout = canned
    orig = mito.subprocess.run
    mito.subprocess.run = lambda *a, **k: _R()
    try:
        calls = mito._parse_calls(pathlib.Path("ignored.vcf.gz"))
    finally:
        mito.subprocess.run = orig
    check("parses one row per well-formed record (skips malformed)", len(calls) == 3)
    c0 = calls[0]
    check("AF parsed as float heteroplasmy fraction", abs(c0["af"] - 0.350) < 1e-9)
    check("DP parsed as int", c0["dp"] == 4200)
    check("FILTER carried through", calls[2]["filter"] == "weak_evidence")


def test_hla_typing() -> None:
    """HLA typing: 4-digit truncation + risk-call re-derivation from arcasHLA genotype JSON."""
    print("hla_typing._four_digit / _risk_calls:")
    from pipeline import hla_typing as ht

    check("truncates 3+ field allele to 2-field",
          ht._four_digit("DRB1*15:01:01:02") == "DRB1*15:01")
    check("leaves a 2-field allele unchanged",
          ht._four_digit("DQA1*05:01") == "DQA1*05:01")

    b27_pos = {"B": ["B*27:05:02", "B*07:02:01"],
               "DQA1": ["DQA1*01:01"], "DQB1": ["DQB1*05:01"]}
    calls = ht._risk_calls(b27_pos)
    check("B27 present when a B*27:xx allele is typed", "present" in calls[0])
    check("celiac absent when neither DQ2.5 nor DQ8 combo is typed", "neither" in calls[1])

    dq25_dq8 = {"B": ["B*44:02"],
                "DQA1": ["DQA1*05:01", "DQA1*02:01"],
                "DQB1": ["DQB1*02:01", "DQB1*03:02"]}
    calls = ht._risk_calls(dq25_dq8)
    check("B27 absent when no B*27:xx allele is typed", "absent" in calls[0])
    check("DQ2.5 needs BOTH DQA1*05:01 and DQB1*02:01",
          "DQ2.5" in calls[1] and "DQ8" in calls[1])

    dq8_only = {"B": [], "DQA1": ["DQA1*03:01"], "DQB1": ["DQB1*03:02"]}
    calls = ht._risk_calls(dq8_only)
    check("DQ8 alone (no DQA1*05:01) doesn't also claim DQ2.5",
          ": DQ8 " in calls[1])
    check("empty B list reported as no B calls", "no B calls" in calls[0])


def _synthetic_family():
    from pipeline.pedigree import Family, Person
    people = {
        "grandfather": Person("grandfather", "Grandfather", "male", None, None, "grandmother"),
        "grandmother": Person("grandmother", "Grandmother", "female", None, None, "grandfather"),
        "father": Person("father", "Father", "male", "grandfather", "grandmother", "mother"),
        "maternal_grandmother": Person("maternal_grandmother", "Maternal Grandmother", "female", None, None, None),
        "mother": Person("mother", "Mother", "female", None, "maternal_grandmother", "father"),
        "daughter": Person("daughter", "Daughter", "female", "father", "mother", None),
    }
    return Family(people, sequenced={"grandfather", "grandmother", "father"})


def test_pedigree_relationships() -> None:
    """Perspective engine: kinship terms across generations, partners and in-laws."""
    print("pedigree.relationship_to / trio / siblings:")
    f = _synthetic_family()
    r = f.relationship_to
    check("parent -> father (sex-aware)", r("father", "grandfather") == "father")
    check("parent -> mother", r("father", "grandmother") == "mother")
    check("partner is symmetric + neutral", r("father", "mother") == "partner"
          and r("grandmother", "grandfather") == "partner")
    check("child -> daughter", r("father", "daughter") == "daughter")
    check("partner's parent -> mother-in-law", r("father", "maternal_grandmother") == "mother-in-law")
    check("reciprocal -> father-in-law", r("mother", "grandfather") == "father-in-law")
    check("grandparent", r("daughter", "grandfather") == "grandfather")
    check("grandchild", r("grandfather", "daughter") == "granddaughter")
    check("self", r("grandfather", "grandfather") == "self")
    check("unrelated/unknown -> 'relative'", r("grandfather", "Nobody") == "relative")
    check("trio resolves both parents", f.trio("daughter") == ("daughter", "father", "mother"))
    check("founder has no trio", f.trio("grandfather") is None)
    check("is_sequenced reflects the registry",
          f.is_sequenced("father") and not f.is_sequenced("mother"))


def test_pedigree_callsets() -> None:
    """Alternate callsets resolve to their person — the bug that made Adam his own
    brother, his daughters his nieces and his partner his sister-in-law."""
    print("pedigree.person_of / callsets / relationship_to across callsets:")
    from pipeline.pedigree import Family, Person
    people = {
        "Jan": Person("Jan", "Jan", "male", None, None, "Maria"),
        "Maria": Person("Maria", "Maria", "female", None, None, "Jan"),
        "Adam": Person("Adam", "Adam", "male", "Jan", "Maria", "Partner"),
        "Partner": Person("Partner", "Partner", "female", None, None, "Adam"),
        "Lena": Person("Lena", "Lena", "female", "Adam", "Partner", None),
    }
    f = Family(people, sequenced={"Jan", "Jan_t2t", "Maria", "Adam", "Adam_vendor",
                                  "Partner", "Lena", "Lena_vendor"})
    check("suffixed id resolves to the person", f.person_of("Adam_vendor") == "Adam")
    check("a real row is never split", f.person_of("Jan") == "Jan")
    check("an unknown id passes through", f.person_of("Someone_else") == "Someone_else")
    check("stacked suffix still resolves", f.person_of("Adam_vendor_rerun") == "Adam")
    check("callsets list the bare id first",
          f.callsets("Adam") == ["Adam", "Adam_vendor"])
    check("callsets are reachable from any of them",
          f.callsets("Adam_vendor") == f.callsets("Adam"))

    r = f.relationship_to
    check("a callset is the person, not a sibling", r("Adam", "Adam_vendor") == "self")
    check("a child's callset is still the child", r("Adam", "Lena_vendor") == "daughter")
    check("and the reverse", r("Lena_vendor", "Adam") == "father")
    check("a parent's callset is still the parent", r("Adam", "Jan_t2t") == "father")
    check("others() lists people, not callsets",
          [pid for pid, _ in f.others("Adam")] == ["Jan", "Maria", "Partner", "Lena"])
    check("a callset id resolves a trio, keeping its own identity",
          f.trio("Lena_vendor") == ("Lena_vendor", "Adam", "Partner"))
    check("a person with only a suffixed callset still counts as sequenced",
          f.is_sequenced("Jan"))
    check("validate() rejects a callset written in as a person",
          any("callset" in i for i in Family(
              {**people, "Adam_vendor": Person("Adam_vendor", "Adam (vendor)", "male",
                                                "Jan", "Maria", None)}).validate()))


def test_inheritance_origin() -> None:
    """Which parent a variant came from, and when the page must refuse to say."""
    print("inheritance.origin / discordant:")
    from pipeline import inheritance as inh
    from pipeline.pedigree import Family, Person
    people = {
        "Dad": Person("Dad", "Dad", "male", None, None, "Mum"),
        "Mum": Person("Mum", "Mum", "female", None, None, "Dad"),
        "Kid": Person("Kid", "Kid", "female", "Dad", "Mum", None),
    }
    seq = {"Dad", "Dad_t2t", "Mum", "Kid"}
    f = Family(people, sequenced=seq)

    def v(*carriers):
        return {"carriers": {c: {"zyg": "HET", "callsets": {c: "HET"}} for c in carriers}}

    check("from father", inh.origin(v("Dad", "Kid"), "Kid", f) == inh.FROM_FATHER)
    check("from mother", inh.origin(v("Mum", "Kid"), "Kid", f) == inh.FROM_MOTHER)
    check("from both", inh.origin(v("Dad", "Mum", "Kid"), "Kid", f) == inh.FROM_BOTH)
    check("neither parent carries it -> de-novo candidate",
          inh.origin(v("Kid"), "Kid", f) == inh.DE_NOVO)
    check("a founder's own variant has no traceable origin",
          inh.origin(v("Dad"), "Dad", f) == inh.FOUNDER)
    check("a non-carrier gets no origin at all", inh.origin(v("Dad"), "Kid", f) == "")

    # The claim "neither parent carries it" requires both parents to have been sequenced.
    unseq = Family(people, sequenced={"Kid", "Dad"})
    check("an unsequenced parent blocks the de-novo claim",
          inh.origin(v("Kid"), "Kid", unseq) == inh.UNRESOLVED)

    partial = {"carriers": {"Dad": {"zyg": "HET", "callsets": {"Dad": "HET"}}}}
    check("one of two callsets calling it is flagged as disagreement",
          inh.discordant(partial, f) == ["Dad"])
    both = {"carriers": {"Dad": {"zyg": "HET",
                                 "callsets": {"Dad": "HET", "Dad_t2t": "HET"}}}}
    check("full callset support is not flagged", inh.discordant(both, f) == [])


def test_pedigree_ped_generation() -> None:
    """family graph -> standard PED (parent columns + sex codes)."""
    print("pedigree.generate_pedigree_ped:")
    import pathlib
    import tempfile
    from pipeline import pedigree
    f = _synthetic_family()
    out = pathlib.Path(tempfile.mktemp(suffix=".ped"))
    pedigree.generate_pedigree_ped(f, out)
    rows = {ln.split("\t")[1]: ln.split("\t") for ln in out.read_text().splitlines()}
    check("founder has 0/0 parents + male code", rows["grandfather"][2:5] == ["0", "0", "1"])
    check("child references parent ids + female code",
          rows["daughter"][2:6] == ["father", "mother", "2", "0"])
    check("one row per person", len(rows) == 6)


def test_family_migrate_plan() -> None:
    """Rename planner: token-aware, deepest-first, applied correctly on a synthetic tree."""
    print("family_migrate.retoken / plan_renames:")
    import pathlib
    import tempfile
    from pipeline import family_migrate as fm

    m = {"father": "SAMPLE001"}
    check("token in filename renamed", fm.retoken("mito_father.md", m) == "mito_SAMPLE001.md")
    check("bare-name file renamed", fm.retoken("father.cram", m) == "SAMPLE001.cram")
    check("substring NOT touched (grandfather)",
          fm.retoken("grandfather.txt", m) == "grandfather.txt")

    root = pathlib.Path(tempfile.mkdtemp())
    (root / "father.cram").write_text("x")
    (root / "father.cram.crai").write_text("x")
    (root / "mito" / "father").mkdir(parents=True)
    (root / "mito" / "father" / "mito_father.md").write_text("x")
    (root / "reports").mkdir()
    (root / "reports" / "index_father.html").write_text("x")

    ops = fm.plan_renames(m, roots=[root])
    names = {(s.name, d.name) for s, d in ops}
    check("plans the bare CRAM + index + nested file + dir", {
        ("father.cram", "SAMPLE001.cram"), ("father.cram.crai", "SAMPLE001.cram.crai"),
        ("mito_father.md", "mito_SAMPLE001.md"), ("father", "SAMPLE001"),
        ("index_father.html", "index_SAMPLE001.html")}.issubset(names))
    # deepest-first: the nested file must be renamed before its parent dir
    paths = [s for s, _ in ops]
    fi = next(i for i, p in enumerate(paths) if p.name == "mito_father.md")
    di = next(i for i, p in enumerate(paths) if p.name == "father" and p.parent.name == "mito")
    check("inner file ordered before its parent dir", fi < di)

    for src, dst in ops:                              # apply, then verify the final tree
        src.rename(dst)
    check("renamed CRAM exists, old gone",
          (root / "SAMPLE001.cram").exists() and not (root / "father.cram").exists())
    check("nested file landed under renamed dir",
          (root / "mito" / "SAMPLE001" / "mito_SAMPLE001.md").exists())


def test_family_migrate_manifests() -> None:
    """Rename must follow path columns too, not just the id column: plan_renames moves the
    files, so a manifest whose raw_path still names the old id points at nothing."""
    print("family_migrate._apply_tsv path-column rewriting:")
    import pathlib
    import tempfile
    from pipeline import family_migrate as fm

    m = {"father": "SAMPLE001", "mother": "SAMPLE002"}
    tsv = pathlib.Path(tempfile.mkdtemp()) / "samples.tsv"
    tsv.write_text(
        "sample_id\trelationship\tsource\traw_path\n"
        "father\tfather\t/g/called/father.GRCh38.vcf.gz\t/g/raw/father/father.GRCh38.vcf.gz\n"
        "SAMPLE_t2t\tfather\t/g/t2t/SAMPLE.t2t.vcf.gz\t/g/raw/SAMPLE_t2t/SAMPLE.t2t.vcf.gz\n")
    fm._apply_tsv(tsv, "sample_id", ("source", "raw_path"), (), m)
    rows = [ln.split("\t") for ln in tsv.read_text().splitlines()[1:]]

    check("id column rewritten", rows[0][0] == "SAMPLE001")
    check("every path occurrence rewritten",
          rows[0][3] == "/g/raw/SAMPLE001/SAMPLE001.GRCh38.vcf.gz")
    check("source column rewritten", rows[0][2] == "/g/called/SAMPLE001.GRCh38.vcf.gz")
    check("already-neutral row untouched (idempotent)",
          rows[1] == ["SAMPLE_t2t", "father", "/g/t2t/SAMPLE.t2t.vcf.gz",
                      "/g/raw/SAMPLE_t2t/SAMPLE.t2t.vcf.gz"])

    cols = fm._column_indexes(["sample_id", "raw_path"], "sample_id", ("raw_path", "absent"))
    check("path column missing from an older header is skipped", cols == (0, (1,), ()))
    check("_rewrite_row reports no-change on a neutral row",
          fm._rewrite_row(["SAMPLE002", "/g/raw/SAMPLE002/SAMPLE002.vcf.gz"], cols, m) is False)

    # family.tsv is the family GRAPH and pedigree.ped is regenerated from it at the end of
    # migrate(), so it must be rewritten too — and its father/mother/partner cells hold
    # OTHER people's ids, which need the same exact-match remap as the id column.
    fam = pathlib.Path(tempfile.mkdtemp()) / "family.tsv"
    fam.write_text("id\tdisplay_name\tsex\tfather\tmother\tpartner\n"
                   "father\tFather\tmale\t\t\tmother\n"
                   "mother\tMother\tfemale\t\t\tfather\n"
                   "kid\tChild\tmale\tfather\tmother\t\n")
    fm._apply_tsv(fam, "id", (), fm._FAMILY_ID_REFS, m)
    frows = [ln.split("\t") for ln in fam.read_text().splitlines()[1:]]
    check("family.tsv id column rewritten",
          [r[0] for r in frows] == ["SAMPLE001", "SAMPLE002", "kid"])
    check("family.tsv father/mother references follow the rename",
          frows[2][3] == "SAMPLE001" and frows[2][4] == "SAMPLE002")
    check("family.tsv partner references follow the rename",
          frows[0][5] == "SAMPLE002" and frows[1][5] == "SAMPLE001")
    check("family.tsv is in the migration manifest list",
          any(path == config.FAMILY_FILE for path, *_ in fm._MANIFEST_COLUMNS))


def test_segregation() -> None:
    """M2 segregation: a cis/trans verdict is built entirely out of 'this parent does NOT
    carry it' claims, so an unconfirmed absence must never be allowed to produce a verdict."""
    print("segregation.transmission / resolve:")
    from pipeline import segregation as sg

    check("only the father carries it -> paternal",
          sg.transmission("carries", "absent") == sg.PATERNAL)
    check("only the mother carries it -> maternal",
          sg.transmission("absent", "carries") == sg.MATERNAL)
    check("neither carries it -> de-novo candidate",
          sg.transmission("absent", "absent") == sg.DE_NOVO)
    check("both carry it -> ambiguous by transmission",
          sg.transmission("carries", "carries") == sg.BOTH)
    check("an unknown parent is UNRESOLVED, never inferred from the other",
          sg.transmission("carries", "unknown") == sg.UNRESOLVED)
    check("unknown on the paternal side too",
          sg.transmission("unknown", "absent") == sg.UNRESOLVED)

    def gene(*tr):
        return [{"transmission": t} for t in tr]

    v, _ = sg.resolve(gene(sg.PATERNAL, sg.MATERNAL))
    check("one from each parent -> in-trans (both copies hit)", v == sg.IN_TRANS)
    v, _ = sg.resolve(gene(sg.PATERNAL, sg.PATERNAL))
    check("both from one parent -> in-cis (one copy intact)", v == sg.IN_CIS)
    v, _ = sg.resolve(gene(sg.PATERNAL, sg.UNRESOLVED))
    check("one unresolved variant makes the gene unresolved", v == sg.UNRESOLVED)
    v, why = sg.resolve(gene(sg.PATERNAL, sg.BOTH))
    check("a variant both parents carry can't be phased by transmission",
          v == sg.UNRESOLVED and "cannot be told apart" in why)
    v, why = sg.resolve(gene(sg.MATERNAL, sg.DE_NOVO))
    check("a de-novo member defers rather than assuming a copy",
          v == sg.UNRESOLVED and "denovo" in why)


def test_segregation_indel_matching() -> None:
    """A candidate is matched back to the child's VCF to know WHERE to force-call the
    parents. Matching indels on "is an indel" alone binds whichever record `bcftools norm
    -m -both` happened to emit first at that position — and the parental absence that the
    cis/trans verdict rests on is then evidence about a different variant."""
    print("segregation._same_event (VEP <-> VCF allele reconstruction):")
    from pipeline import segregation as sg

    # VEP deletion 'T'/'-' at p == VCF AT>A at p-1.
    dele = {"ref": "T", "alt": "-"}
    check("deletion matches its own anchored VCF form", sg._same_event(dele, "AT", "A"))
    check("deletion does NOT match a different deletion at the same position",
          not sg._same_event(dele, "AG", "A"))
    check("deletion does NOT match a longer deletion sharing its first base",
          not sg._same_event(dele, "ATG", "A"))
    check("deletion does NOT match an insertion at the same position",
          not sg._same_event(dele, "A", "AT"))

    # VEP insertion '-'/'TG' at p == VCF A>ATG at p-1.
    ins = {"ref": "-", "alt": "TG"}
    check("insertion matches its own anchored VCF form", sg._same_event(ins, "A", "ATG"))
    check("insertion does NOT match a different inserted sequence",
          not sg._same_event(ins, "A", "ATC"))
    check("insertion does NOT match a deletion at the same position",
          not sg._same_event(ins, "ATG", "A"))

    # delins: VCF ATT>AG at p is VEP TT>G at p+1 — used to take the exact-match branch,
    # which can never match the anchored form, so these were silently dropped.
    check("delins matches its anchored VCF form",
          sg._same_event({"ref": "TT", "alt": "G"}, "ATT", "AG"))
    check("delins does NOT match a different delins",
          not sg._same_event({"ref": "TT", "alt": "G"}, "ATT", "AC"))

    check("case differences don't defeat the match",
          sg._same_event({"ref": "t", "alt": "-"}, "at", "A"))


def test_prs_ref_effect_supplement() -> None:
    """The genome-wide percentile compared a sample against a reference panel scored only at
    the sample's OWN variant sites — ascertainment bias that pushed every sample toward the
    top percentile. Both sides now score the full variant set, which requires adding back the
    hom-ref dosage a variant-only VCF cannot express."""
    print("prs._ref_effect_supplement (variant-only hom-ref dosage):")
    import pathlib
    import tempfile
    from pipeline import prs

    wd = pathlib.Path(tempfile.mkdtemp())
    (wd / "chr1.pvar").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "1\t100\t1:100\tA\tG\n"       # effect allele A == REF
        "1\t200\t1:200\tC\tT\n"       # effect allele T == ALT
        "1\t300\t1:300\tG\tA\n")      # effect allele G == REF, but sample HAS a call here
    variants = [("1", "100", "A", "G", 1.5),      # absent + ref-effect  -> owed 2*1.5
                ("1", "200", "T", "C", 2.0),      # absent + alt-effect  -> owed nothing
                ("1", "300", "G", "A", 5.0),      # covered              -> plink2 has it
                ("1", "400", "A", "G", 9.0)]      # absent, REF unknown  -> skipped
    covered = {("1", "300")}

    supp, n = prs._ref_effect_supplement(variants, covered, wd)
    check("only the absent reference-effect allele is supplemented", n == 1)
    check("supplement is 2 x weight (hom-ref = two copies)", abs(supp - 3.0) < 1e-9)

    supp_all_covered, n2 = prs._ref_effect_supplement(
        variants, {("1", "100"), ("1", "200"), ("1", "300"), ("1", "400")}, wd)
    check("a fully-covered sample is owed nothing",
          n2 == 0 and supp_all_covered == 0.0)

    check("REF lookup only returns positions the score asks for",
          set(prs._ref_alleles(variants, wd)) == {("1", "100"), ("1", "200"), ("1", "300")})


def test_sample_vcf_fallback() -> None:
    """An INGESTED sample (a finished VCF from another lab) has no called/<s>.GRCh38.vcf.gz,
    which used to make traits/HLA/PRS/ancestry refuse to run at all. The fallback must be the
    NORMALIZED vcf, never the raw ingested one: raw carries its source's contig naming, and
    every one of those engines treats an unmatched position as hom-ref — so reading raw would
    report confident wrong genotypes rather than failing."""
    print("call.sample_vcf (ingested-sample fallback):")
    import pathlib
    import tempfile
    from pipeline import call as call_mod, config, normalize as norm_mod

    root = pathlib.Path(tempfile.mkdtemp())
    called, normd = root / "called", root / "normalized"
    called.mkdir(); normd.mkdir()
    o_called, o_norm = config.CALLED_DIR, config.NORMALIZED_DIR
    try:
        config.CALLED_DIR, config.NORMALIZED_DIR = called, normd
        (normd / f"Ingested.norm.vcf.gz").write_text("x")
        check("ingested sample falls back to the normalized VCF",
              call_mod.sample_vcf("Ingested") == norm_mod.normalized_path("Ingested"))

        (called / f"Own.{config.BUILD}.vcf.gz").write_text("x")
        (normd / f"Own.norm.vcf.gz").write_text("x")
        check("a self-called sample still prefers its own called VCF",
              call_mod.sample_vcf("Own") == call_mod.called_vcf("Own"))

        check("neither present -> reports the called path (the actionable error)",
              call_mod.sample_vcf("Absent") == call_mod.called_vcf("Absent"))
    finally:
        config.CALLED_DIR, config.NORMALIZED_DIR = o_called, o_norm


def test_denovo() -> None:
    """M2 trio de-novo: the candidate filter and — the part that actually matters — the
    parental verdict, where 'absent from the parent's VCF' must never by itself pass."""
    print("denovo.child_passes / parent_verdict:")
    from pipeline import config, denovo as dn

    med = 40
    max_dp = int(config.DN_MAX_DP_FACTOR * med)          # 100x at a 40x median
    ok, _ = dn.child_passes(dp=40, gq=50, vaf=0.5, max_dp=max_dp)
    check("balanced het at ordinary depth is a candidate", ok)
    check("low GQ rejected", dn.child_passes(40, 5, 0.5, max_dp) == (False, "low_gq"))
    check("low depth rejected", dn.child_passes(4, 50, 0.5, max_dp) == (False, "low_depth"))
    check("collapsed-paralog depth rejected",
          dn.child_passes(300, 50, 0.5, max_dp) == (False, "excess_depth"))
    check("mosaic-range VAF rejected", dn.child_passes(40, 50, 0.08, max_dp)[1] == "low_vaf")
    check("near-hom VAF rejected", dn.child_passes(40, 50, 0.95, max_dp)[1] == "high_vaf")

    check("AD parses to an alt fraction", abs(dn._vaf("20,20") - 0.5) < 1e-9)
    check("unparseable AD is not silently 0.5", dn._vaf(".") == 0.0)

    # The core rule: only evidenced hom-ref at depth supports a de-novo call.
    check("evidenced hom-ref at depth confirms",
          dn.parent_verdict("0/0", 40, 0)[0] == dn.CONFIRMED)
    check("parental no-call is indeterminate, NOT de-novo",
          dn.parent_verdict("./.", 40, 0)[0] == dn.INDETERMINATE)
    check("missing record (no force-call result) is indeterminate, NOT de-novo",
          dn.parent_verdict("", 0, 0)[0] == dn.INDETERMINATE)
    check("hom-ref below the depth floor is indeterminate",
          dn.parent_verdict("0/0", 3, 0)[0] == dn.INDETERMINATE)
    check("parent carrying the allele is inherited",
          dn.parent_verdict("0/1", 40, 20)[0] == dn.INHERITED)
    check("one stray alt read is tolerated (sequencing error)",
          dn.parent_verdict("0/0", 40, 1)[0] == dn.CONFIRMED)
    check("several alt reads read as transmission, not error",
          dn.parent_verdict("0/0", 40, 9)[0] == dn.INHERITED)

    # Positional AD parsing: a missing per-allele value must not shift later alleles down.
    check("AD with a missing middle allele keeps positions",
          dn._ad_counts("10,.,5") == [10, 0, 5])
    check("VAF from a gappy AD uses the right slots",
          abs(dn._vaf("10,.,5") - 5 / 15) < 1e-9)

    # Natural contig ordering (the report used to list 1, 10, 11, …, 2).
    check("chr10 sorts after chr9",
          dn._chrom_key("9") < dn._chrom_key("10"))
    check("X sorts after the autosomes",
          dn._chrom_key("22") < dn._chrom_key("X"))

    # Cache keying: a stale artefact must be impossible to serve, and must be cleaned up.
    import pathlib, tempfile
    wd = pathlib.Path(tempfile.mkdtemp())
    a = dn._keyed_path(wd, "Dad.pileup", "aaaaaaaaaaaa", ".tsv")
    a.write_text("x")
    (wd / "Dad.pileup.tsv").write_text("pre-keying leftover")
    b = dn._keyed_path(wd, "Dad.pileup", "bbbbbbbbbbbb", ".tsv")
    check("a different key yields a different path", a != b)
    check("the stale keyed artefact is pruned", not a.exists())
    check("the pre-keying unsuffixed artefact is pruned too",
          not (wd / "Dad.pileup.tsv").exists())

    # Clustering: independent point mutations don't arrive in bursts.
    near = [{"chrom": "1", "pos": 1000, "status": dn.CONFIRMED},
            {"chrom": "1", "pos": 1400, "status": dn.CONFIRMED},
            {"chrom": "1", "pos": 900000, "status": dn.CONFIRMED},
            {"chrom": "2", "pos": 1200, "status": dn.CONFIRMED},
            {"chrom": "1", "pos": 1100, "status": dn.INHERITED}]
    dn._flag_clusters(near)
    check("two calls within the window are both demoted",
          [c["status"] for c in near[:2]] == [dn.CLUSTERED, dn.CLUSTERED])
    check("an isolated call survives", near[2]["status"] == dn.CONFIRMED)
    check("same position on another contig is not a cluster",
          near[3]["status"] == dn.CONFIRMED)
    check("an already-inherited call is not reclassified",
          near[4]["status"] == dn.INHERITED)

    # Multi-callset exclusion: one person called twice is one genome, two samples.tsv rows.
    import tempfile, pathlib
    tsv = pathlib.Path(tempfile.mkdtemp()) / "samples.tsv"
    tsv.write_text("sample_id\trelationship\n"
                   "father\tfather\nfather_t2t\tfather\nmother\tmother\nfatherly\tsister\n")
    orig_tsv, orig_np = config.SAMPLES_TSV, dn.normalized_path
    try:
        config.SAMPLES_TSV = tsv
        dn.normalized_path = lambda s: pathlib.Path("/dev/null")   # "exists" for every id
        got = dn.parent_callsets("father")
        check("parent's own callsets collected", set(got) == {"father", "father_t2t"})
        check("a different person with a prefix-sharing name is NOT collected",
              "fatherly" not in got)
    finally:
        config.SAMPLES_TSV, dn.normalized_path = orig_tsv, orig_np

    # A candidate is only de-novo when BOTH parents give positive evidence.
    def status(v1, v2):
        s = {v1[0], v2[0]}
        return (dn.INHERITED if dn.INHERITED in s
                else dn.INDETERMINATE if dn.INDETERMINATE in s else dn.CONFIRMED)
    check("one confirming + one indeterminate parent is NOT de-novo",
          status(dn.parent_verdict("0/0", 40, 0),
                 dn.parent_verdict("./.", 0, 0)) == dn.INDETERMINATE)
    check("both parents evidenced hom-ref is de-novo",
          status(dn.parent_verdict("0/0", 40, 0),
                 dn.parent_verdict("0/0", 35, 0)) == dn.CONFIRMED)


def test_ingest_does_not_add_callsets_as_people() -> None:
    """A callset id must not become a row in family.tsv.

    `_upsert_family` keyed on the exact id, so ingesting `Jan_t2t` appended a second,
    parentless Jan — the duplicate family tree that makes someone their own sibling."""
    print("ingest._upsert_family (callset vs person):")
    import pathlib
    import tempfile
    from pipeline import config, ingest
    from pipeline.util import read_tsv

    root = pathlib.Path(tempfile.mkdtemp())
    saved = {n: getattr(config, n) for n in ("FAMILY_FILE", "PEDIGREE_FILE", "SAMPLES_TSV")}
    try:
        config.FAMILY_FILE = root / "family.tsv"
        config.PEDIGREE_FILE = root / "pedigree.ped"
        config.SAMPLES_TSV = root / "samples.tsv"
        config.SAMPLES_TSV.write_text("sample_id\nJan\n")
        config.FAMILY_FILE.write_text(
            "id\tdisplay_name\tsex\tfather\tmother\tpartner\n"
            "Jan\tJan\tmale\t\t\t\n")

        ingest._upsert_family("Jan_t2t", "male")
        ids = [r["id"] for r in read_tsv(config.FAMILY_FILE)]
        check("a callset of a known person adds no row", ids == ["Jan"])

        ingest._upsert_family("Lena", "female")
        ids = [r["id"] for r in read_tsv(config.FAMILY_FILE)]
        check("a genuinely new person is still added", ids == ["Jan", "Lena"])

        ingest._upsert_family("Jan", "male")
        check("re-ingesting an existing person does not duplicate them",
              [r["id"] for r in read_tsv(config.FAMILY_FILE)] == ["Jan", "Lena"])
    finally:
        for n, v in saved.items():
            setattr(config, n, v)


def test_consent_is_per_person() -> None:
    """A grant belongs to the person, not to one of their callsets.

    Keying on the exact sample_id meant `Jan` was disclosed while `Jan_t2t` — the
    same genome called against T2T — was withheld, so one person's two reports disagreed
    about what he was allowed to be told."""
    print("consent.has_consented / set_consent across callsets:")
    import pathlib
    import tempfile
    from pipeline import config, consent
    from pipeline.util import read_tsv

    root = pathlib.Path(tempfile.mkdtemp())
    saved = {n: getattr(config, n) for n in ("CONSENT_FILE", "FAMILY_FILE",
                                             "SAMPLES_TSV", "REPEATS_DIR", "REPORTS_DIR")}
    try:
        config.CONSENT_FILE = root / "consent.tsv"
        config.FAMILY_FILE = root / "family.tsv"
        config.SAMPLES_TSV = root / "samples.tsv"
        config.REPEATS_DIR = root / "repeats"
        config.REPORTS_DIR = root / "reports"
        config.FAMILY_FILE.write_text(
            "id\tdisplay_name\tsex\tfather\tmother\tpartner\n"
            "Jan\tJan\tmale\t\t\t\n"
            "Maria\tMaria\tfemale\t\t\t\n")
        config.SAMPLES_TSV.write_text(
            "sample_id\n" + "\n".join(["Jan", "Jan_t2t", "Maria"]) + "\n")
        cat = consent.CATEGORY_ADULT_ONSET_UNTREATABLE

        consent.set_consent("Jan", cat, True)
        check("the grant covers the person's other callsets",
              consent.has_consented("Jan_t2t", cat) is True)
        check("and does not leak to anyone else",
              consent.has_consented("Maria", cat) is False)

        # A ledger written before the rule changed keys on a callset id.
        config.CONSENT_FILE.write_text(
            f"sample_id\tcategory\tconsent\nJan_t2t\t{cat}\tyes\n")
        check("a legacy callset row is read as that person's grant",
              consent.has_consented("Jan", cat) is True)

        # Person row and callset row disagree: the person's own row is the decision.
        config.CONSENT_FILE.write_text(
            f"sample_id\tcategory\tconsent\nJan\t{cat}\tno\n"
            f"Jan_t2t\t{cat}\tyes\n")
        check("the person's own row wins over a callset's",
              consent.has_consented("Jan_t2t", cat) is False)

        # Granting through a callset id records the PERSON and folds the stale rows in.
        consent.set_consent("Jan_t2t", cat, True)
        rows = [r for r in read_tsv(config.CONSENT_FILE) if r["category"] == cat]
        check("one row per person, recorded under the person's id",
              [r["sample_id"] for r in rows] == ["Jan"])
        check("granted through a callset id, consented as the person",
              consent.has_consented("Jan", cat) is True)

        # Revocation must reach every callset's disclosing report, not just one.
        from pipeline import reportpaths
        mds = [reportpaths.genome_report("repeats", cs) for cs in ("Jan", "Jan_t2t")]
        check("both callsets' reports are filed under the one person",
              {m.parent.name for m in mds} == {"Jan"} and len(set(mds)) == 2)
        for md in mds:
            md.parent.mkdir(parents=True, exist_ok=True)
            md.write_text("HTT 41 repeats")
        consent.set_consent("Jan", cat, False)
        check("revoking removes the disclosing report of EVERY callset",
              not any(md.exists() for md in mds))
    finally:
        for n, v in saved.items():
            setattr(config, n, v)


def test_consent() -> None:
    """S1: consent.tsv round-trip + repeats.py's withhold-unless-consented gate."""
    print("consent.has_consented / set_consent, repeats withholding:")
    import pathlib
    import tempfile
    from pipeline import config, consent

    from pipeline.util import read_tsv

    orig_file = config.CONSENT_FILE
    config.CONSENT_FILE = pathlib.Path(tempfile.mkdtemp()) / "consent.tsv"
    try:
        cat = consent.CATEGORY_ADULT_ONSET_UNTREATABLE
        check("no file yet -> default not consented", consent.has_consented("SAMPLE001", cat) is False)
        consent.set_consent("SAMPLE001", cat, True)
        check("granted -> consented", consent.has_consented("SAMPLE001", cat) is True)
        check("other sample unaffected", consent.has_consented("SAMPLE002", cat) is False)
        consent.set_consent("SAMPLE001", cat, False)
        check("revoked -> not consented again", consent.has_consented("SAMPLE001", cat) is False)
        check("upsert kept exactly one row for the (sample, category) pair",
              sum(1 for r in read_tsv(config.CONSENT_FILE)
                  if r["sample_id"] == "SAMPLE001" and r["category"] == cat) == 1)
    finally:
        config.CONSENT_FILE = orig_file

    # repeats.py's REAL consent gate (was: an inlined copy of the logic plus a vacuous
    # `r not in []` check that could never fail).
    from pipeline import repeats
    flagged = [("HTT", "CAG", [40], "⚠️ EXPANDED — Huntington disease (pathogenic range, ≥36)", True)]
    normal = [("FMR1", "CGG", [30], "normal (Fragile X syndrome)", False)]
    other = [("XYZ1", "CAG", [12], "—", False)]        # no disease interpretation
    rows = normal + flagged + other
    visible, withheld = repeats._apply_consent(rows, consented=False)
    # The WHOLE incidental category is withheld, normal results included. Withholding only
    # the flagged rows named the finding by elimination: _DISORDERS is a fixed public table,
    # so "HTT missing from a table listing every other locus" IS the diagnosis.
    check("unconsented: flagged incidental locus withheld",
          flagged[0] in withheld and flagged[0] not in visible)
    check("unconsented: NORMAL incidental locus also withheld (no elimination leak)",
          normal[0] in withheld and normal[0] not in visible)
    check("unconsented: non-incidental locus stays visible",
          visible == other)
    visible, withheld = repeats._apply_consent(rows, consented=True)
    check("consented: nothing withheld", visible == rows and withheld == [])

    # The unconsented report must look identical whether or not anything was flagged —
    # otherwise the withheld COUNT alone discloses "you carry an expansion".
    _, withheld_clean = repeats._apply_consent(normal + other, consented=False)
    check("withheld count is independent of the findings",
          len(withheld_clean) == len([r for r in rows if repeats.is_incidental(r[0])]) - 1
          and len(withheld_clean) == 1)

    # ExpansionHunter's GRCh38 catalog spells this locus "C9ORF72"; _DISORDERS keys it
    # "C9orf72". Exact matching left the commonest genetic cause of ALS/FTD genotyped but
    # uninterpreted — status "—", never flagged at any repeat size, and never withheld.
    check("catalog-cased C9ORF72 resolves to its disorder",
          repeats._disorder("C9ORF72") is not None)
    status, flag = repeats._status("C9ORF72", [60])       # well past the ≥30 threshold
    check("a pathogenic C9ORF72 expansion is flagged, not '—'",
          flag is True and "EXPANDED" in status)
    check("catalog-cased C9ORF72 is inside the S1 incidental category",
          repeats.is_incidental("C9ORF72"))
    check("a locus with no interpretation is still not incidental",
          not repeats.is_incidental("XYZ1"))


def test_penetrance() -> None:
    """S2: curated per-gene note surfaces on a finding row; uncurated genes fall back
    to the generic caveat; neither overstates certainty."""
    print("penetrance.note_for + panel._row_md / panel_report caveat wiring:")
    from pipeline import penetrance as pen
    from pipeline.util import read_tsv

    orig = pen._TABLE
    pen._TABLE = {"BRCA1": "~55-72% lifetime breast cancer risk, ~39-44% ovarian — "
                            "elevated, not certain"}
    try:
        check("curated gene -> its specific note",
              pen.note_for(["BRCA1"]) == pen._TABLE["BRCA1"])
        check("uncurated gene -> None (no note fabricated)",
              pen.note_for(["ZZZTOP"]) is None)
        check("multi-gene disease: first curated match wins",
              pen.note_for(["ZZZTOP", "BRCA1"]) == pen._TABLE["BRCA1"])

        dz = panel._Disease("Hereditary breast and ovarian cancer", ["BRCA1"], "AD", "")
        hits = [_v("BRCA1", ("17", 1, "A", "T"), "HET")]
        row = panel._row_md(dz, panel._ATRISK, hits, primary=True)
        check("curated finding row embeds the specific note", pen._TABLE["BRCA1"] in row)
        check("curated note itself doesn't claim certainty",
              "not certain" in pen._TABLE["BRCA1"] or "incomplete" in pen._TABLE["BRCA1"])

        dz2 = panel._Disease("Some uncurated condition", ["ZZZTOP"], "AD", "")
        hits2 = [_v("ZZZTOP", ("1", 1, "A", "T"), "HET")]
        row2 = panel._row_md(dz2, panel._ATRISK, hits2, primary=True)
        check("uncurated finding row defers to the general caveat, not silence",
              "general caveat" in row2)
    finally:
        pen._TABLE = orig

    caveat = panel._GENERIC_PENETRANCE_CAVEAT
    check("generic caveat says elevated risk, not certainty",
          "not certainty" in caveat or "not certain" in caveat)
    check("generic caveat doesn't imply nothing to worry about",
          "discount the finding" in caveat)

    # A quoted lifetime-risk range must say where the number came from. The published
    # BRCA1/2 and Lynch ranges are from cohorts recruited through clinically affected
    # families, so an untagged range would read as this person's risk when, for anyone not
    # referred that way, it is the upper end of the literature.
    table = {r["gene"]: r["note"] for r in read_tsv(pathlib.Path("panels/penetrance.tsv"))}
    quoted = {g: n for g, n in table.items() if "%" in n}
    check("the shipped table still quotes ranges for the genes that have them",
          {"BRCA1", "BRCA2", "MLH1", "MSH2"} <= set(quoted))
    check("every quoted range is tagged with its ascertainment",
          all(pen.is_ascertained(n) for n in quoted.values()))
    check("notes without a number are not tagged",
          not any(pen.is_ascertained(n) for g, n in table.items() if g not in quoted))
    check("is_ascertained is false for an absent note", not pen.is_ascertained(None))
    note = pen.ASCERTAINMENT_NOTE
    check("the explanation says the tagged figure runs high",
          "lower" in note and "not as this person's risk" in note)

    orig2 = pen._TABLE
    pen._TABLE = {"BRCA1": table["BRCA1"], "ZZZTOP2": "qualitative only, no number"}
    try:
        dzs = [(panel._Disease("HBOC", ["BRCA1"], "AD", ""), panel._ATRISK,
                [_v("BRCA1", ("17", 1, "A", "T"), "HET")])]
        notes = [pen.note_for(sorted({v["gene"] for v in h})) for _d, _s, h in dzs]
        check("a tagged gene on the report triggers the explanation",
              any(pen.is_ascertained(n) for n in notes))
        untagged = [pen.note_for(["ZZZTOP2"])]
        check("an untagged-only report does not carry it",
              not any(pen.is_ascertained(n) for n in untagged))
    finally:
        pen._TABLE = orig2


def test_incremental_scan() -> None:
    """Phase 4 incremental scan: byte-identical copy-forward, digest-gated incremental_kind(), and
    the full-anchor decision (forced / baseline / vep-bump / monthly)."""
    print("incremental.copy_forward / incremental_kind / scan_mode:")
    import json
    import tempfile
    from pipeline import incremental as I
    from pipeline.load import TABLE, connect

    with tempfile.TemporaryDirectory() as td:
        root = pathlib.Path(td)
        saved = {k: getattr(config, k) for k in
                 ("DUCKDB_FILE", "DUCKDB_DIR", "SNAPSHOTS_DIR", "NORMALIZED_DIR")}
        config.DUCKDB_DIR = root / "duckdb"
        config.DUCKDB_FILE = config.DUCKDB_DIR / "g.duckdb"
        config.SNAPSHOTS_DIR = root / "snap"
        config.NORMALIZED_DIR = root / "norm"
        config.NORMALIZED_DIR.mkdir(parents=True)
        (config.NORMALIZED_DIR / "dad.norm.vcf.gz").write_text("x" * 100)
        try:
            def mani(sid, clinvar, vep="112"):
                d = config.SNAPSHOTS_DIR / sid
                d.mkdir(parents=True, exist_ok=True)
                (d / "manifest.json").write_text(json.dumps(
                    {"snapshot_id": sid, "clinvar_date": clinvar, "vep_cache_version": vep}))
            mani("2026-06-01", "2026-05-30")   # prev
            mani("2026-06-02", "2026-05-30")   # same-week: ClinVar unchanged
            mani("2026-06-08", "2026-06-06")   # weekly: ClinVar changed
            mani("2026-07-01", "2026-06-06", vep="113")  # VEP cache bump

            con = connect()
            for i, sig in enumerate(["Pathogenic", "Benign"]):
                con.execute(f"INSERT INTO {TABLE} VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                            ["2026-06-01", "dad", "1", 1000 + i, "A", "T", f"G{i}",
                             "missense_variant", sig, "rev", None, "D", 0.01, None, None, "HET"])
            con.close()
            I.record_inputs("dad", "2026-06-01", I.source_digests("dad", "2026-06-01"), 2)

            check("copy_forward when ClinVar identical",
                  I.incremental_kind("dad", "2026-06-01", "2026-06-02") == "copy_forward")
            check("not copy_forward when ClinVar differs",
                  I.incremental_kind("dad", "2026-06-01", "2026-06-08") != "copy_forward")

            n = I.copy_forward("dad", "2026-06-01", "2026-06-02")
            con = connect()
            a = con.execute(f"SELECT chrom,pos,ref,alt,clinvar_sig FROM {TABLE} "
                            f"WHERE snapshot_id='2026-06-01' ORDER BY pos").fetchall()
            b = con.execute(f"SELECT chrom,pos,ref,alt,clinvar_sig FROM {TABLE} "
                            f"WHERE snapshot_id='2026-06-02' ORDER BY pos").fetchall()
            con.close()
            check("copy_forward is byte-identical", n == 2 and a == b)

            check("scan_mode forced -> full", I.scan_mode("2026-06-01", "2026-06-08", True)[0])
            check("scan_mode no-prev -> full", I.scan_mode(None, "2026-06-01", False)[0])
            check("scan_mode vep-bump -> full",
                  I.scan_mode("2026-06-08", "2026-07-01", False)[0])
            check("scan_mode no-anchor -> full",
                  I.scan_mode("2026-06-01", "2026-06-08", False)[0])
            I.mark_anchor("2026-06-01")
            full, _ = I.scan_mode("2026-06-01", "2026-06-08", False)
            check("scan_mode with monthly anchor -> incremental", not full)

            # A second scan on a day that already has a snapshot: create_snapshot() reuses
            # the id, so the predecessor must come from previous_snapshot(curr), never from
            # latest_snapshot(). prev == curr would delete-then-copy-from-itself.
            from pipeline import snapshot as S
            check("previous_snapshot is strictly older",
                  S.previous_snapshot("2026-06-08") == "2026-06-02")
            check("previous_snapshot of the baseline is None",
                  S.previous_snapshot("2026-06-01") is None)
            check("previous_snapshot never returns curr itself",
                  S.previous_snapshot("2026-07-01") == "2026-06-08")
            try:
                I.copy_forward("dad", "2026-06-02", "2026-06-02")
                same_day_refused = False
            except SystemExit:
                same_day_refused = True
            con = connect()
            still = con.execute(f"SELECT count(*) FROM {TABLE} "
                                f"WHERE snapshot_id='2026-06-02'").fetchone()[0]
            con.close()
            check("copy_forward refuses prev == curr", same_day_refused)
            check("...and leaves the snapshot's rows intact", still == 2)
        finally:
            for k, v in saved.items():
                setattr(config, k, v)


def test_incremental_clinvar_delta() -> None:
    """Variant-level incremental scan: detect ClinVar-changed coordinates, intersect with the
    variants a sample carries, and merge a re-annotated subset without touching the rest."""
    print("incremental.clinvar_delta_coords / _carried_changed / load_subset / incremental_kind:")
    import gzip
    import json
    import tempfile
    from pipeline import incremental as I
    from pipeline.load import TABLE, connect, load_subset

    with tempfile.TemporaryDirectory() as td:
        root = pathlib.Path(td)
        saved = {k: getattr(config, k) for k in
                 ("DUCKDB_FILE", "DUCKDB_DIR", "SNAPSHOTS_DIR", "NORMALIZED_DIR", "ANNOTATED_DIR")}
        config.DUCKDB_DIR = root / "duckdb"
        config.DUCKDB_FILE = config.DUCKDB_DIR / "g.duckdb"
        config.SNAPSHOTS_DIR = root / "snap"
        config.NORMALIZED_DIR = root / "norm"
        config.ANNOTATED_DIR = root / "ann"
        config.NORMALIZED_DIR.mkdir(parents=True)
        (config.NORMALIZED_DIR / "dad.norm.vcf.gz").write_text("x" * 100)
        I._DELTA_CACHE.clear()
        try:
            def clinvar(sid, records):
                d = config.SNAPSHOTS_DIR / sid
                d.mkdir(parents=True, exist_ok=True)
                (d / "manifest.json").write_text(json.dumps(
                    {"snapshot_id": sid, "clinvar_date": sid, "vep_cache_version": "112"}))
                lines = ["##fileformat=VCFv4.2",
                         "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
                for chrom, pos, ref, alt, sig in records:
                    lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\tALLELEID=1;CLNSIG={sig}")
                with gzip.open(d / "clinvar.vcf.gz", "wt") as fh:
                    fh.write("\n".join(lines) + "\n")
            # The CFTR-style deletion is the indel case: ClinVar writes it anchored
            # ("2 800 GA G"), VEP/DuckDB write it shifted and unanchored ("2 801 A -").
            clinvar("2026-06-24", [("17", 43000000, "C", "T", "Uncertain_significance"),
                                   ("2", 800, "GA", "G", "Uncertain_significance"),
                                   ("3", 1000, "A", "G", "Benign")])
            clinvar("2026-07-01", [("17", 43000000, "C", "T", "Pathogenic"),        # modified
                                   ("2", 800, "GA", "G", "Pathogenic"),             # modified indel
                                   ("17", 7670000, "G", "A", "Likely_pathogenic")])  # added; '3' removed

            con = connect()
            for chrom, pos, ref, alt, sig in [("17", 43000000, "C", "T", "Uncertain_significance"),
                                              ("17", 7670000, "G", "A", None),
                                              ("2", 801, "A", "-", "Uncertain_significance"),
                                              ("5", 9999, "T", "C", None)]:  # never in ClinVar
                con.execute(f"INSERT INTO {TABLE} VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                            ["2026-06-24", "dad", chrom, pos, ref, alt, "G", "missense_variant",
                             sig, "r", None, "D", 0.001, None, None, "HET"])
            con.close()

            coords = set(I.clinvar_delta_coords("2026-06-24", "2026-07-01"))
            check("Δ detects modified + added + removed ClinVar records",
                  coords == {("17", 43000000, "C", "T"), ("17", 7670000, "G", "A"),
                             ("2", 800, "GA", "G"), ("3", 1000, "A", "G")})
            carried = set(I._carried_changed("dad", "2026-06-24", list(coords)))
            check("carried∩Δ excludes never-in-ClinVar + non-carried",
                  ("5", 9999, "T", "C") not in carried and ("3", 1000, "A", "G") not in carried)
            check("carried∩Δ still matches SNVs exactly",
                  {("17", 43000000, "C", "T"), ("17", 7670000, "G", "A")} <= carried)
            # the regression: an anchored ClinVar deletion vs VEP's shifted dash notation
            check("carried∩Δ matches an indel across the ClinVar/VEP representation gap",
                  ("2", 801, "A", "-") in carried)
            check("a same-length Δ allele does not match a different allele at that position",
                  I.delta_match_rows([("17", 43000000, "C", "A")]) ==
                  [["17", 43000000, "C", "A"]])
            check("an indel Δ expands to the anchored and the shifted coordinate",
                  I.delta_match_rows([("2", 800, "GA", "G")]) ==
                  [["2", 800, None, None], ["2", 801, None, None]])

            # The way back: database keys -> records in the (anchored) normalized VCF.
            check("an insertion key also looks up its anchor base",
                  I.site_regions([("1", 101, "-", "A"), ("17", 43000000, "C", "T")])
                  == [("1", 100), ("1", 101), ("17", 43000000)])
            import shutil
            if shutil.which("bcftools") and shutil.which("tabix") and shutil.which("bgzip"):
                vcf = config.NORMALIZED_DIR / "kid.norm.vcf.gz"
                body = ("##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n"
                        '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tkid\n"
                        "1\t100\t.\tG\tGA\t.\t.\t.\tGT\t0/1\n"       # DB: 1:101 -/A
                        "1\t200\t.\tGAAAAA\tG\t.\t.\t.\tGT\t0/1\n"   # DB: 1:201 AAAAA/-
                        "1\t300\t.\tC\tT\t.\t.\t.\tGT\t0/1\n"
                        "1\t900\t.\tC\tT\t.\t.\t.\tGT\t0/1\n")        # not asked for
                subprocess.run(["bgzip", "-c"], input=body.encode(), check=True,
                               stdout=open(vcf, "wb"))
                subprocess.run(["tabix", "-f", "-p", "vcf", str(vcf)], check=True)
                out = I._extract_sites("kid", "2026-07-01",
                                       [("1", 101, "-", "A"), ("1", 201, "AAAAA", "-"),
                                        ("1", 300, "C", "T")])
                got = [ln.split("\t")[1] for ln in gzip.open(out, "rt")
                       if not ln.startswith("#")]
                check("_extract_sites retrieves the insertion, the deletion and the SNV, "
                      "each once", got == ["100", "200", "300"])
            else:
                print("  SKIP  _extract_sites (bcftools/tabix/bgzip not on PATH)")

            # ledger says genome+cache unchanged, only ClinVar moved -> clinvar_delta
            I.record_inputs("dad", "2026-06-24", I.source_digests("dad", "2026-06-24"), 3)
            check("incremental_kind picks clinvar_delta",
                  I.incremental_kind("dad", "2026-06-24", "2026-07-01") == "clinvar_delta")

            # merge a re-annotated subset: replace the two changed, keep the untouched one
            I.copy_forward("dad", "2026-06-24", "2026-07-01")
            tsv = config.ANNOTATED_DIR / "sub.vep.tsv.gz"
            tsv.parent.mkdir(parents=True, exist_ok=True)

            def vrow(chrom, pos, ref, alt, sig):
                return "\t".join([f"{chrom}_{pos}_{ref}/{alt}", f"{chrom}:{pos}", alt, "G",
                                  "missense_variant", "-", "-", "-", "-", "HET", sig, "crit",
                                  "D", "1"])
            with gzip.open(tsv, "wt") as fh:
                fh.write("#h\n" + vrow("17", 43000000, "C", "T", "Pathogenic") + "\n"
                         + vrow("17", 7670000, "G", "A", "Likely_pathogenic") + "\n")
            load_subset("dad", "2026-07-01", tsv)
            con = connect()
            rows = dict(con.execute(
                f"SELECT pos, clinvar_sig FROM {TABLE} WHERE snapshot_id='2026-07-01'").fetchall())
            total = con.execute(f"SELECT count(*) FROM {TABLE} "
                                f"WHERE snapshot_id='2026-07-01'").fetchone()[0]
            con.close()
            check("load_subset replaced changed rows, no duplicates",
                  total == 4 and rows[43000000] == "Pathogenic"
                  and rows[7670000] == "Likely_pathogenic" and rows[9999] is None)
        finally:
            I._DELTA_CACHE.clear()
            for k, v in saved.items():
                setattr(config, k, v)


def test_validate_outcome() -> None:
    """The retrospective harness's verdict: did a past delta hold up in later ClinVar?"""
    print("validate._outcome (confirmed / reverted / stable_noise / unknown):")
    from pipeline.validate import _outcome, _MISSING
    check("new_pathogenic that stayed pathogenic -> confirmed",
          _outcome("new_pathogenic", "Pathogenic", "Pathogenic") == "confirmed")
    check("new_pathogenic that fell back to VUS -> reverted",
          _outcome("new_pathogenic", "Pathogenic", "Uncertain_significance") == "reverted")
    check("downgraded that stayed benign -> confirmed",
          _outcome("downgraded", "Benign", "Benign") == "confirmed")
    check("af_shift with no pathogenicity move -> stable_noise",
          _outcome("af_shift", "Benign", "Benign") == "stable_noise")
    check("reclassified strengthening further -> confirmed",
          _outcome("reclassified", "Uncertain_significance", "Pathogenic") == "confirmed")
    check("variant absent at truth -> unknown",
          _outcome("new_pathogenic", "Pathogenic", _MISSING) == "unknown")


def test_validate_classify_missing_truth() -> None:
    """Regression for the caller bug: `_classify` (what `validate()` actually runs) must pass the
    _MISSING sentinel through to `_outcome` unconverted. Previously it substituted None first,
    so `_outcome`'s `is _MISSING` branch never fired and a variant absent from the truth snapshot
    fell through to the ordinal comparison — where None ordinal-ranks as 0 (benign), silently
    misclassifying a new_pathogenic/upgraded finding as 'reverted' instead of 'unknown'."""
    print("validate._classify (truth-absent variant -> unknown, excluded from actionable precision):")
    from pipeline.validate import _classify
    findings = [
        {"chrom": "1", "pos": 100, "ref": "A", "alt": "T",
         "change_kind": "new_pathogenic", "curr_sig": "Pathogenic", "tier": "actionable_now"},
    ]
    rows, counts, act = _classify("run1", findings, {})  # empty truth_sigs -> key absent
    check("absent-from-truth finding classified unknown, not reverted",
          rows[0][6] == "unknown")
    check("counts tally the unknown", counts["unknown"] == 1 and counts["reverted"] == 0)
    check("actionable_now tier excludes it from the precision denominator (act unchanged)",
          act["confirmed"] == 0 and act["reverted"] == 0)
    check("stored truth_sig is NULL (sentinel not persisted)", rows[0][7] is None)


def test_triage_model_factor_one_sided() -> None:
    """The in-silico factor may only be raised by evidence FOR pathogenicity: a confidently
    benign AlphaMissense call must not score like a confidently pathogenic one."""
    print("triage._f_model / score_delta one-sidedness:")
    from pipeline import triage as T

    benign_am, path_am = (0.02, "likely_benign"), (0.98, "likely_pathogenic")
    check("confident-benign AM contributes nothing",
          T._f_model(benign_am, None, None, "A", "G") == 0.0)
    check("confident-pathogenic AM still scores ~1.0",
          T._f_model(path_am, None, None, "A", "G") > 0.9)
    check("AM at the 0.5 boundary is 0", T._f_model((0.5, "ambiguous"), None, None, "A", "G") == 0.0)

    # The exact case that reached the actionable_now cut-off on benign evidence: a rare,
    # on-panel missense newly called Pathogenic-adjacent by nothing but a benign AM score.
    kw = dict(prev_sig="Uncertain_significance", curr_sig="Pathogenic",
              prev_rev="no_assertion_criteria_provided",
              curr_rev="no_assertion_criteria_provided",
              prev_af=1e-5, curr_af=1e-5, consequence="missense_variant",
              ref="A", alt="G", acmg_call={"tier": "VUS", "csq": "missense_variant"},
              in_panel=True, constrained=False, stability=0.0)
    v_benign = T.score_delta(am=benign_am, **kw)
    v_path = T.score_delta(am=path_am, **kw)
    check("benign AM no longer promotes to actionable_now",
          v_benign["tier"] != "actionable_now")
    check("pathogenic AM scores strictly higher", v_path["score"] > v_benign["score"])
    check("the model factor itself is 0 on the benign side",
          v_benign["factors"]["model"] == 0.0)


def test_diff_vep_cache_guard() -> None:
    """A DELTA spanning a VEP cache bump must not report cache artefacts as database change."""
    print("diff.cache_mismatch:")
    import json
    import tempfile
    from pipeline import diff as D

    with tempfile.TemporaryDirectory() as td:
        saved = config.SNAPSHOTS_DIR
        config.SNAPSHOTS_DIR = pathlib.Path(td)
        try:
            def mani(sid, vep):
                d = config.SNAPSHOTS_DIR / sid
                d.mkdir(parents=True, exist_ok=True)
                (d / "manifest.json").write_text(json.dumps(
                    {"snapshot_id": sid, "clinvar_date": "2026-06-06",
                     "vep_cache_version": vep}))
            mani("2026-06-01", "112")
            mani("2026-06-08", "112")
            mani("2026-07-01", "113")
            check("same cache version -> no mismatch",
                  D.cache_mismatch("2026-06-01", "2026-06-08") is None)
            check("cache bump is detected",
                  D.cache_mismatch("2026-06-08", "2026-07-01") == ("112", "113"))
            check("missing manifest degrades to 'cannot tell', not a crash",
                  D.cache_mismatch("2026-06-08", "1999-01-01") is None)
        finally:
            config.SNAPSHOTS_DIR = saved


def test_low_severity_verdict_fixes() -> None:
    """Review round 2, low severity: fixes that change what a report says."""
    print("panel XL compound het / HLA subtypes / triage-ACMG AM gating / SpliceAI routing:")
    from pipeline import acmg, hla_typing, panel, spliceai, triage

    # #21 — an XL female with two different pathogenic hets in one gene is a compound
    # heterozygote, not a plain carrier (the check used to be AR-only).
    two_hets = [_v("F8", ("X", 100, "A", "T"), "HET"), _v("F8", ("X", 200, "C", "G"), "HET")]
    one_het = [_v("F8", ("X", 100, "A", "T"), "HET")]
    check("XL female, two different hets -> compound (was: carrier)",
          panel._status("XL", "female", two_hets) == "compound")
    check("XL female, single het -> still carrier",
          panel._status("XL", "female", one_het) == "carrier")
    check("XL male, single het -> still affected (hemizygous)",
          panel._status("XL", "male", one_het) == "affected")
    check("AR compound het unchanged", panel._status("AR", "female", two_hets) == "compound")

    # #25 — B*27:06 / B*27:09 are the non-AS-associated subtypes; DQ2.5 also forms with
    # DQA1*05:05, which encodes the same mature alpha chain as *05:01.
    def b27_line(alleles):
        return next(c for c in hla_typing._risk_calls({"B": alleles}) if "B27" in c)

    check("B*27:05 still reported as present", "present (" in b27_line(["B*27:05:02"]))
    check("B*27:09 flagged as non-associated", "non-associated" in b27_line(["B*27:09"]))
    check("B*27:06 flagged as non-associated", "non-associated" in b27_line(["B*27:06"]))
    check("no B27 still absent", "absent" in b27_line(["B*07:02"]))

    def celiac(dqa1, dqb1):
        return next(c for c in hla_typing._risk_calls({"DQA1": dqa1, "DQB1": dqb1})
                    if "Celiac" in c)

    check("DQ2.5 via DQA1*05:01", "DQ2.5" in celiac(["DQA1*05:01"], ["DQB1*02:01"]))
    check("DQ2.5 via DQA1*05:05 (was missed)", "DQ2.5" in celiac(["DQA1*05:05"], ["DQB1*02:01"]))
    check("DQA1*05:05 without DQB1*02:01 is not DQ2.5",
          "neither" in celiac(["DQA1*05:05"], ["DQB1*03:01"]))

    # #20 — triage and acmg must gate AlphaMissense on the same consequence set.
    check("protein_altering_variant is AM-applicable in triage, as in acmg",
          triage._am_applicable("protein_altering_variant"))
    check("triage AM gate == acmg._MISSENSE",
          all(triage._am_applicable(t) for t in acmg._MISSENSE))
    check("synonymous is not AM-applicable", not triage._am_applicable("synonymous_variant"))

    # #29 — SNV/indel routing must use the normalized alleles the lookup key is built from.
    # _file_for returns None for a table that is not installed, so on a machine without the
    # SpliceAI VCFs both sides collapsed to None and the routing assertion compared
    # None with None. Point the two paths at real temp files: the routing decision is pure
    # logic and should be checked everywhere, not silently skipped where the data is missing.
    import tempfile as _tf
    _saved_sa = (config.SPLICEAI_SNV_FILE, config.SPLICEAI_INDEL_FILE)
    _d = pathlib.Path(_tf.mkdtemp())
    config.SPLICEAI_SNV_FILE = _d / "snv.vcf.gz"
    config.SPLICEAI_INDEL_FILE = _d / "indel.vcf.gz"
    config.SPLICEAI_SNV_FILE.write_bytes(b"")
    config.SPLICEAI_INDEL_FILE.write_bytes(b"")
    try:
        check("padded MNV routes as the SNV it normalizes to",
              spliceai._file_for("AT", "GT") == spliceai._file_for("A", "G"))
        check("a real indel still routes to the indel table",
              spliceai._file_for("A", "AT") != spliceai._file_for("A", "G"))
        check("an uninstalled table routes to None, not to the wrong file",
              (config.SPLICEAI_INDEL_FILE.unlink(),
               spliceai._file_for("A", "AT") is None
               and spliceai._file_for("A", "G") is not None)[1])
    finally:
        config.SPLICEAI_SNV_FILE, config.SPLICEAI_INDEL_FILE = _saved_sa
    check("the prefix cache is per file, not global", isinstance(spliceai._prefix, dict))


HFE_SIG = "Pathogenic/Pathogenic,_low_penetrance|risk_factor"


def test_clinsig_compound() -> None:
    """ClinVar's CLNSIG is a compound expression, and an exact-string list dropped every
    compound record — HFE C282Y homozygous came out of a FULL report as "0 actionable"."""
    print("clinsig: compound ClinVar classifications are matched per term:")
    import json
    import tempfile
    from pipeline import clinsig, panel as P, report as R, triage
    from pipeline.load import TABLE, connect

    yes = ["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic", HFE_SIG,
           "Pathogenic,_low_penetrance", "Likely_pathogenic,_low_penetrance",
           "Pathogenic|other", "Pathogenic|drug_response", "Pathogenic/Likely_risk_allele",
           "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance|other"]
    no = [None, "", "Benign", "Uncertain_significance", "risk_factor", "Likely_risk_allele",
          "Conflicting_classifications_of_pathogenicity",
          "Conflicting_classifications_of_pathogenicity|risk_factor",
          "Uncertain_significance|risk_factor", "Benign/Likely_benign", "Pathogenicish"]
    check("every (Likely) pathogenic compound is actionable",
          all(clinsig.is_actionable(s) for s in yes))
    check("conflicting / uncertain / bare risk terms are not",
          not any(clinsig.is_actionable(s) for s in no))
    con = connect()
    sql = [con.execute(f"SELECT {clinsig.sql_actionable('?')}", [s]).fetchone()[0]
           for s in yes + no]
    con.close()
    check("the SQL predicate agrees with the Python one on every case",
          [bool(x) for x in sql] == [clinsig.is_actionable(s) for s in yes + no])

    check("an uncertain value with an appended assertion is still uncertain",
          clinsig.is_uncertain("Uncertain_significance|risk_factor")
          and clinsig.is_uncertain(None) and not clinsig.is_uncertain(HFE_SIG)
          and not clinsig.is_uncertain("Benign"))
    check("Pathogenic -> compound pathogenic is a reclassification, not a downgrade",
          triage.change_kind("Pathogenic", HFE_SIG, "r", "r", None, None)[0] == "reclassified")
    check("VUS -> compound pathogenic is new_pathogenic",
          triage.change_kind("Uncertain_significance", HFE_SIG, "r", "r", None, None)[0]
          == "new_pathogenic")
    check("compound pathogenic -> VUS is a downgrade",
          triage.change_kind(HFE_SIG, "Uncertain_significance", "r", "r", None, None)[0]
          == "downgraded")

    check("the qualifier reaches the reader",
          clinsig.flags(HFE_SIG) == ["low penetrance (some submitters)", "risk factor"]
          and clinsig.flags("Pathogenic,_low_penetrance") == ["low penetrance"]
          and clinsig.flags("Pathogenic") == [])
    check("a table cell never contains the column separator",
          "|" not in clinsig.cell(HFE_SIG) and "low penetrance" in clinsig.cell(HFE_SIG)
          and clinsig.cell("Pathogenic") == "Pathogenic")
    check("unqualified Pathogenic sorts ahead of the compound",
          clinsig.rank("Pathogenic") < clinsig.rank(HFE_SIG) < clinsig.rank("Likely_pathogenic"))

    with tempfile.TemporaryDirectory() as td:
        root = pathlib.Path(td)
        saved = {k: getattr(config, k) for k in
                 ("DUCKDB_FILE", "DUCKDB_DIR", "SNAPSHOTS_DIR", "REPORTS_DIR")}
        config.DUCKDB_DIR = root / "duckdb"
        config.DUCKDB_FILE = config.DUCKDB_DIR / "g.duckdb"
        config.SNAPSHOTS_DIR = root / "snap"
        config.REPORTS_DIR = root / "reports"
        try:
            d = config.SNAPSHOTS_DIR / "2026-09-28"
            d.mkdir(parents=True)
            (d / "manifest.json").write_text(json.dumps(
                {"snapshot_id": "2026-09-28", "clinvar_date": "2026-09-28",
                 "vep_cache_version": "112", "gnomad_version": "4"}))
            con = connect()
            con.execute(
                f"INSERT INTO {TABLE} VALUES ('2026-09-28', 'kid', '6', 26092913, 'G', 'A', "
                f"'HFE', 'missense_variant', ?, 'criteria_provided,_multiple_submitters,"
                f"_no_conflicts', '15048', 'Hemochromatosis_type_1', 0.03, NULL, NULL, 'HOM')",
                [HFE_SIG])
            by_gene = P._actionable_by_gene(con, "2026-09-28", "kid")
            con.close()
            text = R.full_report("2026-09-28")[0].read_text()
            row = next((ln for ln in text.splitlines() if ln.startswith("| HFE")), "")
            check("FULL report counts the HFE homozygote", "**Actionable variants:** 1" in text)
            check("...in a row with the same column count as the header",
                  row.count("|") == R._HEADER.count("|"))
            check("...carrying the low-penetrance qualification", "low penetrance" in row)
            check("panel query sees it too", list(by_gene) == ["HFE"])
        finally:
            for k, v in saved.items():
                setattr(config, k, v)


def test_diff_change_gate() -> None:
    """Regression for the DELTA stage's cost model, found by running it on whole genomes.

    `score_pair` used to call `triage.model_lookups` (a `tabix` subprocess per missense/splice
    variant) and `acmg.evaluate` on EVERY joined row, then let `score_delta` discard the
    unchanged ones. Two consecutive snapshots join to millions of rows for one sample,
    of which none may have any annotation change while tens of thousands are
    missense/splice: that many subprocesses spent to render an empty report, many minutes. Same question in SQL: 0.1s.

    Two invariants keep the gate honest:
      1. change_kind is None  =>  score_delta is None. This is what licenses skipping the
         expensive work on a row; if score_delta ever started reporting a change that
         change_kind doesn't see, the gate would silently drop findings.
      2. The SQL pre-gate (sig / revstat / gnomad_af IS DISTINCT FROM) is a SUPERSET of
         change_kind, so pushing it into the query cannot drop a finding either."""
    print("diff score_pair change gate (cheap gate before per-variant model/ACMG work):")
    from pipeline import triage

    sigs = [None, "", "Benign", "Uncertain_significance", "Pathogenic", "Likely_pathogenic",
            "Conflicting_classifications_of_pathogenicity",
            "Pathogenic/Pathogenic,_low_penetrance|risk_factor"]
    revs = [None, "no_assertion_criteria_provided", "criteria_provided,_single_submitter",
            "reviewed_by_expert_panel"]
    afs = [None, 0.0, 1e-5, 0.005, 0.02]

    gate_none_but_scored, missed_by_sql = [], []
    n = 0
    for ps in sigs:
        for cs in sigs:
            for prv in revs:
                for crv in revs:
                    for paf in afs:
                        for caf in afs:
                            n += 1
                            ck = triage.change_kind(ps, cs, prv, crv, paf, caf)
                            sd = triage.score_delta(prev_sig=ps, curr_sig=cs, prev_rev=prv,
                                                    curr_rev=crv, prev_af=paf, curr_af=caf)
                            if ck is None and sd is not None:
                                gate_none_but_scored.append((ps, cs, prv, crv, paf, caf))
                            # Mirror of the SQL WHERE clause; DuckDB's IS DISTINCT FROM is
                            # null-safe inequality, which is Python's != for these scalars.
                            sql_passes = (cs != ps) or (crv != prv) or (caf != paf)
                            if ck is not None and not sql_passes:
                                missed_by_sql.append((ps, cs, prv, crv, paf, caf))

    check(f"gate covers score_delta over all {n} input combinations",
          not gate_none_but_scored)
    if gate_none_but_scored:
        print(f"      e.g. {gate_none_but_scored[0]}")
    check("SQL pre-gate is a superset of change_kind (drops no finding)", not missed_by_sql)
    if missed_by_sql:
        print(f"      e.g. {missed_by_sql[0]}")

    # A disease-name-only ClinVar edit must NOT be a delta. Consecutive ClinVar releases can
    # differ at hundreds of carried coordinates in nothing but clinvar_disease (CLNDN), and
    # 0 findings there is correct, not a missed signal. change_kind never reads the disease name, and the SQL gate must not
    # either, or every CLNDN string edit would become a triaged "change".
    check("clinvar_disease is not one of the columns the gate keys on",
          "clinvar_disease" not in _diff_gate_sql())
    for col in ("clinvar_sig", "clinvar_revstat", "gnomad_af"):
        check(f"gate keys on {col}", col in _diff_gate_sql())


def _diff_gate_sql() -> str:
    """The WHERE clause of score_pair's join, as text — asserts the gate stays in the query."""
    import inspect

    from pipeline import diff

    src = inspect.getsource(diff.score_pair)
    start = src.index("WHERE c.clinvar_sig")
    return src[start:src.index('"""', start)]


def test_validate_empty_run() -> None:
    """Regression: a validation run with zero deltas must persist, not crash.

    Found by running `validate` against whole genomes for the first time. DuckDB's
    executemany rejects an empty parameter list outright ("requires a non-empty list of
    parameter sets"), so the run died after its validation_runs row was already inserted —
    leaving a run with no findings and no report. The offline tests never hit it because
    synthetic fixtures always produced at least one finding, and on real data zero is the
    common case: it just means the replayed snapshot pair moved no variant's significance
    (a pair that shares a ClinVar release has nothing to judge).

    `diff._persist` already guarded this; `validate` did not."""
    print("validate() with zero deltas (empty executemany guard):")
    import tempfile
    from pipeline import validate as V
    from pipeline.load import TABLE, connect

    with tempfile.TemporaryDirectory() as td:
        root = pathlib.Path(td)
        saved = {k: getattr(config, k) for k in ("DUCKDB_FILE", "DUCKDB_DIR", "REPORTS_DIR")}
        config.DUCKDB_DIR = root / "duckdb"
        config.DUCKDB_FILE = config.DUCKDB_DIR / "g.duckdb"
        config.REPORTS_DIR = root / "reports"
        try:
            con = connect()
            # Two snapshots + a truth snapshot, all annotation-identical: a real delta of zero.
            for snap in ("2026-01-01", "2026-01-02", "2026-01-03"):
                con.execute(
                    f"INSERT INTO {TABLE} VALUES (?, 'dad', '1', 100, 'A', 'T', 'BRCA1', "
                    f"'missense_variant', 'Benign', 'criteria_provided,_single_submitter', "
                    f"'12345', 'Some_disease', 0.02, NULL, NULL, 'HET')", [snap])
            con.close()

            path = V.validate("2026-01-02", "2026-01-03")
            check("validate() returns a report path instead of raising", path.exists())

            con = connect()
            runs = con.execute(
                "SELECT n_deltas, n_useful, n_noise, precision_actionable FROM validation_runs "
                "WHERE run_id = '2026-01-02_vs_2026-01-03'").fetchall()
            findings = con.execute(
                "SELECT count(*) FROM validation_findings "
                "WHERE run_id = '2026-01-02_vs_2026-01-03'").fetchone()[0]
            con.close()
            check("the zero-delta run is still recorded (audit trail, not a silent gap)",
                  len(runs) == 1)
            check("recorded as 0 deltas / 0 useful / 0 noise",
                  runs and tuple(runs[0][:3]) == (0, 0, 0))
            check("precision is NULL, not 0% (nothing was judged)",
                  runs and runs[0][3] is None)
            check("no findings rows written", findings == 0)
        finally:
            for k, v in saved.items():
                setattr(config, k, v)


def test_enformer() -> None:
    """Enformer: score-table parse + bands/predicate, and interpret's ref/alt window builder
    (SNV, '-'-encoded deletion, ref-mismatch skip, contig-edge skip)."""
    print("enformer.lookup / bands + interpret._windows:")
    import tempfile
    from pipeline import config, enformer as em, interpret as it

    # --- lookup table parse + banding -----------------------------------------
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        fh.write("chrom\tpos\tref\talt\tdelta_max\tl2_center\ttop_track\ttop_desc\n")
        fh.write("3\t15521827\tC\tT\t6.97940\t3.80093\t5212\tCAGE:liver\n")
        fh.write("12\t108880081\tG\tA\t0.40000\t1.10000\t10\tDNASE:lung\n")
        scores = pathlib.Path(fh.name)
    orig_file, em._table = config.ENFORMER_SCORES_FILE, None
    orig_cal, em._calibration = em._calibration, {}     # force uncalibrated for the raw-band checks
    config.ENFORMER_SCORES_FILE = scores
    try:
        hit = em.lookup("3", 15521827, "C", "T")
        check("parses a scored row", hit is not None and abs(hit[0] - 6.9794) < 1e-6)
        check("ref/alt mismatch -> None", em.lookup("3", 15521827, "C", "G") is None)
        check("raw fallback Δ≥3 -> high band, names the assay",
              em.label(hit, "C", "T") == "Δ6.98 high (CAGE)")
        check("is_regulatory_hit True above moderate cut", em.is_regulatory_hit(hit, "C", "T"))
        low = em.lookup("12", 108880081, "G", "A")
        check("raw fallback Δ<1 -> low band", "low" in em.label(low, "G", "A"))
        check("is_regulatory_hit False below cut", not em.is_regulatory_hit(low, "G", "A"))
        check("None -> em dash", em.label(None) == "—")

        # Calibrated: percentile of Δmax vs the benign null, bands key off it (95th/99th).
        em._calibration = {"snv": [float(i) for i in range(1, 101)],     # 1..100
                           "indel": [float(i) for i in range(1, 21)]}    # 1..20
        check("calibration_available True", em.calibration_available())
        check("SNV percentile = fraction of null ≤ Δ", abs(em.percentile(50.5, "C", "T") - 50.0) < 1e-9)
        check("≥99th pctl -> high band", em.band(99.5, "C", "T") == "high")
        check("≥95th pctl -> moderate band", em.band(96.5, "C", "T") == "moderate")
        check("<95th pctl -> low band", em.band(50.5, "C", "T") == "low")
        check("indel scored against the indel null", abs(em.percentile(10.5, "CA", "-") - 50.0) < 1e-9)
        check("label shows the percentile when calibrated",
              "p99.0" in em.label((99.5, 1.0, 5, "CAGE:x"), "C", "T"))
    finally:
        config.ENFORMER_SCORES_FILE, em._table, em._calibration = orig_file, None, orig_cal
        scores.unlink(missing_ok=True)

    # --- window builder: shrink SEQ_LEN and stub the reference for a fast test ----
    REF = "AAAAACGTACGTACGTACGT"           # index 5 == 'C' (the variant's REF base below)
    orig_len, orig_half, orig_faidx = it.SEQ_LEN, it.HALF, it._faidx
    it.SEQ_LEN, it.HALF = 11, 5
    it._faidx = lambda chrom, s, e: REF[s - 1:e]   # 1-based inclusive slice
    try:
        snv = it._windows({"chrom": "1", "pos": 6, "ref": "C", "alt": "T"})
        check("SNV: both windows are SEQ_LEN",
              snv is not None and len(snv[0]) == 11 and len(snv[1]) == 11)
        check("SNV: alt window has the substituted base at centre", snv[1][5] == "T")
        deldash = it._windows({"chrom": "1", "pos": 6, "ref": "C", "alt": "-"})
        check("'-' deletion: alt window stays SEQ_LEN, no literal '-'",
              deldash is not None and len(deldash[1]) == 11 and "-" not in deldash[1])
        check("ref-mismatch -> None (no bogus score)",
              it._windows({"chrom": "1", "pos": 6, "ref": "G", "alt": "T"}) is None)
        check("contig-edge (start<1) -> None",
              it._windows({"chrom": "1", "pos": 3, "ref": "C", "alt": "T"}) is None)
    finally:
        it.SEQ_LEN, it.HALF, it._faidx = orig_len, orig_half, orig_faidx


def test_traits_nocall() -> None:
    """Regression (bug review #3): an explicit no-call in the called VCF must surface as
    indeterminate (None), not be silently reinterpreted as homozygous reference."""
    print("traits._alleles_at (no-call ≠ hom-ref):")
    import subprocess
    from pipeline import traits

    def fake_run(stdout, returncode=0, stderr=""):
        return lambda *a, **kw: type(
            "R", (), {"stdout": stdout, "returncode": returncode, "stderr": stderr})()

    orig = subprocess.run
    try:
        vcf = pathlib.Path("dummy.vcf.gz")
        subprocess.run = fake_run("160524070\tC\tT\t./.\n")
        check("explicit ./. -> None (indeterminate), not hom-ref",
              traits._alleles_at(vcf, "6", "160524070", "C") is None)
        subprocess.run = fake_run("160524070\tC\tT\t0/.\n")
        check("partial no-call 0/. -> None (can't count effect alleles)",
              traits._alleles_at(vcf, "6", "160524070", "C") is None)
        subprocess.run = fake_run("160524070\tC\tT\t0/1\n")
        check("het 0/1 -> (ref, alt)",
              traits._alleles_at(vcf, "6", "160524070", "C") == ("C", "T"))
        subprocess.run = fake_run("160524070\tC\tT\t0/2\n")
        check("allele index beyond ALT list -> None (malformed, don't crash)",
              traits._alleles_at(vcf, "6", "160524070", "C") is None)
        subprocess.run = fake_run("")
        check("absent position -> hom-ref (variant-only VCF convention)",
              traits._alleles_at(vcf, "6", "160524070", "C") == ("C", "C"))
        # A FAILED query produces the same empty stdout as an absent position, and used to
        # fall through to the hom-ref return — a broken lookup reported as a confident
        # genotype. Only a query that succeeded licenses the absent⇒hom-ref convention.
        subprocess.run = fake_run("", returncode=1, stderr="[E::idx_find] no index")
        check("failed query -> None (indeterminate), NOT hom-ref",
              traits._alleles_at(vcf, "6", "160524070", "C") is None)
    finally:
        subprocess.run = orig


def test_faidx_batching() -> None:
    """util.faidx_bases: batches `samtools faidx` calls at _FAIDX_BATCH so a large position
    list (a big PGS, or forcecall's up-to-50k-variant union) never hands a single exec() an
    argument list over the OS limit (bug review #4)."""
    print("util.faidx_bases (batched samtools faidx):")
    import subprocess
    from pipeline import util

    orig_batch = util._FAIDX_BATCH
    orig_run = subprocess.run
    # The whole samtools call is stubbed, so point REF_FASTA at a path with no .fai: with no
    # index to read, faidx_bases does no contig pre-filtering and the synthetic contig names
    # below ("0".."6") pass through. Against the machine's real GRCh38 index they would not.
    from pipeline import config
    orig_ref = config.REF_FASTA
    config.REF_FASTA = pathlib.Path("/nonexistent/ref.fa")
    calls: list[list[str]] = []

    def fake_run(cmd, **kw):
        calls.append(cmd)
        regions = cmd[3:]   # ["samtools", "faidx", ref, *regions]
        out = "".join(f">{r.split(':')[0]}:{r.split(':')[1]}\nA\n" for r in regions)
        return type("R", (), {"stdout": out, "returncode": 0, "stderr": ""})()

    util._FAIDX_BATCH = 3
    subprocess.run = fake_run
    try:
        positions = [(str(i), str(i)) for i in range(7)]   # 7 positions, batch=3 -> 3 calls
        bases = util.faidx_bases(positions)
        check("N > batch size issues more than one samtools call", len(calls) == 3)
        check("every position resolved across the batched calls", len(bases) == 7)
        check("empty input -> no call at all, empty dict",
              util.faidx_bases([]) == {} and len(calls) == 3)

        # samtools' exit code was ignored and stderr discarded, so a missing FASTA or a bad
        # contig became missing dict keys — which every caller reads as "no reference base
        # here", i.e. the absent⇒hom-ref convention applied to a tool failure.
        def fail_run(cmd, **kw):
            return type("R", (), {"stdout": "", "returncode": 1,
                                  "stderr": "[E::fai_build3_core] Failed to open"})()

        subprocess.run = fail_run
        try:
            util.faidx_bases([("1", "100")])
            hard_fail = False
        except SystemExit:
            hard_fail = True
        check("a wholly failed faidx raises instead of returning empty", hard_fail)

        def partial_run(cmd, **kw):
            # one stray unknown contig: samtools exits 1 but still returns the good regions
            good = [r for r in cmd[3:] if not r.startswith("NOPE")]
            out = "".join(f">{r.split(':')[0]}:{r.split(':')[1]}\nA\n" for r in good)
            return type("R", (), {"stdout": out, "returncode": 1,
                                  "stderr": "[faidx] Failed to fetch sequence in NOPE:5-5"})()

        subprocess.run = partial_run
        partial = util.faidx_bases([("1", "100"), ("NOPE", "5")])
        check("a partial faidx keeps the resolved positions", partial.get(("1", "100")) == "A")
        check("...and reports the unresolved one as absent, not as reference",
              ("NOPE", "5") not in partial)
    finally:
        subprocess.run = orig_run
        util._FAIDX_BATCH = orig_batch
        config.REF_FASTA = orig_ref


def test_reference_freqs_partial_cache() -> None:
    """prs._reference_freqs: a PARTIAL 1000G lookup (one contig returns zero rows, e.g. a
    transient remote hiccup) must not be cached as if complete — bug review #6. A fully
    successful lookup still gets cached as before."""
    print("prs._reference_freqs (skip cache write on a partial result):")
    import pathlib
    import subprocess
    import tempfile
    from pipeline import config, prs

    orig_pgs_dir, orig_kg_local = config.PGS_DIR, config.KG_VCF_LOCAL_DIR
    orig_run = subprocess.run
    variants = [("1", "100", "A", "G", 1.0), ("2", "200", "C", "T", 1.0)]

    def reg_chrom(cmd: list[str]) -> str:
        return "1" if "kg_reg_1." in cmd[3] else "2"

    try:
        # --- partial: chr1 returns a row, chr2 returns nothing ---
        config.PGS_DIR = pathlib.Path(tempfile.mkdtemp())
        config.KG_VCF_LOCAL_DIR = config.PGS_DIR / "no_local_mirror"  # exists()==False, no ~/genomes touch
        subprocess.run = lambda cmd, **kw: type("R", (), {"stdout":
            "1\t100\tA\tG\t0.3\t0.25\n" if reg_chrom(cmd) == "1" else "",
            "returncode": 0, "stderr": ""})()
        freqs = prs._reference_freqs("PGS_TESTPARTIAL", variants)
        check("partial result still returned in-memory (chr1 only)",
              freqs == {("1", "100"): ("A", "G", "0.3", "0.25")})
        check("partial result NOT written to the on-disk cache",
              not prs._kg_cache("PGS_TESTPARTIAL").exists())

        # --- complete: both chroms return a row ---
        config.PGS_DIR = pathlib.Path(tempfile.mkdtemp())
        config.KG_VCF_LOCAL_DIR = config.PGS_DIR / "no_local_mirror"
        subprocess.run = lambda cmd, **kw: type("R", (), {"stdout":
            "1\t100\tA\tG\t0.3\t0.25\n" if reg_chrom(cmd) == "1"
            else "2\t200\tC\tT\t0.1\t0.15\n", "returncode": 0, "stderr": ""})()
        freqs = prs._reference_freqs("PGS_TESTFULL", variants)
        check("complete result covers both chroms", len(freqs) == 2)
        check("complete result IS written to the on-disk cache",
              prs._kg_cache("PGS_TESTFULL").exists())
    finally:
        subprocess.run = orig_run
        config.PGS_DIR, config.KG_VCF_LOCAL_DIR = orig_pgs_dir, orig_kg_local


def test_short_line_guards() -> None:
    """prs._present_genotypes / traits._alleles_at: a short/malformed external-tool output
    line must be skipped, not crash the whole report with an uncaught ValueError — bug
    review #8 (matches forcecall.read_genotypes's existing `if len(f) < N: continue` house
    pattern)."""
    print("prs._present_genotypes / traits._alleles_at (malformed line -> skip, not crash):")
    import pathlib
    import subprocess
    from pipeline import prs, traits

    orig_run = subprocess.run
    try:
        # A truncated bcftools line (e.g. missing the trailing GT column) must not raise.
        subprocess.run = lambda *a, **kw: type("R", (), {  # only 4 fields
            "stdout": "1\t100\tA\tG\n", "returncode": 0, "stderr": ""})()
        present = prs._present_genotypes(pathlib.Path("dummy.vcf.gz"), [("1", "100", "A", "G", 1.0)])
        check("prs._present_genotypes: short line skipped, no crash", present == {})

        subprocess.run = lambda *a, **kw: type("R", (), {   # first short, second well-formed
            "stdout": "1\t100\tA\tG\n1\t200\tC\tT\t0/1\n", "returncode": 0, "stderr": ""})()
        present = prs._present_genotypes(pathlib.Path("dummy.vcf.gz"),
                                          [("1", "100", "A", "G", 1.0), ("1", "200", "C", "T", 1.0)])
        check("prs._present_genotypes: well-formed lines still parsed alongside a bad one",
              present == {("1", "200"): ("C", "T", "0/1")})

        subprocess.run = lambda *a, **kw: type("R", (), {   # only 2 fields
            "stdout": "160524070\tC\n", "returncode": 0, "stderr": ""})()
        check("traits._alleles_at: short line skipped -> falls through to absent-> hom-ref",
              traits._alleles_at(pathlib.Path("dummy.vcf.gz"), "6", "160524070", "C") == ("C", "C"))
    finally:
        subprocess.run = orig_run


def test_phasing_windows() -> None:
    """Regression (bug review #9): windows text keys the phased-VCF cache — merged spans,
    0-based BED, and a changed candidate set must change the key."""
    print("phasing._windows_text (cache-invalidation key):")
    from pipeline import phasing

    cands = [{"gene": "CFTR", "chrom": "7",
              "variants": [{"pos": 117559590}, {"pos": 117560100}]}]  # 510bp apart → merge
    txt = phasing._windows_text(cands)
    lines = txt.splitlines()
    check("nearby variants merge into one window", len(lines) == 1)
    chrom, s, e = lines[0].split("\t")
    check("BED is 0-based half-open (start = pos-PAD-1)",
          (chrom, int(s), int(e)) == ("7", 117559590 - 1000 - 1, 117560100 + 1000))
    cands2 = [dict(cands[0], variants=cands[0]["variants"] + [{"pos": 117600000}])]
    check("changed candidate set -> different windows text (re-phase trigger)",
          phasing._windows_text(cands2) != txt)


def test_dashboard_discovery_and_layout() -> None:
    """Regression (bug review #5): mito reports reach the dashboard. Plus the layout
    contract — the per-person page is written into the snapshot's dated directory and
    escapes a hostile id, since the id reaches the page as markup now, not as JS."""
    print("render._discover mito + build_index dated-dir layout:")
    import tempfile
    from pipeline import render

    root = pathlib.Path(tempfile.mkdtemp())
    names = ["REPORTS_DIR", "PGX_DIR", "HLA_DIR", "REPEATS_DIR", "SV_DIR", "MITO_DIR",
             "CALLABLE_DIR", "PHASE_DIR", "TRAITS_DIR", "PRS_DIR", "ANCESTRY_DIR",
             "SNAPSHOTS_DIR", "FAMILY_FILE"]
    saved = {n: getattr(config, n) for n in names}
    try:
        for n in names:
            setattr(config, n, root / n.lower())
        from pipeline import reportpaths
        mito = reportpaths.genome_report("mito", "SAMPLE001")
        mito.parent.mkdir(parents=True)
        mito.write_text("# Mito\n")
        check("a per-genome report is filed under md/genome/<person>/",
              mito.relative_to(config.REPORTS_DIR).parts[:3] == ("md", "genome", "SAMPLE001"))
        labels = [lbl for _g, lbl, _p in render._discover("SAMPLE001")]
        check("mito report discovered on the dashboard",
              any("Mitochondrial" in lbl for lbl in labels))

        snap = reportpaths.html_root() / "2026-01-02"
        out = render.build_index('Bad"Id', out_dir=snap)
        check("page written into the dated directory", out.parent == snap)
        # _esc leaves quotes alone (harmless in a report body, fatal in an href), so ids
        # and paths that reach an attribute must go through _attr instead.
        check("attribute escaping covers quotes", render._attr('a"b') == "a&quot;b")
        shell = out.read_text()
        check("the person page has a reading pane the tree's links open in",
              f'<iframe class="pane" name="{render.PANE}"' in shell
              and f'target="{render.PANE}"' in render._links(
                  [("L", snap / "x.md")], snap, target=render.PANE))
        check("the pane's starting page is carried escaped, not as live markup",
              shell.count("<main") == 0 and "&lt;main class=&quot;content&quot;&gt;" in shell)
        check("a hostile id cannot break out of an href",
              'index_Bad"Id.html' not in out.read_text())

        # An older snapshot must not present the undated per-sample reports as its own:
        # an older page was listing the PGx report of a genome sequenced after it.
        for snap in ("2026-03-04", "2026-05-06"):   # a newer one must exist for the
            d = config.SNAPSHOTS_DIR / snap          # older one to count as older
            d.mkdir(parents=True)
            (d / "manifest.json").write_text("{}")
        check("current_snapshot is the newest one", render.current_snapshot() == "2026-05-06")
        old_dir = reportpaths.html_root() / "2026-03-04"
        groups = {g for g, _l, _p in render._discover(
            "SAMPLE001", "2026-03-04", render.STALE_GROUP)}
        check("undated stage reports are re-grouped on an older snapshot",
              groups == {render.STALE_GROUP})
        check("and stay in place when no group is named",
              "Health" in {g for g, _l, _p in render._discover("SAMPLE001", "2026-03-04")})
        page = render.build_index("SAMPLE001", out_dir=old_dir, snap="2026-03-04").read_text()
        check("an older snapshot's page says which reports are not from it",
              "older" in page and render.STALE_GROUP in page)
    finally:
        for n, v in saved.items():
            setattr(config, n, v)


def test_hla_uses_forcecall() -> None:
    """`hla` must read force-called genotypes like `traits` does. forcecall._collect_sites()
    already force-calls the HLA tag SNPs, so before this the stage computed evidenced
    genotypes and discarded them: a no-call at a DQ2.5/DQ8 tag fell back to absent->hom-ref
    and the report concluded 'celiac genetically very unlikely' — a false negative."""
    print("hla.hla_report (force-called genotypes; no-call != negative):")
    import tempfile, pathlib
    from pipeline import forcecall, hla

    tmp = pathlib.Path(tempfile.mkdtemp())
    saved = {n: getattr(config, n) for n in ("HLA_DIR", "HLA_FILE", "CALLED_DIR")}
    rules = tmp / "hla_tags.tsv"
    rules.write_text("trait\trsid\tchrom\tpos\tref\teffect\tinterp0\tinterp1\tinterp2\n"
                     "Celiac DQ2.5 haplotype\trs2187668\t6\t32605884\tC\tT\tabsent\thet\thom\n")
    called = tmp / "called"; called.mkdir()
    (called / f"x.{config.BUILD}.vcf.gz").write_bytes(b"")   # exists; never read when forced
    orig_read = forcecall.read_genotypes
    try:
        config.HLA_DIR, config.HLA_FILE, config.CALLED_DIR = tmp / "out", rules, called
        # Force-call says: the reads could not resolve this site.
        forcecall.read_genotypes = lambda s: {("6", "32605884"):
                                              {"ref": "C", "alt": "T", "gt": "./.",
                                               "alleles": None, "dp": 2, "nocall": True}}
        text = hla.hla_report("x").read_text()
        check("no-call surfaces as indeterminate, not a genotype", "./." in text)
        check("celiac verdict is indeterminate, not 'very unlikely'",
              "indeterminate" in text.lower() and "very unlikely" not in text)
        check("caveat states genotypes came from the force-called VCF",
              "force-called" in text)
    finally:
        forcecall.read_genotypes = orig_read
        for n, v in saved.items():
            setattr(config, n, v)


def test_serve_extension_gate() -> None:
    """Serving is gated on the file extension as well as on the directory: the old
    octet-stream fallback served any file under an allowed dir over an unauthenticated LAN
    port. Reports now have a tree of their own, but the second gate stays."""
    print("serve.content_type (report extensions served, genotype files refused):")
    import pathlib
    from pipeline.serve import content_type

    for name in ("panel_dad.md", "index_dad.html", "style.css", "notes.txt"):
        check(f"serves {pathlib.Path(name).suffix}", content_type(pathlib.Path(name)) is not None)
    for name in ("dad.GRCh38.vcf.gz", "dad.GRCh38.vcf.gz.tbi", "dad.thinned.pgen",
                 "dad.thinned.psam", "dad.combined.bed", "scores_2026-07-01.tsv",
                 "job_2026-07-01.tsv.gz", "genomes.duckdb"):
        check(f"refuses {name}", content_type(pathlib.Path(name)) is None)


def test_report_dirs_wired() -> None:
    """A report-producing stage must be reachable from everything that needs to know about
    it. That used to mean three hand-kept lists (render's roots, the server allowlist, the
    dashboard's discovery loop) and the same omission was fixed twice. Now every report is
    filed through reportpaths, so the check is: the server exposes the HTML tree and nothing
    else, the dashboard lists every per-genome report kind, and every listed kind has a
    stage that writes it."""
    print("reportpaths wired into serve._ALLOWED + render._discover + the stages:")
    import re
    import tempfile

    from pipeline import render, reportpaths, serve

    check("the server exposes only the HTML reports tree", serve._ALLOWED == ("reports/html",))

    root = pathlib.Path(tempfile.mkdtemp())
    names = ["REPORTS_DIR", "SNAPSHOTS_DIR", "FAMILY_FILE"]
    saved = {n: getattr(config, n) for n in names}
    try:
        for n in names:
            setattr(config, n, root / n.lower())
        kinds = [k for _g, _l, k in render.GENOME_REPORTS]
        for kind in kinds:
            md = reportpaths.genome_report(kind, "SAMPLE001")
            md.parent.mkdir(parents=True, exist_ok=True)
            md.write_text("# Report\n")
        found = [p.stem for _g, _l, p in render._discover("SAMPLE001")]
        check("every per-genome report kind is discovered, in the listed order",
              found == [f"{k}_SAMPLE001" for k in kinds])
        check("an HTML path mirrors the markdown path in its own tree",
              reportpaths.html_for(reportpaths.genome_report("prs", "SAMPLE001"))
              == config.REPORTS_DIR / "html" / "genome" / "SAMPLE001" / "prs_SAMPLE001.html")

        # Every kind the dashboard lists must have a stage that files a report under it,
        # and every stage that files one must be listed — otherwise its report is invisible.
        src = "".join(p.read_text() for p in
                      (pathlib.Path(__file__).resolve().parent.parent / "pipeline").glob("*.py"))
        written = set(re.findall(r'genome_report\("(\w+)"', src))
        check(f"listed kinds == written kinds (listed-only: {sorted(set(kinds) - written)}, "
              f"written-only: {sorted(written - set(kinds))})", set(kinds) == written)

        # Snapshot-scoped reports, and the latest-only rule for md/.
        for snap in ("2026-01-01", "2026-02-01"):
            d = config.SNAPSHOTS_DIR / snap
            d.mkdir(parents=True)
            (d / "manifest.json").write_text("{}")
        old = reportpaths.snapshot_report("full", "SAMPLE001", "2026-01-01")
        new = reportpaths.snapshot_report("full", "SAMPLE001", "2026-02-01")
        check("the latest snapshot is filed under md/, an older one under md-history/",
              new.relative_to(config.REPORTS_DIR).parts[0] == "md"
              and old.relative_to(config.REPORTS_DIR).parts[0] == "md-history")
        check("...and both render to the one html tree",
              reportpaths.html_for(old).relative_to(config.REPORTS_DIR).parts[:2]
              == ("html", "2026-01-01"))
        for md in (old, new):
            md.parent.mkdir(parents=True, exist_ok=True)
            md.write_text("# FULL\n")
        interp = reportpaths.family_report(f"interpret_{config.INTERPRET_MODEL}", "2026-02-01")
        interp.parent.mkdir(parents=True, exist_ok=True)
        interp.write_text("# Interpret\n")
        labels = [lbl for _g, lbl, _p in render._discover("SAMPLE001")]
        check("interpret (current model) discovered",
              any("Regulatory effect" in lbl for lbl in labels))
        interp.rename(interp.with_name("interpret_othermodel.md"))
        check("interpret (other-model fallback glob) discovered",
              any("Regulatory effect" in lbl for _g, lbl, _p in render._discover("SAMPLE001")))

        # A newer snapshot arrives: what was latest leaves md/.
        stale = reportpaths.md_root() / "2026-02-01"
        d = config.SNAPSHOTS_DIR / "2026-03-01"
        d.mkdir(parents=True)
        (d / "manifest.json").write_text("{}")
        moved = reportpaths.archive_older("2026-03-01")
        check("archive_older moves the superseded snapshot out of md/",
              moved == ["2026-02-01"] and not stale.exists()
              and (reportpaths.history_root() / "2026-02-01" / "SAMPLE001"
                   / "full_SAMPLE001.md").exists())
        check("...and never the per-genome reports",
              reportpaths.genome_report("prs", "SAMPLE001").exists())
        check("meta reads a report's identity off its path",
              reportpaths.meta(reportpaths.genome_report("hla_typing", "SAMPLE001"))
              == {"report": "hla_typing", "person": "SAMPLE001", "callset": "SAMPLE001",
                  "scope": "genome"})
    finally:
        for n, v in saved.items():
            setattr(config, n, v)


def test_normalize_force_reindex() -> None:
    """normalize(force=True) over an existing .tbi must not destroy the VCF it just built.

    `tabix -p vcf` refuses to overwrite an existing index ("the index file exists"), which
    is a non-zero exit -> the cleanup handler deleted the freshly normalized VCF. The
    skip-guard means --force is the *only* way to reach the re-index, so every forced
    re-normalization lost the output, and the surviving stale .tbi made each later plain
    `normalize` rebuild-and-delete as well."""
    print("normalize(force=True) with an existing index (rebuild must survive):")
    import shutil
    import tempfile
    from pipeline import config, normalize as norm_mod

    if any(shutil.which(t) is None for t in ("bcftools", "tabix", "bgzip", "samtools")):
        print("  SKIP  bcftools/tabix/bgzip/samtools not on PATH")
        return

    saved = {n: getattr(config, n) for n in ("NORMALIZED_DIR", "SAMPLES_TSV", "REF_FASTA")}
    tmp = pathlib.Path(tempfile.mkdtemp())
    try:
        ref = tmp / "ref.fa"
        ref.write_text(">1\n" + "ACGT" * 25 + "\n")
        subprocess.run(["samtools", "faidx", str(ref)], check=True)

        raw = tmp / "raw.vcf"
        raw.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n"
            "chr1\t5\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
            "chr1\t9\trs334\tA\tT\t50\tPASS\t.\tGT\t0/1\n")

        manifest = tmp / "samples.tsv"
        manifest.write_text("sample_id\traw_path\n" + f"s\t{raw}\n")

        config.NORMALIZED_DIR = tmp / "normalized"
        config.SAMPLES_TSV = manifest
        config.REF_FASTA = ref

        out = norm_mod.normalized_path("s")
        tbi = out.with_suffix(".gz.tbi")

        norm_mod.normalize("s")
        check("first normalize produced a VCF + index", out.exists() and tbi.exists())

        norm_mod.normalize("s", force=True)          # must not raise
        check("forced re-normalization kept the VCF", out.exists())
        check("forced re-normalization kept the index", tbi.exists())
        body = gzip.decompress(out.read_bytes()).decode()
        check("the rebuilt VCF still carries the variant", "\t5\t" in body)
        # VEP echoes a supplied ID into Uploaded_variation, where load.py's
        # "chrom_pos_ref/alt" parse discards it — strip ID so it cannot happen.
        data = [ln.split("\t") for ln in body.splitlines() if not ln.startswith("#")]
        check("both records survived normalization", len(data) == 2)
        check("the rsID was stripped from the ID column",
              all(r[2] == "." for r in data))
    finally:
        for n, v in saved.items():
            setattr(config, n, v)
        shutil.rmtree(tmp, ignore_errors=True)


def test_load_rejects_identifier_rows() -> None:
    """load.load must fail loudly on a VEP TSV whose Uploaded_variation holds an identifier.

    VEP writes the *supplied* ID there when the input VCF's ID column is populated, so
    "rs334" arrives where "11_5227002_A/T" was expected. `_data_select`'s trailing
    `WHERE pos IS NOT NULL` used to drop those rows in silence — an rsID-annotated callset
    could lose most of its variants between annotation and the database with a clean exit
    and a plausible row count."""
    print("load._guard_parsed (identifier rows fail the load instead of vanishing):")
    import shutil
    import tempfile
    from pipeline import config, load as load_mod

    saved = {n: getattr(config, n) for n in ("DUCKDB_DIR", "DUCKDB_FILE", "ANNOTATED_DIR")}
    tmp = pathlib.Path(tempfile.mkdtemp())
    try:
        config.DUCKDB_DIR = tmp / "duckdb"
        config.DUCKDB_FILE = config.DUCKDB_DIR / "g.duckdb"
        config.ANNOTATED_DIR = tmp / "ann"

        def vrow(uploaded, loc):
            return "\t".join([uploaded, loc, "T", "HBB", "missense_variant", "-", "-",
                               "-", "-", "HET", "Pathogenic", "crit", "D", "1"])

        def write(sid, rows):
            tsv = load_mod.annotated_path(sid, "s")
            tsv.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(tsv, "wt") as fh:
                fh.write("#h\n" + "".join(r + "\n" for r in rows))
            return tsv

        write("clean", [vrow("11_5227002_A/T", "11:5227002")])
        check("a coordinate-form TSV loads", load_mod.load("s", "clean") == 1)

        write("ids", [vrow("11_5227002_A/T", "11:5227002"), vrow("rs334", "11:5227002")])
        try:
            load_mod.load("s", "ids")
            check("an identifier row aborts the load", False)
        except SystemExit as e:
            check("an identifier row aborts the load", True)
            check("...naming the offending value", "rs334" in str(e))
            check("...and the count that would have been lost", "1 of 2" in str(e))

        from pipeline.load import TABLE, connect
        con = connect()
        n = con.execute(f"SELECT count(*) FROM {TABLE} WHERE snapshot_id = 'ids'").fetchone()[0]
        con.close()
        check("nothing was written for the rejected snapshot", n == 0)
    finally:
        for n_, v in saved.items():
            setattr(config, n_, v)
        shutil.rmtree(tmp, ignore_errors=True)


def test_prs_genotype_query_failure_and_fallback() -> None:
    """prs._present_genotypes must raise on a failed bcftools query, and prs._score_one must
    treat several VCFs as a fallback chain rather than a substitution.

    (a) The query ignored its exit status, so a broken index gave empty stdout -> every
    position "absent" -> scored as homozygous reference. A failure produced a confident
    percentile instead of an error.

    (b) prs_report swapped the called VCF *out* for the force-called VCF whenever one
    existed. Markers outside the curated force-call union then fell through to "assume
    reference" even though the called VCF holds their real genotype, so adding force-called
    data could make a small score less accurate than not having it."""
    print("prs._present_genotypes (rc guard) / _score_one (VCF fallback chain):")
    import subprocess
    from pipeline import config, prs

    orig_run = subprocess.run
    # PGS_DIR points at a directory that does not exist: these functions write a temp region
    # file there and must create it, which is also what a machine that has never downloaded a
    # scoring file hits in production (config.ALL_DIRS does not include PGS_DIR).
    import tempfile
    orig_pgs = config.PGS_DIR
    config.PGS_DIR = pathlib.Path(tempfile.mkdtemp()) / "never" / "created"
    # (chrom, pos, effect_allele, other_allele, weight) — effect allele is the VCF ALT
    # here, so 1/1 is dosage 2 and 0/1 is dosage 1.
    variants = [("1", "100", "A", "G", 1.0), ("1", "200", "C", "T", 1.0)]
    forced, called = pathlib.Path("forced.vcf.gz"), pathlib.Path("called.vcf.gz")

    def stub(fn):
        def run(cmd, **kw):
            return type("R", (), {"stdout": fn(str(cmd[-1])), "returncode": 0, "stderr": ""})()
        return run

    try:
        # (a) non-zero exit -> raise, never "absent everywhere"
        subprocess.run = lambda *a, **kw: type("R", (), {
            "stdout": "", "returncode": 1,
            "stderr": "[E::idx_find_and_load] Could not retrieve index file"})()
        try:
            prs._present_genotypes(pathlib.Path("broken.vcf.gz"), variants)
            check("a failed bcftools query raises instead of returning {}", False)
        except RuntimeError as e:
            check("a failed bcftools query raises instead of returning {}", True)
            check("...and the error carries bcftools' own stderr", "index" in str(e))

        # (b) forced VCF covers pos 100 only; pos 200 must come from the called VCF, not
        # from the reference base
        subprocess.run = stub(lambda v: "1\t100\tG\tA\t1/1\n" if v == "forced.vcf.gz"
                              else "1\t100\tG\tA\t1/1\n1\t200\tT\tC\t0/1\n")
        _s, _u, _c, contrib = prs._score_one([forced, called], variants)
        # contrib keys are per scoring ROW: (chrom, pos, effect, other)
        check("a site only the called VCF covers keeps its real genotype",
              contrib[("1", "200", "C", "T")] == 1.0)
        check("a site the forced VCF covers uses the forced genotype",
              contrib[("1", "100", "A", "G")] == 2.0)

        _s, _u, _c, sub_only = prs._score_one([forced], variants)
        check("substitution alone would have dropped that site to assumed-reference",
              sub_only.get(("1", "200", "C", "T")) != 1.0)

        # a no-call in the FIRST VCF is authoritative: it must not fall through and then be
        # read as hom-ref from a variant-only callset
        subprocess.run = stub(lambda v: "1\t100\tG\tA\t./.\n" if v == "forced.vcf.gz"
                              else "1\t200\tT\tC\t0/1\n")
        _s, _u, _c, contrib = prs._score_one([forced, called], variants)
        check("a force-called no-call stays out of the score, it does not fall through",
              ("1", "100", "A", "G") not in contrib)
        check("...and PGS_DIR was created rather than crashing on the region file",
              config.PGS_DIR.exists())
    finally:
        subprocess.run = orig_run
        config.PGS_DIR = orig_pgs


def test_phasing_verdicts() -> None:
    """phasing._parse_phased / _verdict: the two ways the cis/trans call could be wrong.

    (a) The phase map is keyed by (chrom,pos) because the lookup side holds VEP-style
    coordinates. `bcftools norm -m -both` splits a 1/2 genotype into two het records at one
    POS, so both landed on the same key and last-wins gave them the same haplotype: two
    alleles on *opposite* copies were reported as "in-cis, other copy intact".

    (b) "in-cis (other intact)" needed only *some* block holding two variants on one side.
    A gene with three damaging hets — two phased together, one unplaced — got the
    reassuring verdict while the third could sit on the other copy."""
    print("phasing._parse_phased / _verdict (multiallelic collision, unplaced third hit):")
    import tempfile
    from pipeline import phasing

    def vcf(records):
        """records: (pos, ref, alt, gt, ps)"""
        lines = ["##fileformat=VCFv4.2",
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts"]
        for pos, ref, alt, gt, ps in records:
            fmt, val = ("GT:PS", f"{gt}:{ps}") if ps else ("GT", gt)
            lines.append(f"7\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\t{fmt}\t{val}")
        f = pathlib.Path(tempfile.mkstemp(suffix=".vcf")[1])
        f.write_text("\n".join(lines) + "\n")
        return f

    def gene(*positions):
        return {"gene": "CFTR", "chrom": "7",
                "variants": [{"pos": p, "ref": "A", "alt": "T", "csq": "missense",
                              "sig": "Pathogenic"} for p in positions]}

    # (a) one position, two ALTs, opposite haplotypes -> must NOT read as in-cis
    ph = phasing._parse_phased(vcf([(100, "A", "T", "1|0", "1"),
                                    (100, "A", "G", "0|1", "1")]))
    check("a conflicting multiallelic position is marked ambiguous",
          ph[("7", "100")].get("ambiguous") is True)
    v = phasing._verdict(gene(100, 100), ph)
    check("...so it is not reported as in-cis", "in-cis" not in v["verdict"])
    check("...and the report says why", v["n_ambiguous"] == 2)

    # a well-behaved pair on one haplotype in one block is still resolved
    ph = phasing._parse_phased(vcf([(100, "A", "T", "1|0", "1"),
                                    (200, "A", "T", "1|0", "1")]))
    check("two variants on one haplotype in one block are still in-cis",
          phasing._verdict(gene(100, 200), ph)["verdict"].startswith("in-cis"))

    # opposite haplotypes at distinct positions are still in-trans
    ph = phasing._parse_phased(vcf([(100, "A", "T", "1|0", "1"),
                                    (200, "A", "T", "0|1", "1")]))
    check("opposite haplotypes in one block are still in-trans",
          phasing._verdict(gene(100, 200), ph)["in_trans"] is True)

    # (b) third damaging variant unplaced -> the other copy is NOT established
    ph = phasing._parse_phased(vcf([(100, "A", "T", "1|0", "1"),
                                    (200, "A", "T", "1|0", "1")]))
    v = phasing._verdict(gene(100, 200, 300), ph)
    check("an unplaced third hit forbids the 'other intact' verdict",
          "other intact" not in v["verdict"])
    check("...and is reported as partially resolved", "partially resolved" in v["verdict"])
    check("...naming how many are unplaced", "1 unplaced" in v["verdict"])

    # two same-side variants in SEPARATE blocks: relative phase across blocks is unknown
    ph = phasing._parse_phased(vcf([(100, "A", "T", "1|0", "1"),
                                    (200, "A", "T", "1|0", "2")]))
    check("one haplotype per block is not evidence of in-cis",
          "other intact" not in phasing._verdict(gene(100, 200), ph)["verdict"])

    # an indel anchor entry must still not displace a real record at that coordinate
    ph = phasing._parse_phased(vcf([(100, "AT", "A", "1|0", "1"),
                                    (101, "A", "G", "0|1", "1")]))
    check("a real record at pos+1 wins over the indel anchor entry",
          ph[("7", "101")]["alt_hap"] == 1 and not ph[("7", "101")].get("ambiguous"))
    check("the indel is still reachable at its shifted coordinate",
          ph[("7", "100")]["alt_hap"] == 0)


def test_scan_completion_marker() -> None:
    """releases.has_new_release must not treat a half-built snapshot as a processed release.

    create_snapshot writes the manifest before annotate/load/report/diff run, so a scan that
    died partway still left a snapshot whose clinvar_date matched upstream. The weekly check
    then said "up to date" and stood down — the failed scan suppressed its own retry, and
    the partial snapshot stayed newest for reports and dashboards."""
    print("snapshot.mark_scan_complete / releases.has_new_release (unfinished scan retries):")
    import json
    import shutil
    import tempfile
    from pipeline import config, releases, snapshot

    saved = config.SNAPSHOTS_DIR
    orig_newest = releases.newest_clinvar_release
    tmp = pathlib.Path(tempfile.mkdtemp())
    try:
        config.SNAPSHOTS_DIR = tmp
        for sid, date in (("2026-08-01", "2026-08-01"), ("2026-09-05", "2026-09-05")):
            d = tmp / sid
            d.mkdir(parents=True)
            (d / "manifest.json").write_text(json.dumps(
                {"snapshot_id": sid, "clinvar_date": date, "vep_cache_version": "112"}))

        releases.newest_clinvar_release = lambda: "20260905"   # upstream == newest pinned

        check("an unmarked newest snapshot is not 'complete'",
              snapshot.scan_completed("2026-09-05") is False)
        check("a pinned-but-unfinished scan is re-run, not counted as processed",
              releases.has_new_release() is True)

        snapshot.mark_scan_complete("2026-09-05")
        check("mark_scan_complete makes it complete",
              snapshot.scan_completed("2026-09-05") is True)
        check("a completed scan at the current release stands down",
              releases.has_new_release() is False)

        releases.newest_clinvar_release = lambda: "20260912"   # upstream moved on
        check("a genuinely new release still triggers a scan",
              releases.has_new_release() is True)

        # require_complete skips a newer-but-unfinished snapshot; the default does not
        d = tmp / "2026-09-12"
        d.mkdir()
        (d / "manifest.json").write_text(json.dumps(
            {"snapshot_id": "2026-09-12", "clinvar_date": "2026-09-12"}))
        check("latest_snapshot() still returns the newest pinned id",
              snapshot.latest_snapshot() == "2026-09-12")
        check("latest_snapshot(require_complete=True) skips the unfinished one",
              snapshot.latest_snapshot(require_complete=True) == "2026-09-05")

        # no markers anywhere (a database predating them) earns exactly one catch-up scan
        for sid in ("2026-09-05",):
            (tmp / sid / ".scan_complete").unlink()
        check("a marker-less database reports its newest snapshot as unfinished",
              releases.has_new_release() is True)
        check("...and require_complete finds nothing rather than guessing",
              snapshot.latest_snapshot(require_complete=True) is None)
    finally:
        config.SNAPSHOTS_DIR = saved
        releases.newest_clinvar_release = orig_newest
        shutil.rmtree(tmp, ignore_errors=True)


def test_scan_baseline_and_skipped_samples() -> None:
    """Two ways a scan could finish "successfully" having told nobody anything.

    1. Scan B dies partway; the retry on day C took B as its baseline. Everyone B had not
       reached looked newly added, and the diff leaves new samples out — so a VUS ->
       Pathogenic change between completed scan A and C produced zero findings.
    2. A sample with no normalized VCF is skipped. The scan was still marked complete, so
       normalizing the sample afterwards never made another scan due."""
    print("scan baseline skips unfinished snapshots; skipped samples are not 'processed':")
    import json
    import tempfile
    from pipeline import diff as D, releases, snapshot as S
    from pipeline.load import TABLE, connect

    with tempfile.TemporaryDirectory() as td:
        root = pathlib.Path(td)
        saved = {k: getattr(config, k) for k in
                 ("DUCKDB_FILE", "DUCKDB_DIR", "SNAPSHOTS_DIR", "NORMALIZED_DIR",
                  "SAMPLES_TSV")}
        orig_newest = releases.newest_clinvar_release
        config.DUCKDB_DIR = root / "duckdb"
        config.DUCKDB_FILE = config.DUCKDB_DIR / "g.duckdb"
        config.SNAPSHOTS_DIR = root / "snap"
        config.NORMALIZED_DIR = root / "norm"
        config.NORMALIZED_DIR.mkdir()
        config.SAMPLES_TSV = root / "samples.tsv"
        try:
            A, B, C = "2026-09-01", "2026-09-08", "2026-09-09"
            for sid in (A, B, C):
                d = config.SNAPSHOTS_DIR / sid
                d.mkdir(parents=True)
                (d / "manifest.json").write_text(json.dumps(
                    {"snapshot_id": sid, "clinvar_date": sid, "vep_cache_version": "112"}))

            def row(snap, sample, sig):
                return [snap, sample, "1", 100, "A", "T", "BRCA1", "missense_variant", sig,
                        "criteria_provided,_single_submitter", "1", "D", 1e-5, None, None,
                        "HET"]
            con = connect()
            for r in (row(A, "dad", "Benign"), row(A, "kid", "Uncertain_significance"),
                      row(B, "dad", "Benign"),              # B died before reaching kid
                      row(C, "dad", "Benign"), row(C, "kid", "Pathogenic")):
                con.execute(f"INSERT INTO {TABLE} VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", r)
            S.mark_scan_complete(A, ["dad", "kid"])

            check("the baseline for the retry is the last COMPLETED scan, not the failed one",
                  S.previous_snapshot(C) == A)
            via_failed = D.score_pair(con, B, C)
            via_done = D.score_pair(con, S.previous_snapshot(C), C)
            con.close()
            check("(the failed scan as baseline really does hide the change)",
                  via_failed == [])
            check("...and the completed baseline delivers the reclassification",
                  [(f["sample"], f["change_kind"]) for f in via_done]
                  == [("kid", "new_pathogenic")])
            S.mark_scan_complete(B, ["dad"])
            check("once B is complete it is the baseline again", S.previous_snapshot(C) == B)

            # --- skipped samples -------------------------------------------------------
            config.SAMPLES_TSV.write_text("sample_id\ndad\nkid\n")
            releases.newest_clinvar_release = lambda: "20260909"
            (config.NORMALIZED_DIR / "dad.norm.vcf.gz").write_text("x")
            S.mark_scan_complete(C, ["dad"], ["kid"])
            check("marker records who was scanned", S.scanned_samples(C) == {"dad"})
            check("a sample that still cannot be scanned does not force a rescan",
                  releases.has_new_release() is False)
            from pipeline.normalize import normalized_path
            normalized_path("kid").write_text("x")
            check("normalizing the skipped sample makes a scan due",
                  releases.unscanned_samples(C) == ["kid"]
                  and releases.has_new_release() is True)
            S.mark_scan_complete(C, ["dad", "kid"])
            check("...and stands down once it has been scanned",
                  releases.has_new_release() is False)
            (config.SNAPSHOTS_DIR / C / ".scan_complete").write_text("2026-09-09\n")
            check("a legacy bare-date marker is complete, with samples unknown",
                  S.scan_completed(C) and S.scanned_samples(C) is None
                  and releases.has_new_release() is False)
            S.mark_scan_complete(C)
            check("a marker written without a sample list claims nothing about samples",
                  S.scanned_samples(C) is None)

            # cmd_scan with nothing normalized: fail, and leave no marker behind.
            import run as RUN
            (config.SNAPSHOTS_DIR / C / ".scan_complete").unlink()
            for f in config.NORMALIZED_DIR.iterdir():
                f.unlink()
            orig_create = S.create_snapshot
            S.create_snapshot = lambda snapshot_id=None: C
            try:
                import argparse
                rc = RUN.cmd_scan(argparse.Namespace(full=False))
            finally:
                S.create_snapshot = orig_create
            check("a scan that processed no sample fails", rc == 1)
            check("...and is not marked complete", not S.scan_completed(C))
        finally:
            releases.newest_clinvar_release = orig_newest
            for k, v in saved.items():
                setattr(config, k, v)


def test_gitignore_covers_data_root() -> None:
    """.gitignore is the safety net for a GENOMES_ROOT pointed inside the checkout, and it
    was written when the pipeline had a dozen stages. `family/` (names, relationships,
    consent decisions) and the newer stage directories were never added, so `git add .`
    would have staged them. Derived from config.py so a new stage cannot be forgotten."""
    print(".gitignore covers every path config.py places under GENOMES_ROOT:")
    import shutil
    repo = pathlib.Path(__file__).resolve().parent.parent
    if not shutil.which("git") or not (repo / ".git").exists():
        print("  SKIP  (not a git checkout)")
        return
    probes = set()
    for name in dir(config):
        value = getattr(config, name)
        if not isinstance(value, pathlib.Path) or value == config.GENOMES_ROOT:
            continue
        try:
            rel = value.relative_to(config.GENOMES_ROOT)
        except ValueError:
            continue
        # A directory is probed through a file inside it; a file through itself.
        probes.add(str(rel) if rel.suffix else str(rel / "probe.json"))
    probes |= {"family/family.tsv", "family/consent.tsv", "traits/kid/traits_kid.json"}
    out = subprocess.run(["git", "check-ignore", "--no-index", "--stdin"], cwd=repo,
                         input="\n".join(sorted(probes)) + "\n", capture_output=True,
                         text=True).stdout.split("\n")
    missed = sorted(probes - set(out))
    check(f"all {len(probes)} data paths are ignored" + (f" (not: {missed})" if missed else ""),
          not missed and len(probes) > 30)
    tracked = subprocess.run(["git", "check-ignore", "--no-index", "traits/traits.tsv",
                              "tests/fixtures/trio/x.vcf"], cwd=repo,
                             capture_output=True, text=True).stdout.split()
    check("the tracked trait panel and the test fixtures stay unignored", tracked == [])


def test_annotate_retry_safety() -> None:
    """annotate must not reuse a truncated or stale TSV, and _run_vep must not leave one.

    VEP wrote straight to the final filename and the skip check tested only existence, so a
    run killed mid-VEP left a partial gzip under exactly the name the next attempt skipped to
    and loaded. Retries were rare while a failed scan stayed stopped; they are the normal
    path now that an unfinished scan re-runs itself, so the corpse must not be reusable.
    Existence is also not identity: a re-normalized sample against the same snapshot id kept
    the previous run's annotations, which is what made `scan --full` not re-annotate."""
    print("annotate._run_vep (atomic) / annotate (input stamp):")
    import shutil
    import tempfile
    from pipeline import annotate as ann
    from pipeline import config

    saved = {n: getattr(config, n) for n in
             ("ANNOTATED_DIR", "NORMALIZED_DIR", "SNAPSHOTS_DIR", "VEP_BIN", "RUNS_TSV")}
    tmp = pathlib.Path(tempfile.mkdtemp())
    calls = []
    orig_run_vep = ann._run_vep
    try:
        config.ANNOTATED_DIR = tmp / "ann"
        config.NORMALIZED_DIR = tmp / "norm"
        config.SNAPSHOTS_DIR = tmp / "snap"
        config.RUNS_TSV = tmp / "runs.tsv"
        config.VEP_BIN = tmp / "vep"
        config.VEP_BIN.write_text("#!/bin/sh\n")

        norm = config.NORMALIZED_DIR / "s.norm.vcf.gz"
        norm.parent.mkdir(parents=True)
        norm.write_text("x" * 50)
        snap = config.SNAPSHOTS_DIR / "2026-09-01"
        snap.mkdir(parents=True)
        (snap / "clinvar.vcf.gz").write_text("cv")
        (snap / "manifest.json").write_text(
            '{"snapshot_id": "2026-09-01", "clinvar_date": "2026-09-01", '
            '"vep_cache_version": "112"}')

        def fake_vep(inp, clinvar, out, cache_version):
            calls.append((out, cache_version))
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text("annotated")
        ann._run_vep = fake_vep

        out = ann.annotate("s", "2026-09-01")
        check("first annotate runs VEP", len(calls) == 1)
        check("...against the version the snapshot pins, not whatever config says",
              calls[0][1] == "112")
        check("...and writes an input stamp", ann._stamp_path(out).exists())

        ann.annotate("s", "2026-09-01")
        check("unchanged inputs skip VEP", len(calls) == 1)

        # Identity is content. A touched-but-identical VCF is the same input...
        import json
        import os
        os.utime(norm, (0, 0))
        ann.annotate("s", "2026-09-01")
        check("an mtime change alone does not re-annotate", len(calls) == 1)
        # ...and different bytes at the same size and mtime (rsync -a, cp -p, a restore) are
        # a different one. The (size, mtime) stamp called these equal and reused the TSV
        # while the scan ledger, which hashes, recorded the new input beside it.
        norm.write_text("y" * 50)
        os.utime(norm, (0, 0))
        ann.annotate("s", "2026-09-01")
        check("a changed normalized VCF re-annotates despite identical size and mtime",
              len(calls) == 2)

        # a stamp from before hashing is honoured once on its own terms, then upgraded
        st = norm.stat()
        legacy = {**json.loads(ann._stamp_path(out).read_text()),
                  "normalized_vcf": f"{st.st_size}:{int(st.st_mtime)}"}
        ann._stamp_path(out).write_text(json.dumps(legacy))
        ann.annotate("s", "2026-09-01")
        check("a matching pre-hash stamp is reused, not mass re-annotated", len(calls) == 2)
        check("...and rewritten by content",
              json.loads(ann._stamp_path(out).read_text())["normalized_vcf"]
              .startswith("sha256:"))

        # runs.tsv records the cache VEP ran with, which for an older snapshot is the
        # pinned one rather than today's config.
        from pipeline.util import read_tsv
        saved_cache = config.VEP_CACHE_VERSION, config.VEP_CACHE_DIR
        config.VEP_CACHE_DIR = tmp / "cache"
        (config.VEP_CACHE_DIR / "homo_sapiens" / f"112_{config.BUILD}").mkdir(parents=True)
        config.VEP_CACHE_VERSION = "113"
        try:
            ann.annotate("s", "2026-09-01", force=True)
        finally:
            config.VEP_CACHE_VERSION, config.VEP_CACHE_DIR = saved_cache
        check("a re-annotation after a cache bump still runs the pinned cache",
              calls[-1][1] == "112")
        check("...and the run ledger records that version, not the configured one",
              read_tsv(config.RUNS_TSV)[-1]["vep_cache_version"] == "112")
        calls.pop()

        # a TSV with no stamp (written before stamping existed) is reused, not silently redone
        ann._stamp_path(out).unlink()
        ann.annotate("s", "2026-09-01")
        check("an unstamped TSV is reused rather than forcing a mass re-annotation",
              len(calls) == 2)
        check("...but --force still redoes it",
              (ann.annotate("s", "2026-09-01", force=True), len(calls) == 3)[1])
    finally:
        ann._run_vep = orig_run_vep
        for n, v in saved.items():
            setattr(config, n, v)
        shutil.rmtree(tmp, ignore_errors=True)


def test_run_vep_atomicity() -> None:
    """_run_vep must leave no file under the final name when VEP fails."""
    print("annotate._run_vep (a failed VEP leaves no loadable output):")
    import shutil
    import subprocess
    import tempfile
    from pipeline import annotate as ann

    tmp = pathlib.Path(tempfile.mkdtemp())
    orig_run = ann.run
    try:
        out = tmp / "s.vep.tsv.gz"

        def dying_vep(cmd, **kw):
            # VEP writes some output, then dies — exactly the truncated-gzip case
            part = pathlib.Path(cmd[cmd.index("-o") + 1])
            part.write_text("half a gzip stream")
            raise subprocess.CalledProcessError(1, cmd)
        ann.run = dying_vep
        try:
            ann._run_vep(tmp / "in.vcf", tmp / "cv.vcf.gz", out, "112")
        except subprocess.CalledProcessError:
            pass
        check("a killed VEP leaves nothing under the final name", not out.exists())
        check("...and cleans up its .part file",
              not list(tmp.glob("*.part")))
    finally:
        ann.run = orig_run
        shutil.rmtree(tmp, ignore_errors=True)


def test_replacement_atomicity() -> None:
    """A DELETE-then-INSERT replacement must not leave the target emptied when the INSERT fails.

    DuckDB auto-commits per statement. `load` wrapped its own replacement by hand, but
    copy_forward and the ledger writers did not — so a scan interrupted between the two
    statements wiped that sample's rows for the snapshot rather than leaving them alone. An
    interrupted scan is exactly when it happens, and the retry then copies forward from a
    snapshot that has since been emptied."""
    print("load.transaction / incremental.copy_forward (atomic replacement):")
    import tempfile
    from pipeline import config
    from pipeline import incremental as I
    from pipeline.load import TABLE, connect

    class FailingCon:
        """Passes SQL through until it sees the statement under test, then fails."""
        def __init__(self, con, fail_on):
            self._con, self._fail_on = con, fail_on
        def execute(self, sql, *a, **kw):
            if self._fail_on in sql:
                raise RuntimeError("simulated failure mid-replacement")
            return self._con.execute(sql, *a, **kw)
        def __getattr__(self, name):
            return getattr(self._con, name)

    with tempfile.TemporaryDirectory() as td:
        root = pathlib.Path(td)
        saved = {k: getattr(config, k) for k in ("DUCKDB_FILE", "DUCKDB_DIR")}
        config.DUCKDB_DIR = root / "duckdb"
        config.DUCKDB_FILE = config.DUCKDB_DIR / "g.duckdb"
        orig_connect = I.connect
        try:
            con = connect()
            for snap, pos in (("2026-06-01", 100), ("2026-07-01", 999)):
                con.execute(f"INSERT INTO {TABLE} VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                            [snap, "dad", "1", pos, "A", "T", "G", "missense_variant",
                             "Pathogenic", "r", None, "D", 0.001, None, None, "HET"])
            con.close()

            def before():
                c = connect()
                n = c.execute(f"SELECT count(*) FROM {TABLE} WHERE snapshot_id='2026-07-01'"
                              ).fetchone()[0]
                pos = c.execute(f"SELECT pos FROM {TABLE} WHERE snapshot_id='2026-07-01'"
                                ).fetchone()
                c.close()
                return n, (pos[0] if pos else None)

            check("curr starts with its own row", before() == (1, 999))

            I.connect = lambda: FailingCon(orig_connect(), "INSERT INTO variants SELECT")
            try:
                I.copy_forward("dad", "2026-06-01", "2026-07-01")
                check("a failed copy_forward propagates the error", False)
            except RuntimeError:
                check("a failed copy_forward propagates the error", True)
            I.connect = orig_connect

            check("...and rolls back, leaving curr's rows intact", before() == (1, 999))

            # the successful path still replaces
            I.copy_forward("dad", "2026-06-01", "2026-07-01")
            check("a successful copy_forward re-stamps prev's rows onto curr",
                  before() == (1, 100))
        finally:
            I.connect = orig_connect
            for k, v in saved.items():
                setattr(config, k, v)


def test_disease_status_per_gene() -> None:
    """panel.disease_status: a multi-gene disease must not propagate one gene's status to another.

    Status used to be computed over the disease's POOLED hits, so `_status` saw has_hom set by
    a homozygote in one gene and returned "affected" — which was then written into the AR
    carrier index for every gene with a hit. A plain heterozygote in a second gene was
    recorded as affected, doubling its per-child transmission probability (0.5 -> 1.0) in the
    family reproductive-risk table."""
    print("panel.disease_status (per-gene status, no cross-gene promotion):")
    from pipeline.panel import (_AFFECTED, _ATRISK, _CARRIER, _COMPOUND, _Disease,
                                disease_status)

    # one AR disease, two genes: hom in the first, single het in the second
    dz = _Disease("Retinitis pigmentosa", ["USH2A", "RPGR"], "AR", "")
    by_gene = {
        "USH2A": [_v("USH2A", ("1", 1, "A", "T"), "HOM")],
        "RPGR": [_v("RPGR", ("X", 5, "C", "G"), "HET")],
    }
    status, per_gene = disease_status(dz, "female", by_gene)
    check("the disease row still shows the most severe status", status == _AFFECTED)
    check("the homozygous gene is affected", per_gene["USH2A"] == _AFFECTED)
    check("the heterozygous gene stays a carrier, not promoted to affected",
          per_gene["RPGR"] == _CARRIER)

    # compound het in one gene must not make a single het in another compound either
    by_gene = {
        "USH2A": [_v("USH2A", ("1", 1, "A", "T"), "HET"),
                  _v("USH2A", ("1", 2, "C", "G"), "HET")],
        "RPGR": [_v("RPGR", ("X", 5, "C", "G"), "HET")],
    }
    status, per_gene = disease_status(dz, "female", by_gene)
    check("the two-allele gene is a possible compound het", per_gene["USH2A"] == _COMPOUND)
    check("the one-allele gene is still just a carrier", per_gene["RPGR"] == _CARRIER)
    check("the row takes the more severe of the two", status == _COMPOUND)

    # single hets in two different genes are two carriers, not a compound het
    by_gene = {
        "USH2A": [_v("USH2A", ("1", 1, "A", "T"), "HET")],
        "RPGR": [_v("RPGR", ("X", 5, "C", "G"), "HET")],
    }
    status, per_gene = disease_status(dz, "female", by_gene)
    check("hets in two genes are two carriers, not a compound het",
          status == _CARRIER and set(per_gene.values()) == {_CARRIER})

    # a gene with no hits contributes nothing; a disease with no hits at all yields None
    check("a disease with no hits yields no status",
          disease_status(dz, "female", {}) == (None, {}))

    # single-gene diseases are unaffected by the change
    ad = _Disease("HBOC", ["BRCA1"], "AD", "")
    status, per_gene = disease_status(
        ad, "female", {"BRCA1": [_v("BRCA1", ("17", 1, "A", "T"), "HET")]})
    check("a single-gene AD het is still at-risk",
          status == _ATRISK and per_gene == {"BRCA1": _ATRISK})


def test_launchd_plist_has_tool_path() -> None:
    """The scheduled job must declare a PATH that actually contains the scan's tools.

    launchd starts jobs with PATH=/usr/bin:/bin:/usr/sbin:/sbin and does not read a shell
    profile, so a plist without EnvironmentVariables dies at the scan's first require_tools()
    with "Missing required tool(s): tabix, bgzip" — the weekly job would look installed and
    never produce a snapshot."""
    print("launchd plist (PATH covers the tools the scan shells out to):")
    import plistlib
    import shutil

    root = pathlib.Path(__file__).resolve().parent.parent
    plist = root / "launchd" / "com.example.genomeledger.plist"
    with plist.open("rb") as fh:
        cfg = plistlib.load(fh)

    path = cfg.get("EnvironmentVariables", {}).get("PATH", "")
    check("the plist declares a PATH", bool(path))
    dirs = path.split(":")
    check("...that still contains the system directories launchd would have given it",
          all(d in dirs for d in ("/usr/bin", "/bin")))

    # the tools cmd_scan reaches for, checked against where they actually are on this machine
    for tool in ("tabix", "bgzip", "bcftools"):
        found = shutil.which(tool)
        if found is None:
            print(f"  SKIP  {tool} not installed here")
            continue
        check(f"...and the directory holding {tool} ({str(pathlib.Path(found).parent)})",
              str(pathlib.Path(found).parent) in dirs)


def test_ingest_build_and_sample_guards() -> None:
    """ingest._looks_like_build / _sample_names: a contig length must veto a build name.

    The build check was a flat OR over the whole header, so any mention of "GRCh38" anywhere
    — a ##bcftools_viewCommand quoting a path, a ##source naming a pipeline — outvoted a
    contig length that said GRCh37 outright. And nothing checked the sample count: `annotate`
    runs VEP with --individual all, which emits a row per sample, and `load` files every row
    under the one sample_id ingest recorded, so a multisample VCF attributes everyone's
    genotypes to one person."""
    print("ingest._looks_like_build (contig length vetoes) / _sample_names:")
    from pipeline import config, ingest

    saved = config.BUILD
    try:
        config.BUILD = "GRCh38"

        def hdr(*lines):
            return list(lines) + ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1"]

        ok, _ = ingest._looks_like_build(hdr("##contig=<ID=1,length=248956422>"))
        check("a GRCh38 contig length passes", ok)
        ok, _ = ingest._looks_like_build(hdr("##contig=<ID=chr1,length=248956422>"))
        check("...with or without the chr prefix", ok)

        # the regression: GRCh37 lengths plus an incidental "GRCh38" mention
        ok, why = ingest._looks_like_build(hdr(
            "##contig=<ID=1,length=249250621>",
            "##bcftools_viewCommand=view -R /refs/GRCh38/panel.bed in.vcf"))
        check("a GRCh37 contig length is refused despite a GRCh38 mention", not ok)
        check("...and the reason names the actual assembly", "GRCh37/hg19" in why)

        ok, why = ingest._looks_like_build(hdr("##contig=<ID=chr1,length=248387328>"))
        check("a native T2T/CHM13 callset is refused for GRCh38", not ok)
        check("...and is identified as CHM13", "CHM13" in why)

        # ID=10/ID=11 and alt contigs must not be mistaken for chromosome 1
        check("ID=10 is not chromosome 1",
              ingest._chr1_length(["##contig=<ID=10,length=249250621>"]) is None)
        check("an alt contig starting with chr1 is not chromosome 1",
              ingest._chr1_length(
                  ["##contig=<ID=chr1_KI270706v1_random,length=175055>"]) is None)

        # no contig lines at all -> fall back to the name
        ok, _ = ingest._looks_like_build(hdr("##reference=GRCh38_full_analysis_set.fa"))
        check("with no contig lines the build name is still accepted", ok)
        ok, why = ingest._looks_like_build(hdr("##source=someCaller"))
        check("...and a header with neither is refused", not ok and "no build name" in why)

        # sample columns
        chrom = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        check("one sample column is counted", ingest._sample_names([chrom + "\ts1"]) == ["s1"])
        check("a multisample VCF is counted as such",
              ingest._sample_names([chrom + "\ts1\ts2\ts3"]) == ["s1", "s2", "s3"])
        check("a sites-only VCF has no sample columns",
              ingest._sample_names(["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]) == [])
    finally:
        config.BUILD = saved


def test_prs_partial_nocall_and_duplicate_coords() -> None:
    """prs._score_one: a half-called genotype must not become a homozygote, and two scoring
    rows at one coordinate must not overwrite each other's contribution.

    "1/." was reduced to its one called allele and then doubled by the haploid rule, scoring
    two effect alleles — indistinguishable from a confirmed "1/1", and always in the
    risk-raising direction. Separately, `contrib` was keyed by (chrom,pos), so a scoring file
    with two rows at one position kept only the last while `total` summed both, leaving the
    headline raw score and the percentile's recomputed score disagreeing."""
    print("prs._score_one (partial no-calls, duplicate coordinates):")
    import subprocess
    import tempfile
    from pipeline import config, prs

    orig_run, orig_pgs = subprocess.run, config.PGS_DIR
    config.PGS_DIR = pathlib.Path(tempfile.mkdtemp())
    vcf = pathlib.Path("s.vcf.gz")

    def stub(out):
        def run(cmd, **kw):
            return type("R", (), {"stdout": out, "returncode": 0, "stderr": ""})()
        return run

    try:
        v = [("1", "100", "A", "G", 1.0)]      # effect allele A == the VCF ALT

        subprocess.run = stub("1\t100\tG\tA\t1/1\n")
        _s, _u, _c, contrib = prs._score_one([vcf], v)
        check("a confirmed 1/1 scores two effect alleles",
              contrib[("1", "100", "A", "G")] == 2.0)

        subprocess.run = stub("1\t100\tG\tA\t0/1\n")
        _s, _u, _c, contrib = prs._score_one([vcf], v)
        check("a het scores one", contrib[("1", "100", "A", "G")] == 1.0)

        # the regression: half-called must be dropped, not doubled
        subprocess.run = stub("1\t100\tG\tA\t1/.\n")
        score, n_used, _c, contrib = prs._score_one([vcf], v)
        check("a partial no-call '1/.' is NOT scored as a homozygote",
              contrib.get(("1", "100", "A", "G")) != 2.0)
        check("...it is dropped from the score entirely", not contrib and score == 0.0)
        check("...and not counted as a scored site", n_used == 0)
        subprocess.run = stub("1\t100\tG\tA\t./1\n")
        _s, _u, _c, contrib = prs._score_one([vcf], v)
        check("...in either order ('./1')", not contrib)

        # a genuinely haploid call (one GT field) still doubles
        subprocess.run = stub("1\t100\tG\tA\t1\n")
        _s, _u, _c, contrib = prs._score_one([vcf], v)
        check("a haploid '1' still scales to two copies",
              contrib[("1", "100", "A", "G")] == 2.0)

        # two scoring rows at one coordinate must both survive
        dup = [("1", "100", "A", "G", 1.0), ("1", "100", "T", "G", 2.0)]
        subprocess.run = stub("1\t100\tG\tA,T\t1/2\n")
        score, _u, _c, contrib = prs._score_one([vcf], dup)
        check("both scoring rows at one coordinate keep their own contribution",
              len(contrib) == 2)
        check("...and the raw score is their sum, matching contrib",
              abs(score - sum(contrib.values())) < 1e-9 and abs(score - 3.0) < 1e-9)
    finally:
        subprocess.run = orig_run
        config.PGS_DIR = orig_pgs


def test_reference_freqs_failed_query_not_cached() -> None:
    """prs._reference_freqs must not cache a query that reported an error, even when that
    query still emitted rows.

    bcftools can stream rows and then fail — a truncated remote read, a timeout partway. The
    existing guard only asked whether a chromosome produced ZERO rows, so a non-zero exit
    that still produced output looked complete and got written to the per-PGS cache, which
    every later run then read as authoritative."""
    print("prs._reference_freqs (a failed-but-partial query is never cached):")
    import subprocess
    import tempfile
    from pipeline import config, prs

    orig_run = subprocess.run
    orig_pgs, orig_kg = config.PGS_DIR, config.KG_VCF_LOCAL_DIR
    variants = [("1", "100", "A", "G", 1.0), ("2", "200", "C", "T", 1.0)]
    try:
        config.PGS_DIR = pathlib.Path(tempfile.mkdtemp())
        config.KG_VCF_LOCAL_DIR = config.PGS_DIR / "no_local_mirror"

        # both chromosomes return rows, but chr2's query exits non-zero
        def partial_fail(cmd, **kw):
            c2 = "kg_reg_2." in str(cmd[3])
            return type("R", (), {
                "stdout": "2\t200\tC\tT\t0.1\t0.15\n" if c2 else "1\t100\tA\tG\t0.3\t0.25\n",
                "returncode": 1 if c2 else 0,
                "stderr": "[E::bcf_sr_next_line] Error reading the remote index" if c2 else ""})()
        subprocess.run = partial_fail
        freqs = prs._reference_freqs("PGS_TESTFAILROWS", variants)
        check("rows from the failed query are still usable in-memory this run",
              len(freqs) == 2)
        check("...but the result is NOT cached",
              not prs._kg_cache("PGS_TESTFAILROWS").exists())

        # a clean run of the same shape still caches
        config.PGS_DIR = pathlib.Path(tempfile.mkdtemp())
        config.KG_VCF_LOCAL_DIR = config.PGS_DIR / "no_local_mirror"
        subprocess.run = lambda cmd, **kw: type("R", (), {
            "stdout": ("2\t200\tC\tT\t0.1\t0.15\n" if "kg_reg_2." in str(cmd[3])
                       else "1\t100\tA\tG\t0.3\t0.25\n"),
            "returncode": 0, "stderr": ""})()
        freqs = prs._reference_freqs("PGS_TESTCLEAN", variants)
        check("a clean complete run is still cached",
              len(freqs) == 2 and prs._kg_cache("PGS_TESTCLEAN").exists())
    finally:
        subprocess.run = orig_run
        config.PGS_DIR, config.KG_VCF_LOCAL_DIR = orig_pgs, orig_kg


def test_ref_effect_supplement_sources() -> None:
    """prs._ref_alleles / _ref_effect_supplement: the hom-ref supplement must not depend on
    whether the 1000G reference run succeeded, and must not break its symmetry when it did.

    A variant-only VCF omits sites where the sample is hom-ref; where the effect allele IS
    the reference, the sample carries two copies plink2 cannot see. That correction used to
    be applied only inside `if ref:`, so a failed reference lookup silently left the score
    unsupplemented and the same genome reported a different raw number depending on the
    network. But it cannot simply resolve everything either: when the pvars DO exist,
    positions they omit were not scored on the reference side, so crediting the sample for
    them would inflate the percentile."""
    print("prs._ref_alleles (pvar when present, FASTA when there is no panel at all):")
    import shutil
    import tempfile
    from pipeline import config, prs

    tmp = pathlib.Path(tempfile.mkdtemp())
    saved_ref = config.REF_FASTA
    try:
        variants = [("1", "100", "A", "G", 1.0), ("1", "200", "C", "T", 2.0)]

        # --- pvars present: they are authoritative, and what they omit stays unresolved ---
        wd = tmp / "with_pvar"
        wd.mkdir()
        (wd / "chr1.pvar").write_text(
            "#CHROM\tPOS\tID\tREF\tALT\n1\t100\t.\tA\tG\n")
        refs = prs._ref_alleles(variants, wd)
        check("pvar supplies the REF it carries", refs.get(("1", "100")) == "A")
        check("a position the pvar omits stays unresolved (keeps both sides symmetric)",
              ("1", "200") not in refs)
        supp, n = prs._ref_effect_supplement(variants, set(), wd)
        check("only the pvar-backed effect==REF site contributes", (supp, n) == (2.0, 1))

        # --- a genome-wide-sized set with no pvars must NOT be resolved from the FASTA ---
        # The supplement exists to make the sample comparable with the reference panel. With
        # no panel there is nothing to compare against, and resolving REF for millions of
        # positions cost 29.3M faidx lookups per sample (3.5 h for one) to adjust a number
        # reported with no percentile at all.
        seen = []
        orig_fb = prs.faidx_bases
        prs.faidx_bases = lambda pos: (seen.append(len(list(pos))), {})[1]
        try:
            big = [("1", str(i), "A", "G", 1.0) for i in range(config.PRS_PLINK_MIN + 5)]
            check("a genome-wide set with no pvars skips the FASTA pass entirely",
                  prs._ref_alleles(big, tmp / "no_pvar_here") == {} and not seen)
            small = [("1", str(i), "A", "G", 1.0) for i in range(50)]
            prs._ref_alleles(small, tmp / "no_pvar_here")
            check("...while a small score still resolves from the FASTA", seen == [50])
        finally:
            prs.faidx_bases = orig_fb

        # --- no pvars at all: fall back to the local FASTA so the raw score is deterministic
        fa = tmp / "ref.fa"
        fa.write_text(">1\n" + "A" * 120 + "\n")
        if shutil.which("samtools"):
            import subprocess
            subprocess.run(["samtools", "faidx", str(fa)], check=True)
            config.REF_FASTA = fa
            refs = prs._ref_alleles(variants, tmp / "no_pvar_here")
            check("with no pvars the FASTA resolves REF", refs.get(("1", "100")) == "A")
            supp2, n2 = prs._ref_effect_supplement(variants, set(), tmp / "no_pvar_here")
            check("...so the supplement is computed rather than silently zero", n2 >= 1)
            check("...and a covered site is still excluded",
                  prs._ref_effect_supplement(variants, {("1", "100")},
                                             tmp / "no_pvar_here")[1] == 0)
        else:
            print("  SKIP  samtools not on PATH")
    finally:
        config.REF_FASTA = saved_ref
        shutil.rmtree(tmp, ignore_errors=True)


def test_forcecall_sites_follow_the_scored_panel() -> None:
    """forcecall._prs_sites must draw from the same panel prs_report scores.

    It read config.PRS_SCORES — the 3-score fallback list — while the report had moved to the
    44-score COMPLEX_PANEL. The two share no PGS ids, so the force-called VCF covered
    essentially none of the markers actually scored (a fraction of a percent), and while prs_report was substituting that VCF for the sample's own callset every
    small score came out as an assumed-reference genome."""
    print("forcecall._prs_sites (sources the scored panel, not the fallback list):")
    from pipeline import forcecall, prs

    orig_panel, orig_load, orig_file = (prs._panel_scores, prs._load_variants,
                                        prs._scoring_file)
    seen = []
    try:
        prs._panel_scores = lambda: [{"pgs_id": "PGS_PANEL_A"}, {"pgs_id": "PGS_PANEL_B"}]
        prs._scoring_file = lambda pid: pathlib.Path(f"/nonexistent/{pid}.txt.gz")
        def fake_load(path):
            pid = path.name.split(".")[0]
            seen.append(pid)
            return [("1", "100", "A", "G", 1.0)]
        prs._load_variants = fake_load

        sites = forcecall._prs_sites()
        check("every panel score is consulted", seen == ["PGS_PANEL_A", "PGS_PANEL_B"])
        check("...and none of the config.PRS_SCORES fallback ids are",
              not any(p in seen for p, _t in __import__("pipeline.config",
                      fromlist=["x"]).PRS_SCORES))
        check("their markers land in the site list", len(sites) == 2)

        # a genome-wide score is still left to plink2, by the same threshold prs uses
        from pipeline import config
        prs._load_variants = lambda path: [("1", str(i), "A", "G", 1.0)
                                           for i in range(config.FC_MAX_LIST + 2)]
        check("a score above FC_MAX_LIST is excluded, as prs routes it to plink2",
              forcecall._prs_sites() == [])
    finally:
        prs._panel_scores, prs._load_variants, prs._scoring_file = (orig_panel, orig_load,
                                                                    orig_file)


def test_faidx_bad_region_does_not_truncate_batch() -> None:
    """util.faidx_bases: one unfetchable region must not take the rest of its batch with it.

    `samtools faidx` ABORTS at the first region it cannot fetch — everything requested after
    it is never emitted, while exit 1 plus partial stdout looks exactly like "those positions
    have no reference base". At 20,000 regions per batch, a handful of positions on alt
    contigs silently cost tens of thousands of real lookups, and every caller reads the gap
    through the pipeline's absent-means-hom-ref convention. Measured on the force-call site
    union: 7 alt-contig positions took 44,217 canonical ones with them (29,599 of 73,886
    sites built)."""
    print("util.faidx_bases (a bad region costs only itself):")
    import shutil
    import subprocess
    import tempfile
    from pipeline import config, util

    if shutil.which("samtools") is None:
        print("  SKIP  samtools not on PATH")
        return

    saved_ref, saved_batch = config.REF_FASTA, util._FAIDX_BATCH
    tmp = pathlib.Path(tempfile.mkdtemp())
    try:
        fa = tmp / "ref.fa"
        fa.write_text(">1\n" + "ACGT" * 50 + "\n>2\n" + "TTTT" * 50 + "\n")
        subprocess.run(["samtools", "faidx", str(fa)], check=True)
        config.REF_FASTA = fa

        # sanity: samtools really does abort at the bad region (the premise of the fix)
        raw = subprocess.run(["samtools", "faidx", str(fa), "1:5-5",
                              "NOPE:1-1", "2:9-9"], capture_output=True, text=True)
        check("premise: samtools drops what follows a bad region",
              ">2:9-9" not in raw.stdout)

        # the fix: an unknown contig is filtered up front, the good ones all resolve
        got = util.faidx_bases([("1", "5"), ("NOPE", "1"), ("2", "9")])
        check("a position after an unknown contig is still resolved",
              got.get(("2", "9")) == "T")
        check("...as is the one before it", got.get(("1", "5")) == "A")
        check("...and the unfetchable one is simply absent", ("NOPE", "1") not in got)

        # a bad region on a KNOWN contig (past its end) must also cost only itself
        got = util.faidx_bases([("1", "5"), ("1", "999999"), ("2", "9")])
        check("an out-of-range position on a known contig costs only itself",
              got.get(("1", "5")) == "A" and got.get(("2", "9")) == "T")

        # batching boundary: the retry must work inside every batch, not just the first
        util._FAIDX_BATCH = 2
        got = util.faidx_bases([("1", "1"), ("1", "5"), ("NOPE", "7"), ("2", "9"),
                                ("1", "9"), ("2", "1")])
        check("with a tiny batch size every good position still resolves",
              len([k for k in got if k[0] in ("1", "2")]) == 5)
    finally:
        config.REF_FASTA, util._FAIDX_BATCH = saved_ref, saved_batch
        shutil.rmtree(tmp, ignore_errors=True)


def test_kg_ref_pgen_names_its_outputs_correctly() -> None:
    """prs._kg_ref_pgen must look for the files plink2 actually writes.

    The build prefix carries a per-process tag as its own dot-segment
    (".build_chr1.68067"), so Path.with_suffix(".pgen") produced ".build_chr1.pgen" — a name
    plink2 never writes. Every successful build was therefore reported as
    "plink2 reference pgen build failed", every genome-wide percentile came back n/a, and
    the cleanup missed the real files too (gigabytes of orphaned builds)."""
    print("prs._kg_ref_pgen (per-process build prefix keeps its tag):")
    import shutil
    import subprocess
    import tempfile
    from pipeline import prs

    tmp = pathlib.Path(tempfile.mkdtemp())
    orig_run, orig_bin = subprocess.run, prs._plink2_bin
    try:
        # the trap itself, stated so the test explains the bug even out of context
        pref = tmp / ".build_chr1.68067"
        check("with_suffix would strip the per-process tag",
              str(pref.with_suffix(".pgen")) != str(pref) + ".pgen")

        prs._plink2_bin = lambda: "plink2-stub"
        written = {}

        def fake(cmd, **kw):
            cmd = [str(c) for c in cmd]
            if "--make-pgen" in cmd:                       # emulate plink2's real naming
                out = pathlib.Path(cmd[cmd.index("--out") + 1])
                for e in (".pgen", ".pvar", ".psam", ".log"):
                    pathlib.Path(str(out) + e).write_text("x")
                    written[e] = str(out) + e
            else:                                          # the bcftools extract
                pathlib.Path(cmd[cmd.index("-o") + 1]).write_text("x")
            return type("R", (), {"stdout": "", "returncode": 0, "stderr": ""})()
        subprocess.run = fake

        wd = tmp / "PGS_TEST_1kg"
        got = prs._kg_ref_pgen("1", [("1", "100", "A", "G", 1.0)], wd)
        check("a successful build is recognised, not reported as failed", got is not None)
        if got:
            check("...and renamed to the cached per-chromosome prefix",
                  (wd / "chr1.pgen").exists() and (wd / "chr1.pvar").exists())
        check("no tagged build files are left behind",
              not list(wd.glob(".build_chr*")))
    finally:
        subprocess.run, prs._plink2_bin = orig_run, orig_bin
        shutil.rmtree(tmp, ignore_errors=True)


def test_write_score_file_drops_duplicate_coordinates() -> None:
    """prs._write_score_file must not emit two rows at one coordinate.

    plink2 keys variants as chrom:pos (--set-all-var-ids @:#), so two scoring rows at one
    position are indistinguishable to it and it aborts the whole chromosome with
    "--score: ALT1 allele for variant '19:39893126' appears multiple times". One such
    position voids that score's genome-wide percentile: PGS003446 has exactly one in 538,001
    variants, and its percentile was n/a because of it. The reference panel already drops
    these itself via --rm-dup exclude-all, so the score file matches that."""
    print("prs._write_score_file (duplicate coordinates dropped, as the panel drops them):")
    import tempfile
    from pipeline import prs

    tmp = pathlib.Path(tempfile.mkdtemp())
    out = tmp / "score.txt"
    prs._write_score_file(out, [
        ("1", "100", "A", "G", 1.0),
        ("1", "200", "C", "T", 2.0),
        ("1", "200", "G", "T", 3.0),      # same coordinate, different allele
        ("2", "300", "T", "C", 4.0),
    ])
    lines = [ln for ln in out.read_text().splitlines()[1:] if ln.strip()]
    ids = [ln.split("\t")[0] for ln in lines]
    check("unambiguous rows are kept", set(ids) == {"1:100", "2:300"})
    check("both rows at the duplicated coordinate are dropped, not one kept arbitrarily",
          "1:200" not in ids)
    check("no ID appears twice (what plink2 aborts on)", len(ids) == len(set(ids)))
    check("the header is still written", out.read_text().startswith("ID\tA1\tW"))

    # BOTH plink2 score paths must go through this writer. _score_plink2 (the sample side)
    # used to write its file inline and keep duplicates, so plink2 aborted the whole score
    # whenever the sample carried a position the scoring file lists twice with the allele
    # plink2 knows as ALT1 — sample-dependent, so the same score failed for some samples
    # and succeeded for others.
    src = (pathlib.Path(__file__).resolve().parent.parent / "pipeline" / "prs.py").read_text()
    i = src.index("def _score_plink2")
    body = src[i:i + 1800]
    check("_score_plink2 writes its score file via _write_score_file",
          "_write_score_file(sf, variants)" in body)
    check("...and no longer writes one inline", 'fh.write("ID\\tA1\\tW' not in body)

    # the common case must be untouched
    out2 = tmp / "clean.txt"
    prs._write_score_file(out2, [("1", "1", "A", "G", 1.0), ("1", "2", "C", "T", 2.0)])
    check("a file with no duplicates keeps every row",
          len([l for l in out2.read_text().splitlines()[1:] if l.strip()]) == 2)
    import shutil
    shutil.rmtree(tmp, ignore_errors=True)


def test_reference_scores_parallel_and_targets() -> None:
    """prs._reference_scores must fan out across chromosomes and aggregate correctly, and the
    1000G extract must use -T rather than -R.

    Serially, building the 1000G panels for the 22 genome-wide panel scores was projected at
    ~38 h. Two causes: the per-chromosome work ran one at a time although each chromosome is
    fully independent, and the extract used -R, which seeks the tabix index per region. On
    chr22 (0.5 GB, 86,160 positions) -R had not finished after 13 minutes where -T took 30 s."""
    print("prs._reference_scores (parallel chromosomes) / _kg_ref_pgen (-T not -R):")
    import shutil
    import tempfile
    import threading
    import time
    from pipeline import config, prs

    src = (pathlib.Path(__file__).resolve().parent.parent / "pipeline" / "prs.py").read_text()
    i = src.index("def _kg_ref_pgen")
    body = src[i:i + 4000]
    check("the 1000G extract streams targets (-T), not index regions (-R)",
          '"-T", str(reg)' in body and '"-R", str(reg)' not in body)
    check("...and hands bcftools its thread count",
          "PRS_BCFTOOLS_THREADS" in body)

    # aggregation across concurrently-scored chromosomes
    saved = prs._ref_score_chrom
    concurrent, lock, peak = [0], threading.Lock(), [0]

    def fake(chrom, chrom_vars, workdir):
        with lock:
            concurrent[0] += 1
            peak[0] = max(peak[0], concurrent[0])
        try:
            time.sleep(0.05)     # long enough that genuine overlap is observable
            if chrom == "3":
                return chrom, None                 # one chromosome fails
            return chrom, {"HG1": float(chrom), "HG2": 10.0}
        finally:
            with lock:
                concurrent[0] -= 1

    tmp = pathlib.Path(tempfile.mkdtemp())
    saved_dir, saved_popmap = config.PGS_DIR, prs._kg_population_map
    try:
        config.PGS_DIR = tmp
        prs._ref_score_chrom = fake
        prs._kg_population_map = lambda: {"HG1": "EUR", "HG2": "EUR"}

        variants = [(c, "1", "A", "G", 1.0) for c in ("1", "2", "4", "5")]
        got = prs._reference_scores("PGS_PAR_OK", variants)
        check("scores from every chromosome are summed per individual",
              got is not None and sorted(got) == [("EUR", 12.0), ("EUR", 40.0)])
        check("chromosomes ran concurrently, not one at a time", peak[0] > 1)

        # a failed chromosome must void the whole panel, not silently shrink it
        variants = [(c, "1", "A", "G", 1.0) for c in ("1", "2", "3")]
        check("one failed chromosome voids the reference distribution",
              prs._reference_scores("PGS_PAR_BAD", variants) is None)
    finally:
        prs._ref_score_chrom, prs._kg_population_map = saved, saved_popmap
        config.PGS_DIR = saved_dir
        shutil.rmtree(tmp, ignore_errors=True)


def test_consent_revocation_and_per_person_gate() -> None:
    """Revoking consent must reach the report that disclosed it, and the PRS gate must be
    per person rather than one global switch for the whole family.

    set_consent updated the ledger and stopped there, so the previously generated report sat
    on disk with the findings still in it — readable by the dashboard, by render and by anything reading the reports tree.
    And repeats reused that report on existence alone, so a consent change in either
    direction left the old disclosure in place until someone remembered --force. Separately,
    prs._panel_filters consulted only config.PRS_S1_GATE: one flag decided disclosure for
    everyone, although consent.tsv records it per (sample, category)."""
    print("consent.set_consent (revocation reaches reports) / prs._panel_filters (per person):")
    import shutil
    import tempfile
    from pipeline import config, consent, prs

    saved = {n: getattr(config, n) for n in ("CONSENT_FILE", "REPEATS_DIR", "PRS_S1_GATE",
                                             "REPORTS_DIR", "FAMILY_FILE", "PRS_DIR")}
    tmp = pathlib.Path(tempfile.mkdtemp())
    try:
        config.CONSENT_FILE = tmp / "consent.tsv"
        config.REPEATS_DIR = tmp / "repeats"
        config.PRS_DIR = tmp / "prs"
        config.PRS_S1_GATE = True
        config.REPORTS_DIR = tmp / "reports"
        config.FAMILY_FILE = tmp / "family.tsv"
        cat = consent.CATEGORY_ADULT_ONSET_UNTREATABLE

        # a granted consent, and a report generated under it
        consent.set_consent("SAMPLE001", cat, True)
        check("consent is recorded", consent.has_consented("SAMPLE001", cat) is True)
        gated = consent.gated_outputs("SAMPLE001", cat)
        md, html = gated[0], gated[1]
        for f in gated:
            f.parent.mkdir(parents=True, exist_ok=True)
            f.write_text("FMR1 premutation, 92 repeats\n")
        check("gated_outputs names the files the decision governs: the repeats and PRS "
              "reports in each format tree, and their stamps",
              gated[0] == md and len(gated) == 2 * (2 + len(render_formats_count())) + 1
              and "md" in md.parts and "html" in html.parts)

        consent.set_consent("SAMPLE001", cat, False)
        check("revocation is recorded", consent.has_consented("SAMPLE001", cat) is False)
        check("...and the report that disclosed it is gone",
              not md.exists())
        check("...including its rendered HTML", not html.exists())

        # an unrelated category must not touch it
        md.write_text("x")
        consent.set_consent("SAMPLE001", "some_other_category", False)
        check("revoking a different category leaves it alone", md.exists())
        # nor another person's report
        consent.set_consent("SAMPLE002", cat, False)
        check("revoking for another person leaves it alone", md.exists())

        # --- the PRS gate, per person ---
        rows = [{"pgs_id": "PGS1", "s1_gate": "yes", "n_variants": "10"},
                {"pgs_id": "PGS2", "s1_gate": "no", "n_variants": "10"}]
        config.PRS_S1_GATE = True
        consent.set_consent("SAMPLE001", cat, True)
        consent.set_consent("SAMPLE002", cat, False)
        kept_m, _ = prs._panel_filters(list(rows), "SAMPLE001")
        kept_a, notes_a = prs._panel_filters(list(rows), "SAMPLE002")
        check("a consenting person sees the gated trait", len(kept_m) == 2)
        check("a non-consenting person does not", len(kept_a) == 1)
        check("...and the note names who, not just the flag",
              any("SAMPLE002" in n for n in notes_a))

        config.PRS_S1_GATE = False
        kept_a2, _ = prs._panel_filters(list(rows), "SAMPLE002")
        check("the global switch off ungates everyone, as documented", len(kept_a2) == 2)
        check("...and with it off, revocation has no PRS output to remove",
              not any("prs_" in f.name for f in consent.gated_outputs("SAMPLE001", cat)))

        # --- revocation reaches the PRS report ---
        # Scored while consent stood, so it lists the gated trait. It used to survive the
        # revocation in every format, the dashboard kept quoting its percentile, and a
        # plain `prs` re-run returned it before the consent filter was ever consulted.
        from pipeline import dashboard, reportpaths, util
        config.PRS_S1_GATE = True
        consent.set_consent("SAMPLE001", cat, True)
        table = ["# Polygenic risk scores — SAMPLE001", "", "## Not actionable", "",
                 "| Trait | PGS | Score | EUR percentile |", "|---|---|---|---|",
                 "| Alzheimer disease | PGS1 | 1.2 | 97% |"]

        def scored():
            report = util.write_report(reportpaths.genome_report("prs", "SAMPLE001"), table)
            prs.consent_stamp("SAMPLE001").parent.mkdir(parents=True, exist_ok=True)
            prs.consent_stamp("SAMPLE001").write_text(prs.s1_token("SAMPLE001") + "\n")
            prs.panel_stamp("SAMPLE001").write_text(prs.panel_token() + "\n")
            return report
        pmd = scored()
        twins = [reportpaths.derived(pmd, name, ext)
                 for name, (ext, _r) in render_formats_count().items()]
        check("a consented PRS report is on disk in every format, and reusable",
              pmd.exists() and all(t.exists() for t in twins)
              and prs.s1_token("SAMPLE001") == "s1=shown"
              and prs._report_is_current("SAMPLE001"))
        check("...and the dashboard reads its percentile",
              dashboard.prs_percentiles("SAMPLE001") == [("Alzheimer disease", 97.0)])
        consent.set_consent("SAMPLE001", cat, False)
        check("revoking removes the PRS report in every format",
              not pmd.exists() and not any(t.exists() for t in twins))
        check("...and its stamps", not prs.consent_stamp("SAMPLE001").exists()
              and not prs.panel_stamp("SAMPLE001").exists())
        check("...so the dashboard has no percentile left to show",
              dashboard.prs_percentiles("SAMPLE001") == [])

        # A report that outlives the change some other way (restored, synced from another
        # machine) must still not be reused under the new state.
        consent.set_consent("SAMPLE001", cat, True)
        scored()
        (config.CONSENT_FILE).write_text(
            config.CONSENT_FILE.read_text().replace("SAMPLE001\t" + cat + "\tyes",
                                                    "SAMPLE001\t" + cat + "\tno"))
        check("a report stamped 'shown' is not reused once consent reads 'no'",
              consent.has_consented("SAMPLE001", cat) is False
              and not prs._report_is_current("SAMPLE001"))
        prs.consent_stamp("SAMPLE001").unlink()
        check("an unstamped report is regenerated while traits must be withheld",
              not prs._report_is_current("SAMPLE001"))
        consent.set_consent("SAMPLE001", cat, True)
        scored()
        prs.consent_stamp("SAMPLE001").unlink()
        check("...and reused while nothing is withheld (it cannot over-disclose)",
              prs._report_is_current("SAMPLE001"))
    finally:
        for n, v in saved.items():
            setattr(config, n, v)
        shutil.rmtree(tmp, ignore_errors=True)


def test_prs_panel_stamp_and_rescoring() -> None:
    """A PRS report is only as current as the panel it was scored against.

    Registering a score is an edit to panels/complex_traits.tsv. Nothing noticed: `prs`
    returned the cached report on existence alone, and the scheduled job never looked at
    polygenic scores at all, so a new score reached nobody until someone remembered
    --force for every callset. The report now carries a stamp of the panel, and
    `check-releases` re-scores the callsets whose stamp no longer matches."""
    print("prs panel stamp / check-releases re-scores stale reports:")
    import argparse
    import shutil
    import tempfile
    import run as RUN
    from pipeline import prs, releases, reportpaths, util

    names = ("REPORTS_DIR", "PRS_DIR", "COMPLEX_PANEL", "PRS_MAX_VARIANTS", "PRS_S1_GATE",
             "SAMPLES_TSV", "CONSENT_FILE", "FAMILY_FILE", "PRS_AUTO_RESCORE", "NOTIFY_CMD")
    saved = {n: getattr(config, n) for n in names}
    orig = (releases.has_new_release, prs.prs_report, RUN._build_dashboards)
    tmp = pathlib.Path(tempfile.mkdtemp())
    try:
        config.REPORTS_DIR = tmp / "reports"
        config.PRS_DIR = tmp / "prs"
        config.COMPLEX_PANEL = tmp / "panel.tsv"
        config.SAMPLES_TSV = tmp / "samples.tsv"
        config.CONSENT_FILE = tmp / "consent.tsv"
        config.FAMILY_FILE = tmp / "family.tsv"
        config.PRS_MAX_VARIANTS = None
        config.PRS_S1_GATE = False
        config.PRS_AUTO_RESCORE = True
        config.NOTIFY_CMD = ""
        config.SAMPLES_TSV.write_text("sample_id\nJan\nMaria\nLena\n")
        head = "pgs_id\ttrait\ts1_gate\tn_variants\n"
        config.COMPLEX_PANEL.write_text(head + "PGS1\tType 2 diabetes\tno\t10\n")

        def scored(sample):
            util.write_report(reportpaths.genome_report("prs", sample),
                              [f"# Polygenic risk scores — {sample}"])
            prs.consent_stamp(sample).parent.mkdir(parents=True, exist_ok=True)
            prs.consent_stamp(sample).write_text(prs.s1_token(sample) + "\n")
            prs.panel_stamp(sample).write_text(prs.panel_token() + "\n")
        scored("Jan")
        scored("Maria")      # Lena has never been scored
        t0 = prs.panel_token()
        check("a report scored against the current panel is current",
              prs._report_is_current("Jan") and prs.stale_reports(["Jan", "Maria", "Lena"]) == [])

        config.COMPLEX_PANEL.write_text(head + "PGS1\tType 2 diabetes\tno\t10\n"
                                        "PGS2\tCoronary artery disease\tno\t10\n")
        check("registering a score changes the panel token", prs.panel_token() != t0)
        check("...which makes every scored callset stale, and only those",
              prs.stale_reports(["Jan", "Maria", "Lena"]) == ["Jan", "Maria"])
        check("a callset that was never scored is not picked up",
              "Lena" not in prs.stale_reports(["Lena"]))

        scored("Jan")
        prs.panel_stamp("Jan").unlink()
        check("a report with no panel stamp is stale (scored against an unknown panel)",
              prs.stale_reports(["Jan"]) == ["Jan"])
        config.PRS_MAX_VARIANTS = 50000
        scored("Jan")
        config.PRS_MAX_VARIANTS = None
        check("a capped (partial) report is stale for an uncapped run",
              prs.stale_reports(["Jan"]) == ["Jan"])

        # check-releases: nothing new upstream, but two stale PRS reports
        calls = []

        def fake_prs(sample, force=False):
            calls.append(sample)
            if sample == "Maria":
                raise SystemExit("scoring file download failed")
            scored(sample)
        releases.has_new_release = lambda: False
        prs.prs_report = fake_prs
        RUN._build_dashboards = lambda snap=None: None
        rc = RUN.cmd_check_releases(argparse.Namespace(full=False))
        check("check-releases re-scores the stale callsets without a ClinVar release",
              calls == ["Jan", "Maria"])
        check("one failure does not stop the others, and is reported in the exit code",
              rc == 1 and prs.stale_reports(["Jan", "Maria"]) == ["Maria"])
        calls.clear()
        config.PRS_AUTO_RESCORE = False
        check("PRS_AUTO_RESCORE=0 leaves them alone",
              RUN.cmd_check_releases(argparse.Namespace(full=False)) == 0 and calls == [])
        config.PRS_AUTO_RESCORE = True
        prs.prs_report = lambda sample, force=False: (calls.append(sample), scored(sample))
        check("...and the retry clears the backlog",
              RUN.cmd_check_releases(argparse.Namespace(full=False)) == 0
              and calls == ["Maria"]
              and RUN.cmd_check_releases(argparse.Namespace(full=False)) == 0
              and calls == ["Maria"])
    finally:
        releases.has_new_release, prs.prs_report, RUN._build_dashboards = orig
        for n, v in saved.items():
            setattr(config, n, v)
        shutil.rmtree(tmp, ignore_errors=True)


def test_cache_keys_are_content_based() -> None:
    """Identity for cache decisions must be content, not (size, mtime), and --force must
    actually discard resumable shards.

    The old proxy argued a content change at identical size and whole-second mtime was
    unrealizable. True of an edit in place; false of `rsync -a` and `cp -p`, which preserve
    mtime by design — a file restored from backup or copied between hosts presents the
    previous token with different bytes. Separately, force_call keyed only on the site set,
    so a re-aligned sample kept the genotypes called off its previous alignment, and
    call_at_sites returned existing shards even under --force, so the CRAM was never re-read."""
    print("util.file_token / incremental.source_digests (content) / --force discards shards:")
    import shutil
    import tempfile
    from pipeline import forcecall, util

    tmp = pathlib.Path(tempfile.mkdtemp())
    try:
        # same size, same mtime, different bytes — the case the proxy could not see
        a, b = tmp / "a.bin", tmp / "b.bin"
        a.write_bytes(b"A" * 4096)
        b.write_bytes(b"B" * 4096)
        import os
        os.utime(b, (1_700_000_000, 1_700_000_000))
        os.utime(a, (1_700_000_000, 1_700_000_000))
        check("same size + same mtime, different content -> different tokens",
              util.file_token(a) != util.file_token(b))
        check("identical content -> identical tokens",
              util.file_token(a) == util.file_token(a))

        # source_digests must be content-based too
        src = (pathlib.Path(__file__).resolve().parent.parent
               / "pipeline" / "incremental.py").read_text()
        i = src.index("def source_digests")
        body = src[i:i + 1600]
        check("source_digests hashes the normalized VCF",
              "sha256(normalized_path(sample))" in body)
        check("...and no longer fingerprints it by size:mtime",
              "st.st_size" not in body)

        # --force must clear the shard directory
        fsrc = (pathlib.Path(__file__).resolve().parent.parent
                / "pipeline" / "forcecall.py").read_text()
        j = fsrc.index("def call_at_sites")
        cbody = fsrc[j:j + 1400]
        check("call_at_sites takes a force flag", "force: bool = False" in cbody)
        check("...and clears resumable shards when it is set",
              "if force or not parts_key.exists()" in cbody)
        check("force_call keys on the CRAM as well as the site set",
              '"|cram=" + file_token(cram)' in fsrc)

        ccall = (pathlib.Path(__file__).resolve().parent.parent
                 / "pipeline" / "call.py").read_text()
        check("call() also discards its shards under --force",
              "if force:" in ccall and "for stale in parts_dir.glob" in ccall)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_pinned_vep_version_and_snapshot_lock() -> None:
    """annotate must use the snapshot's pinned VEP cache version, and snapshot creation must
    serialise concurrent creators of the same id.

    _run_vep always passed config.VEP_CACHE_VERSION, so re-annotating an older snapshot after
    a cache bump produced rows disagreeing with the version recorded in that snapshot's own
    manifest and in runs.tsv. And create_snapshot's per-process .part stopped two runs
    interleaving bytes but not racing on everything after: both renamed onto the final
    clinvar.vcf.gz, both ran tabix (whose .tbi write is not atomic), both wrote the manifest."""
    print("annotate._pinned_cache_version / snapshot creation lock:")
    import json
    import shutil
    import tempfile
    from pipeline import annotate as ann
    from pipeline import config, snapshot

    saved = {n: getattr(config, n) for n in
             ("SNAPSHOTS_DIR", "VEP_CACHE_VERSION", "VEP_CACHE_DIR", "BUILD")}
    tmp = pathlib.Path(tempfile.mkdtemp())
    try:
        config.SNAPSHOTS_DIR = tmp / "snap"
        config.VEP_CACHE_DIR = tmp / "vep_cache"
        config.VEP_CACHE_VERSION = "112"
        config.BUILD = "GRCh38"

        def mk(sid, ver):
            d = config.SNAPSHOTS_DIR / sid
            d.mkdir(parents=True, exist_ok=True)
            payload = {"snapshot_id": sid}
            if ver is not None:
                payload["vep_cache_version"] = ver
            (d / "manifest.json").write_text(json.dumps(payload))

        mk("2026-09-01", "112")
        check("a snapshot pinned to the configured version uses it",
              ann._pinned_cache_version("2026-09-01") == "112")

        # pinned to an older cache that IS installed -> use the pinned one
        mk("2026-06-01", "110")
        (config.VEP_CACHE_DIR / "homo_sapiens" / "110_GRCh38").mkdir(parents=True)
        check("an older pinned version is used when its cache is installed",
              ann._pinned_cache_version("2026-06-01") == "110")

        # pinned to a cache that is NOT installed -> refuse rather than substitute
        mk("2026-05-01", "108")
        try:
            ann._pinned_cache_version("2026-05-01")
            check("a missing pinned cache refuses rather than substituting", False)
        except SystemExit as e:
            check("a missing pinned cache refuses rather than substituting", True)
            check("...naming both versions", "108" in str(e) and "112" in str(e))

        mk("2026-04-01", None)
        check("a manifest with no recorded version falls back to config",
              ann._pinned_cache_version("2026-04-01") == "112")

        # the creation lock exists and is exclusive per snapshot id
        src = (pathlib.Path(__file__).resolve().parent.parent
               / "pipeline" / "snapshot.py").read_text()
        check("create_snapshot takes an exclusive lock", "fcntl.flock(lock, fcntl.LOCK_EX)" in src)
        check("...and re-checks the manifest under it, so the loser reuses",
              "was created by a concurrent run" in src)
    finally:
        for n, v in saved.items():
            setattr(config, n, v)
        shutil.rmtree(tmp, ignore_errors=True)


def test_demo_trio_fixtures() -> None:
    """Validate that the synthetic 3-generation family cohort fixtures are well-formed VCFs."""
    print("demo_trio_fixtures (well-formed VCFs, headers, traits, pedigree):")
    import pathlib
    from pipeline import ingest, pedigree, traits

    fixtures = pathlib.Path(__file__).resolve().parent / "fixtures" / "trio"
    check("fixture directory exists", fixtures.exists())
    samples = ("father", "mother", "daughter", "partner", "grandson", "granddaughter", "child")
    for s in samples:
        vcf_gz = fixtures / f"{s}.vcf.gz"
        tbi = fixtures / f"{s}.vcf.gz.tbi"
        check(f"{s}.vcf.gz exists", vcf_gz.exists())
        check(f"{s}.vcf.gz.tbi exists", tbi.exists())
        hdr = ingest._header_lines(vcf_gz)
        ok, _ = ingest._looks_like_build(hdr)
        check(f"{s} passes GRCh38 build check", ok)
        check(f"{s} has exactly 1 sample column", ingest._sample_names(hdr) == [s])

    # Single-SNP trait queries
    check("father has lactase persistence allele",
          traits._alleles_at(fixtures / "father.vcf.gz", "chr2", "135851076", "G") == ("G", "A"))
    check("mother is homozygous reference (lactose intolerant)",
          traits._alleles_at(fixtures / "mother.vcf.gz", "chr2", "135851076", "G") == ("G", "G"))
    check("daughter inherited lactase persistence allele",
          traits._alleles_at(fixtures / "daughter.vcf.gz", "chr2", "135851076", "G") == ("G", "A"))
    check("child alias inherited lactase persistence allele",
          traits._alleles_at(fixtures / "child.vcf.gz", "chr2", "135851076", "G") == ("G", "A"))

    check("mother has blue eye alleles",
          traits._alleles_at(fixtures / "mother.vcf.gz", "chr15", "28120472", "A") == ("G", "G"))
    check("daughter has blue eye alleles",
          traits._alleles_at(fixtures / "daughter.vcf.gz", "chr15", "28120472", "A") == ("G", "G"))

    # APOB (Familial Hypercholesterolemia, chr2:21002409 G>A) AD transmission
    check("father carries APOB mutation (het)",
          traits._alleles_at(fixtures / "father.vcf.gz", "chr2", "21002409", "G") == ("G", "A"))
    check("mother is homozygous reference for APOB",
          traits._alleles_at(fixtures / "mother.vcf.gz", "chr2", "21002409", "G") == ("G", "G"))
    check("daughter inherited APOB mutation (het)",
          traits._alleles_at(fixtures / "daughter.vcf.gz", "chr2", "21002409", "G") == ("G", "A"))
    check("partner founder is homozygous reference for APOB",
          traits._alleles_at(fixtures / "partner.vcf.gz", "chr2", "21002409", "G") == ("G", "G"))
    check("grandson inherited APOB mutation from daughter (het)",
          traits._alleles_at(fixtures / "grandson.vcf.gz", "chr2", "21002409", "G") == ("G", "A"))
    check("granddaughter inherited reference allele for APOB (spared)",
          traits._alleles_at(fixtures / "granddaughter.vcf.gz", "chr2", "21002409", "G") == ("G", "G"))

    fam = pedigree.load(fixtures / "family.tsv")
    check("pedigree resolves trio for daughter",
          fam.trio("daughter") == ("daughter", "father", "mother"))
    check("pedigree resolves trio for grandson",
          fam.trio("grandson") == ("grandson", "partner", "daughter"))
    check("pedigree resolves trio for granddaughter",
          fam.trio("granddaughter") == ("granddaughter", "partner", "daughter"))
    check("pedigree relationship daughter -> father",
          fam.relationship_to("daughter", "father") == "father")
    check("pedigree relationship father -> daughter",
          fam.relationship_to("father", "daughter") == "daughter")
    check("pedigree relationship father -> grandson",
          fam.relationship_to("father", "grandson") == "grandson")
    check("pedigree relationship father -> granddaughter",
          fam.relationship_to("father", "granddaughter") == "granddaughter")
    check("pedigree relationship granddaughter -> father",
          fam.relationship_to("granddaughter", "father") == "grandfather")
    check("pedigree relationship partner -> father",
          fam.relationship_to("partner", "father") == "father-in-law")
    check("pedigree relationship father -> partner",
          fam.relationship_to("father", "partner") == "son-in-law")
    check("pedigree relationship grandson -> granddaughter",
          fam.relationship_to("grandson", "granddaughter") == "sister")


def test_cmd_demo() -> None:
    """Validate that `run.py demo` runs cleanly end-to-end with the 3-generation fixtures."""
    print("cmd_demo (zero-data end-to-end pipeline run):")
    import argparse
    import pathlib
    import tempfile
    import run

    with tempfile.TemporaryDirectory(prefix="genome_demo_test_") as tmp:
        args = argparse.Namespace(out_dir=tmp, clean=True)
        ret = run.cmd_demo(args)
        check("cmd_demo exits with 0", ret == 0)

        tmp_p = pathlib.Path(tmp)
        # Check that reports were created
        check("demo family pedigree.ped exists", (tmp_p / "pedigree.ped").exists())
        check("demo samples.tsv exists", (tmp_p / "manifest" / "samples.tsv").exists())
        check("demo DuckDB database exists", (tmp_p / "duckdb" / "genomes.duckdb").exists())

        # Check trait reports and HTML dashboards for all 6 members
        cohort = ("father", "mother", "daughter", "partner", "grandson", "granddaughter")
        for s in cohort:
            check(f"{s} traits report exists",
                  (tmp_p / "reports" / "md" / "genome" / s / f"traits_{s}.md").exists())

        md_root, html_root = tmp_p / "reports" / "md", tmp_p / "reports" / "html"
        genome_md = md_root / "genome"
        check("md/ holds the latest snapshot only; the earlier one is in md-history/",
              [d.name for d in md_root.iterdir() if d.name != "genome"]
              == [max(d.name for d in html_root.glob("*-demo"))]
              and len(list((tmp_p / "reports" / "md-history").glob("*-demo"))) == 1)
        check("stage directories hold no reports",
              not list(tmp_p.glob("*/*/*.md")) or all(
                  "reports" in p.parts for p in tmp_p.glob("*/*/*.md")))
        check("every report opens with front matter, then its title",
              all(p.read_text().startswith("---\ntitle: ") and "\n---\n\n# " in p.read_text()
                  for p in (tmp_p / "reports").rglob("*.md")))

        def _last_h2(p):
            return [ln for ln in p.read_text().splitlines() if ln.startswith("## ")][-1]
        from pipeline.util import CLOSING_HEADINGS
        bad = [p.name for p in (tmp_p / "reports").rglob("*.md")
               if _last_h2(p) not in CLOSING_HEADINGS]
        check(f"every report ends with a closing caveats section (not: {bad})", not bad)

        # Check panel reports
        report_dirs = list(html_root.glob("*-demo"))
        check("demo report dirs exist: the snapshot and its earlier one", len(report_dirs) == 2)
        if report_dirs:
            snap_rep, earlier = max(report_dirs), min(report_dirs)
            now_page = (snap_rep / "index_father.html").read_text()
            then_page = (earlier / "index_father.html").read_text()
            check("the tree offers the other snapshot, each way",
                  f'href="../{earlier.name}/index_father.html"' in now_page
                  and f'href="../{snap_rep.name}/index_father.html"' in then_page)
            check("...and so do the snapshot index and the inheritance page",
                  f"../{earlier.name}/index.html" in (snap_rep / "index.html").read_text()
                  and f"../{snap_rep.name}/inheritance.html"
                  in (earlier / "inheritance.html").read_text())
            check("the earlier snapshot shows the earlier classification",
                  "Pathogenic" not in "".join(
                      ln for ln in (tmp_p / "reports" / "md-history" / earlier.name / "father"
                                    / "full_father.md").read_text().splitlines()
                      if "APOB" in ln))
            for s in cohort:
                check(f"index_{s}.html exists", (snap_rep / f"index_{s}.html").exists())
            check("family inheritance page exists", (snap_rep / "inheritance.html").exists())
            check("reports landing page exists", (html_root / "index.html").exists())
            snap_md = md_root / snap_rep.name
            for who in ("daughter", "grandson", "granddaughter"):
                check(f"panel_{who}.md exists", (snap_md / who / f"panel_{who}.md").exists())
            check("panel_family.md exists", (snap_md / "_family" / "panel_family.md").exists())
            check("the inheritance findings exist as markdown too",
                  "| APOB |" in (snap_md / "_family" / "inheritance.md").read_text())

            daughter_content = (snap_md / "daughter" / "panel_daughter.md").read_text()
            check("daughter panel reports Hemochromatosis", "Hemochromatosis" in daughter_content)
            check("daughter panel reports Familial hypercholesterolemia (APOB)",
                  "Familial hypercholesterolemia" in daughter_content or "APOB" in daughter_content)

            grandson_content = (snap_md / "grandson" / "panel_grandson.md").read_text()
            check("grandson panel reports APOB", "APOB" in grandson_content)

            check("every demo report says it is demo output",
                  "Demo only" in daughter_content.split("\n---\n", 1)[1].splitlines()[3])
            check("so do the dashboards",
                  "Demo only" in (snap_rep / "index_daughter.html").read_text())

        # HLA is read from the demo genomes for real; the rest are bundled sample results
        # run through the real report code and labelled as illustrative.
        hla = (genome_md / "grandson" / "hla_grandson.md").read_text()
        check("HLA report is computed from the fixture genotype", "DQ2.5 present (1 copy)" in hla)
        check("...and is not labelled illustrative", "illustrative values" not in hla)
        for stage in ("repeats", "pgx", "ancestry", "prs"):
            for s in cohort:
                md = genome_md / s / f"{stage}_{s}.md"
                check(f"{stage} report for {s} exists and is labelled illustrative",
                      md.exists() and "Demo only — illustrative values" in md.read_text())
        rep_m = (genome_md / "mother" / "repeats_mother.md").read_text()
        check("toy repeat sizes go through the real thresholds",
              "| FMR1 | CGG | 29/58 | premutation" in rep_m)
        prs_f = (genome_md / "father" / "prs_father.md").read_text()
        check("toy scores go through the real banding", "91.0% | high (top 10%)" in prs_f)
        check("the demo family sees the consent-gated trait", "multiple sclerosis" in prs_f)
        anc = (genome_md / "grandson" / "ancestry_grandson.md").read_text()
        check("haplogroups follow the lineages", "H1a1" in anc and "I-M253" in anc)
        check("no Y haplogroup for a woman", "not applicable to females" in
              (genome_md / "mother" / "ancestry_mother.md").read_text())
        check("the dashboard links the new reports",
              "Polygenic risk scores" in (snap_rep / "index_father.html").read_text())
        deltas = list(md_root.glob("*-demo/_family/delta_*.md"))
        check("a DELTA report is produced against an earlier demo snapshot", len(deltas) == 1)
        if deltas:
            delta = deltas[0].read_text()
            check("the real diff finds the reclassification",
                  "'Uncertain_significance' → 'Pathogenic'" in delta and "APOB" in delta)
            check("...for each carrier and nobody else",
                  delta.count("| APOB |") == 3 and "| partner |" not in delta)
            check("the DELTA says both snapshots are invented", "synthetic history" in delta)
            check("the dashboard links it",
                  "Latest changes (DELTA)" in (snap_rep / "index_father.html").read_text())
        check("the demo leaves the notice off for the rest of the process",
              config.DEMO_NOTICE == "")

        if report_dirs:
            family_content = (md_root / snap_rep.name / "_family" / "panel_family.md").read_text()
            check("family panel reports reproductive risk",
                  "Hemochromatosis" in family_content or "HFE" in family_content)


def test_demo_notice_and_callouts() -> None:
    """A `> ` line renders as a callout, and DEMO_NOTICE lands under a report's title."""
    print("render callouts + demo notice:")
    import tempfile
    from pipeline import render, util

    html = render._convert_body("# T\n\n> **Note** one\n> two\n\nbody")
    check("quoted lines become one callout",
          '<div class="disclaimer"><strong>Note</strong> one two</div>' in html)
    check("...and no literal '>' leaks into the page", "&gt;" not in html)
    h = render._convert_body("## Two\n\n### Three\n\ntext")
    check("a third-level heading renders as one, not as literal hashes",
          "<h3>Three</h3>" in h and "<h2>Two</h2>" in h and "###" not in h)
    saved = config.DEMO_NOTICE
    try:
        with tempfile.TemporaryDirectory() as d:
            p = pathlib.Path(d) / "r.md"
            config.DEMO_NOTICE = ""
            util.write_report(p, ["# Title", "", "body"])
            check("no notice unless the demo sets one", p.read_text() == "# Title\n\nbody")
            config.DEMO_NOTICE = "**Demo only.**"
            util.write_report(p, ["# Title", "", "body"])
            check("notice sits directly under the title",
                  p.read_text().splitlines()[:3] == ["# Title", "", "> **Demo only.**"])
            check("...and reaches the HTML as a callout",
                  'class="disclaimer"><strong>Demo only.' in p.with_suffix(".html").read_text())
    finally:
        config.DEMO_NOTICE = saved


def test_prs_guide_lists_the_panel() -> None:
    """docs/prs-guide.md carries a table of the panel; it must not drift from the TSV."""
    print("docs/prs-guide.md mirrors panels/complex_traits.tsv:")
    import csv
    root = pathlib.Path(__file__).resolve().parent.parent
    guide = (root / "docs" / "prs-guide.md").read_text()
    with open(root / "panels" / "complex_traits.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    listed = [ln for ln in guide.splitlines() if ln.startswith("| ") and "pgscatalog.org/score/" in ln]
    check("one table row per panel row", len(listed) == len(rows))
    missing = [r["pgs_id"] for r in rows
               if not any(f"[{r['pgs_id']}]" in ln and f"| {r['trait']} |" in ln
                          and f"{int(r['n_variants']):,}" in ln for ln in listed)]
    check(f"every score is listed with its trait and size (missing: {missing})", not missing)
    gated = {r["pgs_id"] for r in rows if r["s1_gate"] == "yes"}
    check("consent-gated flags match",
          {ln.split("[")[1].split("]")[0] for ln in listed if ln.rstrip().endswith("| yes |")} == gated)


def test_perspective_switch() -> None:
    """The "Viewing as" picker: relationships for every pair, embedded for a static page."""
    print("render perspective picker (static, browser-side ego):")
    import json
    from pipeline import render
    from pipeline.pedigree import Family, Person

    fam = Family({
        "gp": Person("gp", "Gran </script>", "female", None, None, None),
        "mum": Person("mum", "Mum", "female", None, "gp", None),
        "kid": Person("kid", "Kid", "male", None, "mum", None),
    }, sequenced={"gp", "mum", "kid"})
    script = render.perspective_script(fam)
    rel = json.loads(script.split("const REL = ", 1)[1].split(", NAMES = ", 1)[0])
    check("a relationship for every ordered pair",
          all(set(rel[e]) == set(fam.people) - {e} for e in fam.people))
    check("...computed from each person's own point of view",
          rel["kid"]["gp"] == "grandmother" and rel["gp"]["kid"] == "grandson"
          and rel["mum"]["kid"] == "son")
    check("a name cannot close the script element", script.count("</script>") == 1)
    control = render.perspective_control(fam)
    check("the picker offers everyone, and no one",
          control.count("data-ego-pick=") == 4 and 'data-ego-pick=""' in control)
    page = render._landing_page(fam, "2026-01-02", "", "Genome reports", "sub")
    check("the landing roster has a relationship slot per person",
          all(f'data-rel-of="{p}"' in page for p in fam.people))
    check("...and carries the picker and its script",
          "Viewing as:" in page and "localStorage" in page)
    from pipeline import serve
    check("the server no longer builds its own landing page", not hasattr(serve, "_landing"))


def test_report_formats_and_migration() -> None:
    """A registered format lands in its own tree; the old layout migrates by moving files."""
    print("render.register_format / migrate_reports (old layout -> reports/md):")
    import json
    import tempfile
    from pipeline import migrate_reports, render, reportpaths, util

    root = pathlib.Path(tempfile.mkdtemp())
    names = ["GENOMES_ROOT", "REPORTS_DIR", "SNAPSHOTS_DIR", "FAMILY_FILE", "PGX_DIR",
             "PRS_DIR", "INTERPRET_DIR"]
    saved = {n: getattr(config, n) for n in names}
    saved_formats = dict(render.FORMATS)
    try:
        config.GENOMES_ROOT = root
        for n in names[1:]:
            setattr(config, n, root / n.lower().replace("_dir", ""))
        config.REPORTS_DIR = root / "reports"
        for snap in ("2026-01-01", "2026-02-01"):
            d = config.SNAPSHOTS_DIR / snap
            d.mkdir(parents=True)
            (d / "manifest.json").write_text(json.dumps({"clinvar_date": snap}))

        # --- a second format, registered the way docs/report-formats.md describes
        render.register_format("outline", "json", lambda body, fields, title: json.dumps(
            {**fields, "title": title, "sections": [ln[3:] for ln in body.splitlines()
                                                    if ln.startswith("## ")]}))
        md = util.write_report(reportpaths.genome_report("prs", "S1"),
                               ["# PRS — S1", "", "- a result", "", "## Notes", "- n"])
        out = config.REPORTS_DIR / "outline" / "genome" / "S1" / "prs_S1.json"
        check("a registered format is written to reports/<format>/<same path>", out.exists())
        doc = json.loads(out.read_text())
        check("...and receives the front matter and the body",
              doc["report"] == "prs" and doc["person"] == "S1" and doc["sections"] == ["Caveats"])
        check("'Notes' is filed as the shared closing section", "## Caveats" in md.read_text())
        check("html is still produced beside it",
              (config.REPORTS_DIR / "html" / "genome" / "S1" / "prs_S1.html").exists())
        html = (config.REPORTS_DIR / "html" / "genome" / "S1" / "prs_S1.html").read_text()
        check("front matter is shown as a line, never as raw YAML",
              "report: prs" not in html and "generated " in html)
        render.FORMATS.pop("outline")

        # --- the old layout: dated dir with md+html side by side, reports in stage dirs
        old = config.REPORTS_DIR / "2026-01-01"
        old.mkdir(parents=True)
        for f in ("full_S1.md", "panel_S1.md", "panel_family.md",
                  "delta_2025-12-01_to_2026-01-01.md"):
            (old / f).write_text("# Old report\n\nbody\n")
            (old / f).with_suffix(".html").write_text("<p>old</p>")
        (old / "index_S1.html").write_text("x")
        new = config.REPORTS_DIR / "2026-02-01"
        new.mkdir()
        (new / "full_S1.md").write_text("# New report\n")
        (config.PGX_DIR / "S1").mkdir(parents=True)
        (config.PGX_DIR / "S1" / "pgx_S1.md").write_text("# PGx\n")
        (config.PGX_DIR / "S1" / "pgx_S1.html").write_text("<p>old</p>")
        (config.PGX_DIR / "S1" / "S1.chr.report.json").write_text("{}")
        config.INTERPRET_DIR.mkdir(parents=True)
        (config.INTERPRET_DIR / "interpret_enformer_2026-02-01.md").write_text("# I\n")

        moves, deletions = migrate_reports.plan()
        check("the plan finds every old report", len(moves) == 7)
        migrate_reports.migrate(apply=False)
        check("a dry run moves nothing", (old / "full_S1.md").exists())
        migrate_reports.migrate(apply=True)
        check("the latest snapshot lands in md/, the older one in md-history/",
              (config.REPORTS_DIR / "md" / "2026-02-01" / "S1" / "full_S1.md").exists()
              and (config.REPORTS_DIR / "md-history" / "2026-01-01" / "S1" / "full_S1.md").exists())
        check("family-level reports go under _family/",
              (config.REPORTS_DIR / "md-history" / "2026-01-01" / "_family"
               / "panel_family.md").exists()
              and (config.REPORTS_DIR / "md" / "2026-02-01" / "_family"
                   / "interpret_enformer.md").exists())
        moved = config.REPORTS_DIR / "md" / "genome" / "S1" / "pgx_S1.md"
        check("a stage report moves out of its working directory, with front matter added",
              moved.exists() and moved.read_text().startswith("---\ntitle: \"PGx\"")
              and moved.read_text().endswith("# PGx\n"))
        check("the stage's working files stay where they were",
              (config.PGX_DIR / "S1" / "S1.chr.report.json").exists())
        check("old derived HTML and the emptied dated directories are gone",
              not old.exists() and not (config.PGX_DIR / "S1" / "pgx_S1.html").exists())
        check("a second run finds nothing to do", migrate_reports.plan() == ([], []))
    finally:
        render.FORMATS.clear()
        render.FORMATS.update(saved_formats)
        for n, v in saved.items():
            setattr(config, n, v)


def test_notify_hook() -> None:
    """NOTIFY_CMD: off by default, message passed as one argv entry, failure never raises."""
    print("notify (optional NOTIFY_CMD hook):")
    import tempfile
    from pipeline import notify as notify_mod

    saved = config.NOTIFY_CMD
    try:
        config.NOTIFY_CMD = ""
        check("unset -> nothing sent", notify_mod.notify("x") is False)
        with tempfile.TemporaryDirectory() as d:
            out = pathlib.Path(d) / "got.txt"
            script = pathlib.Path(d) / "n.py"
            script.write_text("import sys, pathlib\n"
                              "pathlib.Path(sys.argv[1]).write_text(sys.argv[2])\n")
            config.NOTIFY_CMD = f"{sys.executable} {script} {out}"
            msg = "scan 2026-01-02: 3 change(s); $(touch x) `id` \"quoted\""
            check("configured command runs", notify_mod.notify(msg) is True)
            check("message arrives verbatim as a single argument", out.read_text() == msg)
            config.NOTIFY_CMD = f"{sys.executable} -c 'import sys; sys.exit(3)'"
            check("non-zero exit -> False, no raise", notify_mod.notify("x") is False)
            config.NOTIFY_CMD = str(pathlib.Path(d) / "no-such-notifier")
            check("missing command -> False, no raise", notify_mod.notify("x") is False)
    finally:
        config.NOTIFY_CMD = saved


def test_allele_balance() -> None:
    """Low-allele-fraction het flag: VEP-key mapping, threshold, and an AD round-trip."""
    print("allele_balance (low-VAF review flag):")
    import tempfile
    from unittest import mock
    from pipeline import allele_balance as ab

    # VEP trims the shared anchor base of an indel and shifts POS by one.
    check("insertion G>GT -> -/T at pos+1",
          ab.vep_key("6", 107893550, "G", "GT") == ("6", 107893551, "-", "T"))
    check("deletion AC>A -> C/- at pos+1",
          ab.vep_key("13", 20189546, "AC", "A") == ("13", 20189547, "C", "-"))
    check("SNV untouched", ab.vep_key("12", 102840493, "G", "A") == ("12", 102840493, "G", "A"))
    check("equal-length MNV untouched", ab.vep_key("1", 5, "AC", "GT") == ("1", 5, "AC", "GT"))

    check("het 3/18 is low", ab.is_low("HET", (3, 18)))
    check("het at exactly the threshold is not low", not ab.is_low("HET", (1, 4)))
    check("het 20/45 is not low", not ab.is_low("HET", (20, 45)))
    check("HOM is never flagged by this rule", not ab.is_low("HOM", (3, 18)))
    check("missing AD is not low", not ab.is_low("HET", None))
    check("label keeps a normal het as-is", ab.label("HET", (20, 45)) == "HET")
    check("label flags a low het with the counts",
          "3/18 (17%)" in ab.label("HET", (3, 18)) and "review" in ab.label("HET", (3, 18)))

    with tempfile.TemporaryDirectory() as d:
        raw = pathlib.Path(d) / "s.vcf"
        raw.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=6,length=200000000>\n"
            "##contig=<ID=12,length=200000000>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="a">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
            "6\t107893550\t.\tG\tGT\t39\tPASS\t.\tGT:AD\t0/1:15,3\n"
            "12\t102840493\t.\tG\tA\t900\tPASS\t.\tGT:AD\t0/1:33,28\n")
        subprocess.run(["bgzip", "-f", str(raw)], check=True)
        gz = pathlib.Path(str(raw) + ".gz")
        subprocess.run(["tabix", "-f", "-p", "vcf", str(gz)], check=True)
        with mock.patch.object(ab, "normalized_path", lambda s: gz):
            got = ab.fractions("S", [("6", 107893551, "-", "T"), ("12", 102840493, "G", "A"),
                                     ("12", 999, "C", "T")])
    check("indel AD found via its VEP key", got.get(("6", 107893551, "-", "T")) == (3, 18))
    check("SNV AD found", got.get(("12", 102840493, "G", "A")) == (28, 61))
    check("absent key simply missing", ("12", 999, "C", "T") not in got)


def test_json_format() -> None:
    """The JSON report format (docs/report-formats.md, "Adding a format"): the twin carries
    the front matter and title, content is structured by heading, and a markdown table is a
    list of row objects keyed by the column headers — never raw table text."""
    print("formats_json.to_json (structured sections + tables as row objects):")
    import json
    import tempfile
    from pipeline import formats_json, render, reportpaths, util

    check("importing formats_json registers the format", "json" in render.FORMATS)

    body = "\n".join([
        "# Panel — S1",
        "",
        "> demo notice",
        "- two findings",
        "",
        "## Findings",
        "",
        "| Gene | Variant | ClinVar |",
        "|------|---------|---------|",
        "| BRCA1 | 17:43000000 T>G | **Pathogenic** |",
        "| CFTR | 7:117548628 G>A | Pathogenic |",
        "",
        "## Caveats",
        "- Research-grade, not a clinical diagnosis.",
    ])
    doc = json.loads(formats_json.to_json(body, {"report": "panel", "person": "S1"},
                                          "Panel — S1"))
    check("front matter carried", doc["report"] == "panel" and doc["person"] == "S1")
    check("title carried", doc["title"] == "Panel — S1")
    check("content before the first heading is kept as the summary",
          doc["summary"] == ["> demo notice", "- two findings"] and doc["tables"] == [])
    sec = doc["sections"]
    check("sections are an ordered list with their level",
          [(x["heading"], x["level"]) for x in sec] == [("Findings", 2), ("Caveats", 2)])
    tbl = sec[0]["tables"][0]
    check("table headers carried", tbl["headers"] == ["Gene", "Variant", "ClinVar"])
    check("a table is row objects keyed by the column headers, not raw text",
          tbl["rows"] == [{"Gene": "BRCA1", "Variant": "17:43000000 T>G",
                           "ClinVar": "**Pathogenic**"},
                          {"Gene": "CFTR", "Variant": "7:117548628 G>A",
                           "ClinVar": "Pathogenic"}])
    check("the |---| separator row is not a data row",
          all(set(r) == {"Gene", "Variant", "ClinVar"} for r in tbl["rows"]))
    check("non-table content is the section's text",
          sec[1]["text"] == ["- Research-grade, not a clinical diagnosis."]
          and sec[1]["tables"] == [])

    # A table is structured wherever it stands. Some reports put one before any heading
    # (HLA, traits, the family panel) and the polygenic report groups its tables under ###
    # headings with no ## above them; left as raw lines, a third of the reports would be
    # unparsed. A repeated heading must not merge or overwrite.
    odd = json.loads(formats_json.to_json("\n".join([
        "# T", "", "- lead", "", "| A | B |", "|---|---|", "| 1 | 2 |", "",
        "### Group one", "", "| Trait | Pct |", "|---|---:|", "| x | 91.0% |", "",
        "## Notes", "- n", "", "### Same", "- first", "", "### Same", "- second",
    ]), {"title": "stale quoted title", "report": "prs"}, "T"))
    check("a table before any heading is structured, not left in the summary",
          odd["summary"] == ["- lead"]
          and odd["tables"] == [{"headers": ["A", "B"], "rows": [{"A": "1", "B": "2"}]}])
    check("a table under a ### heading is structured too",
          odd["sections"][0] == {"heading": "Group one", "level": 3, "text": [],
                                 "tables": [{"headers": ["Trait", "Pct"],
                                             "rows": [{"Trait": "x", "Pct": "91.0%"}]}]})
    check("a repeated heading stays two sections, in order",
          [(x["heading"], x["text"]) for x in odd["sections"][2:]]
          == [("Same", ["- first"]), ("Same", ["- second"])])
    check("the title argument wins over a title in the front matter", odd["title"] == "T")
    fields, _body = render.split_front_matter(
        '---\ntitle: "Panel \\"x\\" — S1"\nreport: panel\n---\n\n# P\n')
    check("front matter hands back the title unquoted",
          fields == {"title": 'Panel "x" — S1', "report": "panel"})

    root = pathlib.Path(tempfile.mkdtemp())
    saved = {n: getattr(config, n) for n in ("REPORTS_DIR", "FAMILY_FILE")}
    try:
        config.REPORTS_DIR = root / "reports"
        config.FAMILY_FILE = root / "family.tsv"   # absent: the id is its own person
        md = util.write_report(reportpaths.genome_report("prs", "S1"),
                               ["# PRS — S1", "", "- a result"])
        out = root / "reports" / "json" / "genome" / "S1" / "prs_S1.json"
        check("write_report renders the json twin at reports/json/<same path>", out.exists())
        doc = json.loads(out.read_text())
        check("the twin carries the front matter + the enforced closing section",
              doc["report"] == "prs" and doc["scope"] == "genome"
              and doc["person"] == "S1" and doc["callset"] == "S1"
              and [x["heading"] for x in doc["sections"]] == ["Caveats"])
        check("the markdown and HTML are written as before",
              md.exists()
              and (root / "reports" / "html" / "genome" / "S1" / "prs_S1.html").exists())
    finally:
        for n, v in saved.items():
            setattr(config, n, v)


def test_json_consent_revocation() -> None:
    """Revoking consent deletes the report from every registered format's tree — the JSON
    twin discloses exactly what the markdown did (docs/report-formats.md, "Consent")."""
    print("consent revocation removes the JSON rendering too:")
    import tempfile
    from pipeline import consent, formats_json, reportpaths, util  # noqa: F401  (registers json)

    root = pathlib.Path(tempfile.mkdtemp())
    saved = {n: getattr(config, n) for n in ("CONSENT_FILE", "FAMILY_FILE",
                                             "SAMPLES_TSV", "REPEATS_DIR", "REPORTS_DIR")}
    try:
        config.CONSENT_FILE = root / "consent.tsv"
        config.FAMILY_FILE = root / "family.tsv"
        config.SAMPLES_TSV = root / "samples.tsv"
        config.REPEATS_DIR = root / "repeats"
        config.REPORTS_DIR = root / "reports"
        config.FAMILY_FILE.write_text(
            "id\tdisplay_name\tsex\tfather\tmother\tpartner\nJan\tJan\tmale\t\t\t\n")
        config.SAMPLES_TSV.write_text("sample_id\nJan\n")
        cat = consent.CATEGORY_ADULT_ONSET_UNTREATABLE

        md = util.write_report(reportpaths.genome_report("repeats", "Jan"),
                               ["# Repeats — Jan", "", "- HTT 41 repeats"])
        twin = root / "reports" / "json" / "genome" / "Jan" / "repeats_Jan.json"
        check("the repeats report has a json twin while consent stands", twin.exists())
        consent.set_consent("Jan", cat, False)
        check("revoking removes the markdown report", not md.exists())
        check("...and its JSON rendering — a rendering left behind keeps disclosing",
              not twin.exists())
    finally:
        for n, v in saved.items():
            setattr(config, n, v)


if __name__ == "__main__":
    test_panel_status()
    test_normalize_force_reindex()
    test_load_rejects_identifier_rows()
    test_ingest_build_and_sample_guards()
    test_family_risk_rows()
    test_ar_status_index_severity()
    test_disease_status_per_gene()
    test_consent_revocation_and_per_person_gate()
    test_prs_panel_stamp_and_rescoring()
    test_scatter_tasks()
    test_read_platform_detection()
    test_fastq_naming_and_lanes()
    test_share_rebuilds_headers()
    test_diff_classify()
    test_incremental_scan()
    test_scan_completion_marker()
    test_scan_baseline_and_skipped_samples()
    test_launchd_plist_has_tool_path()
    test_gitignore_covers_data_root()
    test_annotate_retry_safety()
    test_pinned_vep_version_and_snapshot_lock()
    test_run_vep_atomicity()
    test_incremental_clinvar_delta()
    test_replacement_atomicity()
    test_cache_keys_are_content_based()
    test_validate_classify_missing_truth()
    test_validate_outcome()
    test_render()
    test_forcecall_decode()
    test_forcecall_sites()
    test_forcecall_sites_follow_the_scored_panel()
    test_acmg_constraint_pvs1()
    test_gnomad_constraint()
    test_spliceai()
    test_traits_nocall()
    test_faidx_batching()
    test_faidx_bad_region_does_not_truncate_batch()
    test_reference_freqs_partial_cache()
    test_reference_freqs_failed_query_not_cached()
    test_ref_effect_supplement_sources()
    test_kg_ref_pgen_names_its_outputs_correctly()
    test_write_score_file_drops_duplicate_coordinates()
    test_reference_scores_parallel_and_targets()
    test_short_line_guards()
    test_prs_genotype_query_failure_and_fallback()
    test_prs_partial_nocall_and_duplicate_coords()
    test_phasing_windows()
    test_phasing_verdicts()
    test_dashboard_discovery_and_layout()
    test_report_dirs_wired()
    test_serve_extension_gate()
    test_hla_uses_forcecall()
    test_prs_genomewide_percentile()
    test_mito_heteroplasmy()
    test_hla_typing()
    test_enformer()
    test_pedigree_relationships()
    test_pedigree_callsets()
    test_inheritance_origin()
    test_pedigree_ped_generation()
    test_family_migrate_plan()
    test_family_migrate_manifests()
    test_segregation()
    test_segregation_indel_matching()
    test_prs_ref_effect_supplement()
    test_sample_vcf_fallback()
    test_denovo()
    test_consent()
    test_consent_is_per_person()
    test_ingest_does_not_add_callsets_as_people()
    test_penetrance()
    test_low_severity_verdict_fixes()
    test_triage_model_factor_one_sided()
    test_diff_vep_cache_guard()
    test_clinsig_compound()
    test_diff_change_gate()
    test_validate_empty_run()
    test_demo_trio_fixtures()
    test_cmd_demo()
    test_demo_notice_and_callouts()
    test_prs_guide_lists_the_panel()
    test_perspective_switch()
    test_report_formats_and_migration()
    test_notify_hook()
    test_allele_balance()
    test_json_format()
    test_json_consent_revocation()
    print()
    if _failures:
        print(f"{len(_failures)} test(s) FAILED: {_failures}")
        sys.exit(1)
    print("ALL TESTS PASSED")
