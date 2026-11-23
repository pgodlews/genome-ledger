#!/usr/bin/env python3
"""Regression tests for data-operation bugs: family-migrate, share, normalize, setup, config.

    GENOMES_ROOT=/some/temp/dir uv run python tests/test_review_dataops.py

Every test works in its own temporary root with invented people and invented sequence.
"""

from __future__ import annotations

import gzip
import logging
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import traceback
import types

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from test_units import _failures, check  # noqa: E402

import duckdb  # noqa: E402

from pipeline import config  # noqa: E402

BASE = pathlib.Path(tempfile.mkdtemp(prefix="genome_review_dataops_")).resolve()
config.set_genomes_root(BASE / "start")    # nothing below can reach a real root
_n = 0


def fresh_root() -> pathlib.Path:
    global _n
    _n += 1
    root = BASE / f"r{_n}"
    config.set_genomes_root(root)
    config.REPORT_DIRS = (config.REPORTS_DIR,)
    return root


def touch(path: pathlib.Path, text: str = "x") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def have(*tools: str) -> bool:
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        print(f"  SKIP  {', '.join(missing)} not on PATH")
    return not missing


def seed_family(ids: list[str], people: list[str] | None = None) -> None:
    for s in ids:
        touch(config.ALIGNED_DIR / f"{s}.cram", f"cram-of-{s}")
        touch(config.NORMALIZED_DIR / f"{s}.norm.vcf.gz", f"vcf-of-{s}")
        touch(config.FORCECALL_DIR / f"{s}.{config.BUILD}.vcf.gz", f"forced-{s}")
        touch(config.PGX_DIR / s / f"pgx_{s}.md", f"pgx-of-{s}")
    touch(config.SAMPLES_TSV, "sample_id\trelationship\tsex\tbuild\tsource\traw_path\n" + "".join(
        f"{s}\tself\tmale\tGRCh38\t{config.RAW_DIR}/{s}.vcf.gz\t{config.RAW_DIR}/{s}.vcf.gz\n"
        for s in ids))
    touch(config.FAMILY_FILE, "id\tdisplay_name\tsex\tfather\tmother\tpartner\n" + "".join(
        f"{p}\t{p}\tmale\t\t\t\n" for p in (people if people is not None else ids)))
    config.DUCKDB_DIR.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(config.DUCKDB_FILE))
    con.execute("CREATE TABLE variants(snapshot_id VARCHAR, sample_id VARCHAR, pos BIGINT)")
    for i, s in enumerate(ids):
        con.execute("INSERT INTO variants VALUES ('snap', ?, ?)", [s, 100 + i])
    con.close()


def db_ids() -> dict[str, int]:
    con = duckdb.connect(str(config.DUCKDB_FILE), read_only=True)
    try:
        return dict(con.execute(
            "SELECT sample_id, min(pos) FROM variants GROUP BY 1").fetchall())
    finally:
        con.close()


def snapshot(root: pathlib.Path) -> dict[str, bytes]:
    return {str(p.relative_to(root)): p.read_bytes()
            for p in sorted(root.rglob("*")) if p.is_file()}


def test_migrate_refuses_unsafe_mappings() -> None:
    """A mapping that could merge people must stop before the first change, dry-run or not."""
    print("family_migrate: unsafe mappings are refused before anything changes:")
    from pipeline import family_migrate as fm

    cases = {
        "chain": (["Adam", "Jan"], {"Adam": "Jan", "Jan": "Maria"}),
        "swap": (["Adam", "Jan"], {"Adam": "Jan", "Jan": "Adam"}),
        "target already exists": (["Adam", "Jan"], {"Adam": "Jan"}),
        "two ids onto one": (["Adam", "Jan"], {"Adam": "Lena", "Jan": "Lena"}),
        "person without their callset": (["Adam", "Adam_vendor"], {"Adam": "Jan"}),
        "unusable target": (["Adam"], {"Adam": "../Jan"}),
    }
    for label, (ids, mapping) in cases.items():
        root = fresh_root()
        seed_family(ids, people=[i for i in ids if "_" not in i])
        before = snapshot(root)
        try:
            fm.migrate(apply=True, mapping=mapping)
            refused = False
        except SystemExit as e:
            refused = "refused" in str(e)
        check(f"{label}: refused", refused)
        check(f"{label}: nothing changed on disk", snapshot(root) == before)

    root = fresh_root()
    seed_family(["Adam"])
    touch(config.PGX_DIR / "Jan" / "stale.md")     # a directory no manifest knows about
    before = snapshot(root)
    try:
        fm.migrate(apply=True, mapping={"Adam": "Jan"})
        refused = False
    except SystemExit as e:
        refused = "existing file" in str(e)
    check("a rename onto an existing directory is refused", refused)
    check("...and nothing changed", snapshot(root) == before)


def test_migrate_keeps_callsets_with_their_ids() -> None:
    print("family_migrate: callset files follow their own ids:")
    from pipeline import family_migrate as fm

    fresh_root()
    seed_family(["Adam", "Adam_vendor", "Adam_Nowak"], people=["Adam", "Adam_Nowak"])
    fm.migrate(apply=True, mapping={"Adam": "Jan", "Adam_vendor": "Jan_vendor"})
    crams = sorted(p.name for p in config.ALIGNED_DIR.iterdir())
    check("each callset's file got its own new name",
          crams == ["Adam_Nowak.cram", "Jan.cram", "Jan_vendor.cram"])
    check("file contents stayed with the right callset",
          (config.ALIGNED_DIR / "Jan_vendor.cram").read_text() == "cram-of-Adam_vendor")
    check("the database agrees with the files",
          set(db_ids()) == {"Jan", "Jan_vendor", "Adam_Nowak"})
    rows = [ln.split("\t") for ln in config.SAMPLES_TSV.read_text().splitlines()[1:]]
    check("a different person sharing the prefix is untouched",
          any(r[0] == "Adam_Nowak" and r[5].endswith("/Adam_Nowak.vcf.gz") for r in rows))

    fresh_root()
    seed_family(["Lena"])
    mapping = {"Lena": "Lena_vendor"}
    fm.migrate(apply=True, mapping=mapping)
    first = sorted(p.name for p in config.ALIGNED_DIR.iterdir())
    fm.migrate(apply=True, mapping=mapping)
    check("re-running a suffixing rename changes nothing",
          first == ["Lena_vendor.cram"]
          and sorted(p.name for p in config.ALIGNED_DIR.iterdir()) == first)


def test_migrate_moves_stage_outputs() -> None:
    print("family_migrate: per-sample stage outputs follow the rename:")
    from pipeline import family_migrate as fm, forcecall

    fresh_root()
    seed_family(["Adam"])
    stage_dirs = (config.ANCESTRY_DIR, config.TRAITS_DIR, config.HLA_DIR, config.HLA_TYPE_DIR,
                  config.PRS_DIR, config.REPEATS_DIR, config.SV_DIR, config.CALLABLE_DIR,
                  config.PHASE_DIR, config.MITO_DIR, config.DENOVO_DIR, config.SEGREGATION_DIR)
    for d in stage_dirs:
        touch(d / "Adam" / "out_Adam.txt")
    touch(config.SHARE_DIR / "Adam" / "Adam.GRCh38.vcf.gz")
    fm.migrate(apply=True, mapping={"Adam": "Jan"})
    check("force-called VCF is still found under the new id", forcecall.available("Jan"))
    check("PGx output moved", (config.PGX_DIR / "Jan" / "pgx_Jan.md").exists())
    check("every other stage directory moved",
          all((d / "Jan" / "out_Jan.txt").exists() and not (d / "Adam").exists()
              for d in stage_dirs))
    check("copies already shared are left alone",
          (config.SHARE_DIR / "Adam" / "Adam.GRCh38.vcf.gz").exists())


def test_migrate_paths_and_atomic_db() -> None:
    print("family_migrate: path cells and the database update:")
    from pipeline import family_migrate as fm

    root = fresh_root()
    outside = BASE / "Adam" / "Adam.vcf.gz"            # exists outside the genomes root
    touch(outside)
    seed_family(["Adam"])
    config.SAMPLES_TSV.write_text(
        "sample_id\tsource\traw_path\n"
        f"Adam\t{outside}\t{root}/raw/Adam/Adam.vcf.gz\n")
    fm.migrate(apply=True, mapping={"Adam": "Jan"})
    row = config.SAMPLES_TSV.read_text().splitlines()[1].split("\t")
    check("a file outside the tree keeps its path", row[1] == str(outside))
    check("a path inside the tree follows the rename", row[2] == f"{root}/raw/Jan/Jan.vcf.gz")

    fresh_root()
    seed_family(["Adam"])
    con = duckdb.connect(str(config.DUCKDB_FILE))
    con.execute("CREATE TABLE scan_inputs(sample_id VARCHAR CHECK (sample_id <> 'Jan'))")
    con.execute("INSERT INTO scan_inputs VALUES ('Adam')")
    con.close()
    try:
        fm._apply_duckdb({"Adam": "Jan"})
    except duckdb.Error:
        pass
    check("a failed table rolls every table back", set(db_ids()) == {"Adam"})


def test_migrate_follows_set_genomes_root() -> None:
    print("family_migrate: manifests are read from the current root:")
    from pipeline import family_migrate as fm

    first = fresh_root()
    touch(config.SAMPLES_TSV, "sample_id\tsource\traw_path\nAdam\tx\ty\n")
    first_tsv = config.SAMPLES_TSV
    fresh_root()
    seed_family(["Adam"])
    fm.migrate(apply=True, mapping={"Adam": "Jan"})
    check("active root migrated",
          config.SAMPLES_TSV.read_text().splitlines()[1].startswith("Jan\t"))
    check("previous root untouched",
          first_tsv.read_text().splitlines()[1].startswith("Adam\t"))
    check("no file created in the previous root", (first / "duckdb").exists() is False)


def test_set_genomes_root_is_complete() -> None:
    print("config.set_genomes_root repoints every derived path:")
    old = fresh_root()
    names = [k for k, v in vars(config).items()
             if k.isupper() and isinstance(v, pathlib.Path) and str(v).startswith(str(old))]
    new = BASE / "moved"
    config.set_genomes_root(new)
    stale = [k for k in names if not str(getattr(config, k)).startswith(str(new))]
    check(f"no path left on the old root ({stale})", not stale)
    check("the tool paths are among those checked",
          {"VEP_BIN", "HAPLOGREP_JAR", "PHARMCAT_JAR", "PHARMCAT_PIPELINE",
           "KG_PANEL_FILE"} <= set(names))


def test_share_rejects_dot_names() -> None:
    print("share: a name cannot leave SHARE_DIR:")
    from pipeline import share

    root = fresh_root()
    for name in ("..", ".", ".hidden"):
        try:
            share.share("Adam", "vcf", name)
            msg = ""
        except SystemExit as e:
            msg = str(e)
        check(f"{name!r} rejected by the name check", "must start with" in msg)
    check("no directory was created", not (root / "share").exists())


def _write_vcf(path: pathlib.Path, contigs: list[str], records: list[tuple]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["##fileformat=VCFv4.2", '##FILTER=<ID=PASS,Description="ok">',
             '##FILTER=<ID=low,Description="low">',
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">']
    lines += [f"##contig=<ID={c},length=60>" for c in contigs]
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAdam")
    lines += ["\t".join(map(str, (c, p, ".", r, a, 50, f, ".", "GT", "0/1")))
              for c, p, r, a, f in records]
    plain = path.with_suffix("")
    plain.write_text("\n".join(lines) + "\n")
    subprocess.run(["bgzip", "-f", str(plain)], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(path)], check=True)


class _Capture(logging.Handler):
    def __init__(self) -> None:
        super().__init__(logging.WARNING)
        self.messages: list[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.messages.append(record.getMessage())


def test_normalize_handles_foreign_contigs() -> None:
    print("normalize: records on contigs the reference lacks:")
    if not have("bcftools", "tabix", "bgzip", "samtools"):
        return
    from pipeline import normalize as nm

    alt = "chr1_KI270706v1_random"

    def prepare(contigs, records):
        root = fresh_root()
        config.REFS_DIR.mkdir(parents=True, exist_ok=True)
        config.REF_FASTA.write_text(">1\n" + "ACGTACGTAC" * 6 + "\n")
        subprocess.run(["samtools", "faidx", str(config.REF_FASTA)], check=True)
        raw = root / "in" / "adam.vcf.gz"
        _write_vcf(raw, contigs, records)
        touch(config.SAMPLES_TSV, f"sample_id\traw_path\nAdam\t{raw}\n")

    def positions():
        out = subprocess.run(["bcftools", "query", "-f", "%CHROM:%POS\n",
                              str(nm.normalized_path("Adam"))],
                             capture_output=True, text=True, check=True).stdout
        return out.split()

    cap = _Capture()
    logging.getLogger("genome").addHandler(cap)
    try:
        prepare(["chr1", alt], [("chr1", 5, "A", "T", "PASS"), ("chr1", 9, "A", "G", "PASS"),
                                ("chr1", 13, "G", "C", "PASS"), (alt, 7, "G", "A", "PASS"),
                                (alt, 9, "G", "A", "low")])
        try:
            nm.normalize("Adam")
            ok = True
        except BaseException:  # noqa: BLE001 — the old code died here
            ok = False
        check("a sample with an alt-contig call normalizes", ok)
        check("reference contigs are all kept", ok and positions() == ["1:5", "1:9", "1:13"])
        check("the dropped PASS record is counted in the log",
              any("dropping 1 record" in m and f"{alt}=1" in m for m in cap.messages))

        cap.messages.clear()
        prepare(["chr1", alt], [("chr1", 5, "A", "T", "PASS"), (alt, 9, "G", "A", "low")])
        nm.normalize("Adam")
        check("a foreign contig with no PASS record drops nothing",
              positions() == ["1:5"] and not cap.messages)
        check("...and the bcftools command is unchanged", nm._reference_targets(
            "Adam", str(config.GENOMES_ROOT / "in" / "adam.vcf.gz")) == [])

        prepare(["chr1", "NC_000000.1"], [("chr1", 5, "A", "T", "PASS"),
                                           ("NC_000000.1", 7, "G", "A", "PASS"),
                                           ("NC_000000.1", 9, "G", "A", "PASS")])
        try:
            nm.normalize("Adam")
            msg = ""
        except SystemExit as e:
            msg = str(e)
        except BaseException:  # noqa: BLE001
            msg = ""
        check("a file mostly on foreign contigs is refused", "naming or build mismatch" in msg)
        check("...and leaves no output", not nm.normalized_path("Adam").exists())
    finally:
        logging.getLogger("genome").removeHandler(cap)


def test_query_quotes_params() -> None:
    print("run.py query: a quote in a value stays inside the literal:")
    import run as cli

    root = fresh_root()
    config.DUCKDB_DIR.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(config.DUCKDB_FILE))
    con.execute("CREATE TABLE people(name VARCHAR)")
    con.execute("INSERT INTO people VALUES ('O''Lena'), ('Adam')")
    con.close()
    sql = root / "q.sql"
    sql.write_text("SELECT name FROM people WHERE name = '${who}';\n")

    def ask(value: str) -> list[str]:
        out = root / "out.txt"
        saved = sys.stdout
        try:
            with open(out, "w") as fh:
                sys.stdout = fh
                cli.cmd_query(types.SimpleNamespace(sqlfile=str(sql), param=[f"who={value}"]))
        except Exception:  # noqa: BLE001 — a parser error is the old behaviour
            return ["<error>"]
        finally:
            sys.stdout = saved
        return out.read_text().splitlines()[1:]

    check("an apostrophe matches literally", ask("O'Lena") == ["O'Lena"])
    check("a quoted condition is only a string", ask("x' OR '1'='1") == [])


def test_clinvar_delta_quoted_path() -> None:
    print("incremental: a quote in GENOMES_ROOT:")
    from pipeline import incremental as inc
    from pipeline.annotate import clinvar_vcf

    root = BASE / "Lena's genomes"
    config.set_genomes_root(root)
    for snap, sig in (("2026-01-01", "Benign"), ("2026-02-01", "Pathogenic")):
        path = clinvar_vcf(snap)
        path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, "wt") as fh:
            fh.write("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                     f"1\t5\t1\tA\tT\t.\t.\tCLNSIG={sig};ALLELEID=9\n")
    inc._DELTA_CACHE.clear()
    try:
        got = inc.clinvar_delta_coords("2026-01-01", "2026-02-01")
    except Exception:  # noqa: BLE001 — the old code raised a parser error
        got = None
    check("the ClinVar delta is computed", got == [("1", 5, "A", "T")])


def test_site_regions_stay_one_based() -> None:
    print("incremental.site_regions: no position 0:")
    from pipeline import incremental as inc

    got = inc.site_regions([("1", 1, "-", "A"), ("1", 8, "-", "A")])
    check("position 1 asks for no base before it", got == [("1", 1), ("1", 7), ("1", 8)])


def test_consent_short_row() -> None:
    print("consent: a hand-edited row with no decision:")
    from pipeline import consent

    fresh_root()
    cat = consent.CATEGORY_ADULT_ONSET_UNTREATABLE
    touch(config.CONSENT_FILE, f"sample_id\tcategory\tconsent\nAdam\t{cat}\tyes\nLena\t{cat}\n")
    try:
        adam, lena = consent.has_consented("Adam", cat), consent.has_consented("Lena", cat)
    except AttributeError:
        adam, lena = None, None
    check("the short row reads as not consented", lena is False)
    check("another person's grant still works", adam is True)


def test_setup_interrupted_steps() -> None:
    print("setup: an interrupted step cannot be mistaken for a finished one:")
    if not have("samtools", "bgzip", "tabix", "gzip"):
        return
    from pipeline import setup as st

    fresh_root()
    config.REFS_DIR.mkdir(parents=True, exist_ok=True)
    names = [n for n in dir(st) if n.startswith("_install_")]
    saved = {n: getattr(st, n) for n in (*names, "_system_tools", "_download")}
    real_alphamissense = st._install_alphamissense
    try:
        for n in names:
            setattr(st, n, lambda *a, **k: None)
        st._system_tools = lambda: None
        st._download = lambda url, dest, label: None

        fa_gz = config.REFS_DIR / f"{config.BUILD}.primary_assembly.fa.gz"
        full = "".join(f">{c}\n" + ("ACGTTGCA" * 7 + "\n") * 700000 for c in ("1", "2", "3"))
        blob = gzip.compress(full.encode(), compresslevel=1)
        fa_gz.write_bytes(blob[:len(blob) // 2])
        try:
            st.setup()
        except EOFError:
            pass
        check("a failed decompression leaves no reference", not config.REF_FASTA.exists())
        fa_gz.write_bytes(blob)
        st.setup()
        check("the next run builds the complete reference",
              config.REF_FASTA.read_text() == full)

        raw = config.REFS_DIR / "AlphaMissense_hg38.tsv.gz"
        body = "".join(f"chr1\t{i}\tA\tT\n" for i in range(1, 200000))
        data = gzip.compress(body.encode())
        raw.write_bytes(data[:len(data) // 2])
        try:
            real_alphamissense()
        except SystemExit:
            pass
        check("a failed recompression leaves no AlphaMissense file",
              not config.ALPHAMISSENSE_FILE.exists())
        raw.write_bytes(data)
        real_alphamissense()
        tbi = pathlib.Path(f"{config.ALPHAMISSENSE_FILE}.tbi")
        check("the next run produces the file and its index",
              config.ALPHAMISSENSE_FILE.exists() and tbi.exists())
        out = subprocess.run(["tabix", str(config.ALPHAMISSENSE_FILE), "chr1:7-7"],
                             capture_output=True, text=True, check=True).stdout
        check("the installed index works", out == "chr1\t7\tA\tT\n")
        tbi.unlink()
        raw.write_bytes(data)
        real_alphamissense()
        check("a file without its index is rebuilt", tbi.exists())
    finally:
        for n, v in saved.items():
            setattr(st, n, v)


def test_demo_logs_under_its_own_root() -> None:
    print("run.py demo: the log file stays in the demo directory:")
    import run as cli

    original = fresh_root()
    out = BASE / "demo_out"
    touch(out / "stale.txt")
    seen = {}
    saved = cli.cmd_demo
    root_logger = logging.getLogger()
    handlers = list(root_logger.handlers)
    root_logger.handlers = []

    def stub(args):
        seen["clean"] = args.clean
        return 0

    cli.cmd_demo = stub
    try:
        cli.main(["demo", "--out-dir", str(out), "--clean"])
    finally:
        cli.cmd_demo = saved
        for h in root_logger.handlers:
            h.close()
        root_logger.handlers = handlers
    check("nothing is logged under the original root", not (original / "logs").exists())
    check("the log is in the demo directory", any((out / "logs").glob("*.log")))
    check("--clean ran before the log was opened", not (out / "stale.txt").exists())
    check("and will not run again over the open log", seen.get("clean") is False)


TESTS = [
    test_migrate_refuses_unsafe_mappings,
    test_migrate_keeps_callsets_with_their_ids,
    test_migrate_moves_stage_outputs,
    test_migrate_paths_and_atomic_db,
    test_migrate_follows_set_genomes_root,
    test_set_genomes_root_is_complete,
    test_share_rejects_dot_names,
    test_normalize_handles_foreign_contigs,
    test_query_quotes_params,
    test_clinvar_delta_quoted_path,
    test_site_regions_stay_one_based,
    test_consent_short_row,
    test_setup_interrupted_steps,
    test_demo_logs_under_its_own_root,
]

if __name__ == "__main__":
    failed_tests = []
    try:
        for test in TESTS:
            before = len(_failures)
            try:
                test()
            except BaseException as exc:  # noqa: BLE001 — report it and keep going
                traceback.print_exc()
                check(f"{test.__name__} ran to completion ({type(exc).__name__})", False)
            if len(_failures) > before:
                failed_tests.append(test.__name__)
    finally:
        shutil.rmtree(BASE, ignore_errors=True)
    print()
    if _failures:
        print(f"{len(_failures)} check(s) FAILED in {len(failed_tests)} of {len(TESTS)} tests: "
              f"{failed_tests}")
        sys.exit(1)
    print(f"ALL {len(TESTS)} TESTS PASSED")
