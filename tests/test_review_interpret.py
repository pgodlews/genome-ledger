#!/usr/bin/env python3
"""Regression tests for the interpret, calibration and triage review fixes.

    uv run python tests/test_review_interpret.py

Offline: a fake client stands in for the AlphaGenome SDK and a fake curl for the download.
All people, genes, positions and scores are invented. Exits non-zero on any failure.
"""

from __future__ import annotations

import io
import json
import os
import pathlib
import stat
import subprocess
import sys
import tempfile
import traceback
import types
import urllib.error

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import test_units  # noqa: E402  (also points the reports tree at a tempdir)
from test_units import check  # noqa: E402

import duckdb  # noqa: E402

from pipeline import acmg, config, consent, enformer as em, interpret as it  # noqa: E402
from pipeline import alphagenome_client as ag, load, provenance, spliceai, triage  # noqa: E402

GENE = "INVG1"
COLS = ",".join("?" * 16)
_NAMES = ("INTERPRET_DIR", "ALPHAGENOME_LONG_FILE", "ALPHAGENOME_EGRESS_FILE",
          "ALPHAGENOME_EPOCH_FILE", "ENFORMER_SCORES_FILE", "ENFORMER_CALIBRATION_FILE",
          "DUCKDB_DIR", "DUCKDB_FILE", "LOGS_DIR", "INTERPRET_MODEL", "ALPHAGENOME_MAX_CALLS",
          "INTERPRET_CALIBRATION_N", "INTERPRET_CACHE_DIR", "INTERPRET_REMOTE_DIR",
          "ENFORMER_TARGETS_URL")


class Sandbox:
    """Point every path these tests touch at a fresh tempdir, and restore them afterwards."""

    def __enter__(self):
        self.saved = {n: getattr(config, n) for n in _NAMES}
        self.hooks = (ag._create_client, ag.time, consent.has_consented, spliceai.lookup,
                      spliceai.available, acmg.panel_context)
        self.key = os.environ.get("ALPHAGENOME_API_KEY")
        self.root = root = pathlib.Path(tempfile.mkdtemp(prefix="review_interpret_"))
        idir = root / "interpret"
        idir.mkdir()
        config.INTERPRET_DIR = idir
        config.ALPHAGENOME_LONG_FILE = idir / "long.tsv"
        config.ALPHAGENOME_EGRESS_FILE = idir / "egress.tsv"
        config.ALPHAGENOME_EPOCH_FILE = idir / "epoch.json"
        config.ENFORMER_SCORES_FILE = idir / "scores.tsv"
        config.ENFORMER_CALIBRATION_FILE = idir / "calibration.json"
        config.DUCKDB_DIR = root / "duckdb"
        config.DUCKDB_DIR.mkdir()
        config.DUCKDB_FILE = config.DUCKDB_DIR / "test.duckdb"
        config.LOGS_DIR = root / "logs"
        config.LOGS_DIR.mkdir()
        config.INTERPRET_MODEL = "alphagenome2"
        os.environ["ALPHAGENOME_API_KEY"] = "invented-test-key"
        ag.time = types.SimpleNamespace(sleep=lambda s: None)
        consent.has_consented = lambda *a: True
        ag._cache.update(sig=None, rows=[], by_key={})
        em._table, em._calibration = None, None
        return self

    def __exit__(self, *exc):
        for n, v in self.saved.items():
            setattr(config, n, v)
        (ag._create_client, ag.time, consent.has_consented, spliceai.lookup,
         spliceai.available, acmg.panel_context) = self.hooks
        if self.key is None:
            os.environ.pop("ALPHAGENOME_API_KEY", None)
        else:
            os.environ["ALPHAGENOME_API_KEY"] = self.key
        ag._cache.update(sig=None, rows=[], by_key={})
        em._table, em._calibration = None, None


def _row(snap, sample, chrom, pos, ref, alt, csq, af=0.3, sig=None, rev=None, gene=GENE):
    return (snap, sample, chrom, pos, ref, alt, gene, csq, sig, rev, None, None, af, None, None,
            "HET")


def _seed_common(snap, base, n_snv, n_indel, sample="Adam"):
    con = load.connect()
    csq = config.INTERPRET_CONSEQUENCES[0]
    rows = [_row(snap, sample, "3", base + i * 7, "A", "G", csq) for i in range(n_snv)]
    rows += [_row(snap, sample, "4", base + i * 7, "A", "AT", csq) for i in range(n_indel)]
    con.executemany(f"INSERT INTO variants VALUES ({COLS})", rows)
    con.close()


class Client:
    """Stands in for the SDK adapter; `down_after` makes it fail from that call onward."""

    class Unavailable(Exception):
        def code(self):
            return "StatusCode.UNAVAILABLE"

    def __init__(self, down_after=None):
        self.calls, self.down_after = 0, down_after

    def __call__(self, api_key):
        return self

    def score(self, v, condition):
        self.calls += 1
        if self.down_after is not None and self.calls > self.down_after:
            raise Client.Unavailable()
        value = 0.1 + (int(v["pos"]) % 97) / 10
        return [{"output_type": "OutputType.DNASE", "gene_name": None, "track_name": "d1",
                 "raw_score": value}]


def _calibrate(snap="snapA"):
    """(exited non-zero, parsed calibration file or None)."""
    failed = False
    try:
        it.calibrate(snap)
    except SystemExit:
        failed = True
    em._calibration = None
    path = config.ENFORMER_CALIBRATION_FILE
    return failed, json.loads(path.read_text()) if path.exists() else None


def test_splice_label_lookups() -> None:
    print("#8 splice_label does no per-row subprocess or file read:")
    with Sandbox():
        calls = []
        spliceai.lookup = lambda *a: calls.append(a) or None
        config.ALPHAGENOME_EPOCH_FILE.write_text(json.dumps({"id": "e1"}))
        cond = ag.conditions_for(GENE)[0]["condition"]
        ag._write_tsv(config.ALPHAGENOME_LONG_FILE, [
            {"chrom": "3", "pos": 700100, "ref": "A", "alt": "G", "condition": cond,
             "scorer": "SPLICE_SITES", "gene": GENE, "track": "t", "raw_score": 0.42,
             "epoch": "e1"}], ag.LONG_COLS)

        opens = {"epoch": 0, "tissue": 0}
        counting = [True]

        def hook(event, args):
            if event == "open" and counting[0]:
                name = str(args[0])
                if name == str(config.ALPHAGENOME_EPOCH_FILE):
                    opens["epoch"] += 1
                elif name.endswith("gene_tissue.tsv"):
                    opens["tissue"] += 1

        sys.addaudithook(hook)
        ag.current_rows()
        opens.update(epoch=0, tissue=0)
        for i in range(50):
            ag.splice_label("3", 800000 + i, "A", "G", GENE, "intron_variant")
        counting[0] = False
        check("no SpliceAI query for a variant without AlphaGenome rows", not calls)
        check("no tissue-table read for a variant without rows", opens["tissue"] == 0)
        check("no epoch-file read per variant", opens["epoch"] == 0)

        got = ag.splice_label("3", 700100, "A", "G", GENE, "intron_variant")
        check("a cached AlphaGenome splice score is still shown", "SPLICE_SITES 0.42" in got)
        spliceai.lookup = lambda *a: (0.1, "donor loss", GENE)
        check("a SpliceAI record still supersedes it",
              ag.splice_label("3", 700100, "A", "G", GENE, "intron_variant") == "—")

        ag._atomic(config.ALPHAGENOME_EPOCH_FILE, json.dumps({"id": "e2"}))
        check("a replaced epoch file is picked up", ag.epoch().get("id") == "e2")
        config.ALPHAGENOME_EPOCH_FILE.unlink()
        check("a removed epoch file reads as no epoch", ag.epoch() == {})


def test_error_code_forms() -> None:
    print("#28 Session.score handles a non-callable .code:")
    with Sandbox():
        class Raises:
            def __init__(self, exc):
                self.exc = exc

            def score(self, v, condition):
                raise self.exc

        secret = "invented-credential-text"
        exc = urllib.error.HTTPError("https://example.invalid/", 403, f"Forbidden {secret}", {},
                                     io.BytesIO())
        ag._create_client = lambda key: Raises(exc)
        session = ag.Session("invented-test-key", max_calls=5)
        kind, text = None, ""
        try:
            session.score({"chrom": "3", "pos": 700000, "ref": "A", "alt": "G"}, "c")
        except Exception as raised:  # noqa: BLE001
            kind, text = type(raised), "".join(traceback.format_exception(raised))
        check("an HTTPError becomes Stopped", kind is ag.Stopped)
        check("its status is reported", "403" in text)
        check("the raw exception text is not printed", secret not in text)


def test_calibration_never_partial() -> None:
    print("#29 an incomplete calibration is not published:")
    with Sandbox():
        config.INTERPRET_CALIBRATION_N = 120
        _seed_common("snapA", 1_000_000, 120, 120)
        client = Client()
        ag._create_client = client

        config.ALPHAGENOME_MAX_CALLS = 200       # 5 sentinels + 120 SNV + 75 indel
        failed, data = _calibrate()
        check("a run stopped by the cap exits non-zero", failed)
        check("and writes no calibration", data is None)
        check("and leaves bands uncalibrated", not em.calibration_available())

        failed, data = _calibrate()
        check("the rerun resumes and completes", not failed and data is not None)
        check("with both strata full", bool(data) and (data["n_snv"], data["n_indel"]) == (120, 120))
        check("stamped with the current epoch", bool(data) and data["epoch"] == ag.epoch()["id"])
        good = config.ENFORMER_CALIBRATION_FILE.read_bytes()

        # A later run on a snapshot whose variants are not cached, with the service down.
        _seed_common("snapB", 50_000_000, 120, 120)
        ag._create_client = Client(down_after=0)
        failed, _ = _calibrate("snapB")
        check("an outage exits non-zero", failed)
        check("and keeps the good null byte for byte",
              config.ENFORMER_CALIBRATION_FILE.read_bytes() == good)

        # The service dies part-way through the same kind of run.
        ag._create_client = Client(down_after=40)
        failed, _ = _calibrate("snapB")
        check("a mid-run outage exits non-zero", failed)
        check("and keeps the good null byte for byte",
              config.ENFORMER_CALIBRATION_FILE.read_bytes() == good)
        check("the kept null still grades", em.calibration_available())

    with Sandbox():
        config.INTERPRET_CALIBRATION_N = 120
        _seed_common("snapC", 1_000_000, 120, 3)     # a stratum too small to resolve p99
        ag._create_client = Client()
        failed, data = _calibrate("snapC")
        check("a tiny stratum is refused", failed and data is None)

    with Sandbox():
        config.INTERPRET_CALIBRATION_N = 120
        _seed_common("snapD", 1_000_000, 120, 120)
        ag._create_client = Client()
        consent.has_consented = lambda *a: False     # nobody eligible: nothing is scored
        failed, data = _calibrate("snapD")
        check("no eligible variants publishes nothing", failed and data is None)


def test_local_calibration_never_partial() -> None:
    print("#29 the local GPU path has the same guard:")
    with Sandbox():
        config.INTERPRET_MODEL = "borzoi"
        config.INTERPRET_CALIBRATION_N = 120
        _seed_common("snapE", 1_000_000, 120, 120)
        good = json.dumps({"snv": [1.0, 2.0], "indel": [3.0, 4.0]})
        config.ENFORMER_CALIBRATION_FILE.write_text(good)
        ref = config.REF_FASTA
        config.REF_FASTA = config.INTERPRET_DIR / "invented.fa"
        config.REF_FASTA.write_text(">3\nA\n")
        saved = it.require_tools, it._write_job, it._run_remote
        variants = it._select_common("snapE")
        placed = {f"v{i}": v for i, v in enumerate(variants)}

        def worker(job, out, keep):
            rows = ["id\tdelta_max\tl2_center\ttop_track\ttop_delta\ttop_desc"]
            rows += [f"{vid}\t{1 + n / 10}\t0.1\t0\t0.1\tDNASE:x"
                     for n, vid in enumerate(list(placed)[:keep])]
            pathlib.Path(out).write_text("\n".join(rows) + "\n")

        it.require_tools = lambda *a: None
        it._write_job = lambda v, job: placed
        try:
            it._run_remote = lambda job, out: worker(job, out, 50)   # short worker output
            failed = False
            try:
                it.calibrate("snapE", force=True)
            except SystemExit:
                failed = True
            check("a short worker result exits non-zero", failed)
            check("and keeps the existing null",
                  config.ENFORMER_CALIBRATION_FILE.read_text() == good)

            it._run_remote = lambda job, out: worker(job, out, len(placed))
            it.calibrate("snapE", force=True)
            data = json.loads(config.ENFORMER_CALIBRATION_FILE.read_text())
            check("a full worker result is published",
                  (data["n_snv"], data["n_indel"]) == (120, 120))
        finally:
            it.require_tools, it._write_job, it._run_remote = saved
            config.REF_FASTA = ref


def _load_worker():
    """interpret_worker imports torch at module level; stub it where it is not installed."""
    stubbed = "torch" not in sys.modules
    try:
        import torch  # noqa: F401
        stubbed = False
    except ImportError:
        torch_mod, nn, func = (types.ModuleType(n) for n in
                               ("torch", "torch.nn", "torch.nn.functional"))
        torch_mod.nn, nn.functional, torch_mod.Tensor = nn, func, object
        sys.modules.update({"torch": torch_mod, "torch.nn": nn, "torch.nn.functional": func})
    from pipeline import interpret_worker
    if stubbed:
        for name in ("torch", "torch.nn", "torch.nn.functional"):
            sys.modules.pop(name, None)
    return interpret_worker


def test_strand_pair_validation() -> None:
    print("#30 a truncated targets file is rejected:")
    worker = _load_worker()
    n = 40
    pair = [i + 1 if i % 2 == 0 else i - 1 for i in range(n)]
    head = "\tidentifier\tstrand_pair\tdescription"
    lines = [head] + [f"{i}\tINV{i}\t{pair[i]}\tinvented track {i}" for i in range(n)]
    with tempfile.TemporaryDirectory() as d:
        path = pathlib.Path(d) / "targets.txt"

        def load(text):
            path.write_text(text)
            return worker._load_strand_pair(str(path), n)

        check("an intact table loads", load("\n".join(lines) + "\n") == pair)
        check("a file cut at a line boundary is refused",
              load("\n".join(lines[:n // 2 + 1]) + "\n") is None)
        check("a file cut inside the description is refused",
              load("\n".join(lines[:n // 2]) + f"\n{n // 2}\tINV\t{pair[n // 2]}\tinv") is None)
        damaged = lines[:-1] + [f"{n - 1}\tINV{n - 1}\t{str(pair[n - 1])[:1]}\tinvented"]
        check("a partner number cut in half is refused", load("\n".join(damaged) + "\n") is None)
        identity = [head] + [f"{i}\tINV{i}\t{i}\tinvented" for i in range(n)]
        check("an all-unstranded table still loads",
              load("\n".join(identity) + "\n") == list(range(n)))


def test_targets_download_atomic() -> None:
    print("#30 an interrupted targets download leaves no file:")
    with Sandbox() as box:
        cache = box.root / "cache"
        cache.mkdir()
        fake_bin = box.root / "bin"
        fake_bin.mkdir()
        curl = fake_bin / "curl"
        # Writes a partial body to the -o path, then fails like a dropped transfer.
        curl.write_text('#!/bin/sh\nprintf "partial" > "$3"\nexit "${FAKE_CURL_RC:-18}"\n')
        curl.chmod(curl.stat().st_mode | stat.S_IEXEC)
        config.INTERPRET_CACHE_DIR, config.INTERPRET_REMOTE_DIR = str(cache), str(box.root / "jobs")
        config.ENFORMER_TARGETS_URL = "https://example.invalid/targets.txt"
        config.INTERPRET_MODEL = "borzoi"
        job, out = box.root / "job.tsv.gz", box.root / "out.tsv"
        commands = []
        saved = it.run, it.require_tools
        it.run = lambda argv, **kw: commands.append(argv)
        it.require_tools = lambda *a: None
        try:
            it._run_remote(job, out)
        finally:
            it.run, it.require_tools = saved
        shell = next(c[2] for c in commands if c[0] == "ssh" and "curl" in c[2])
        target = cache / "targets_borzoi.txt"
        env = {**os.environ, "PATH": f"{fake_bin}:{os.environ['PATH']}"}
        rc = subprocess.run(["sh", "-c", shell], env=env).returncode
        check("the step stays best-effort", rc == 0)
        check("a failed download leaves no targets file", not target.exists())
        check("nor a stray temp file", not list(cache.glob("*.part")))
        rc = subprocess.run(["sh", "-c", shell], env={**env, "FAKE_CURL_RC": "0"}).returncode
        check("a good download lands at the real name", rc == 0 and target.exists())


def test_missense_spliceai() -> None:
    print("#36 triage and provenance see SpliceAI on a missense:")
    with Sandbox():
        spliceai.available = lambda: True
        spliceai.lookup = lambda *a: (0.90, "acceptor gain", GENE)
        am, sa, en = triage.model_lookups("3", 700200, "A", "G", "missense_variant",
                                          False, True, False)
        check("model_lookups returns SpliceAI for a missense", bool(sa) and sa[0] == 0.90)
        _, sa, _ = triage.model_lookups("3", 700200, "A", "G", "protein_altering_variant",
                                        False, True, False)
        check("and for a protein-altering variant", bool(sa))
        _, sa, _ = triage.model_lookups("3", 700200, "A", "G", "synonymous_variant",
                                        False, True, False)
        check("but not for a synonymous variant", sa is None)

        acmg.panel_context = lambda: ({GENE}, {GENE: {"AD"}})
        con = load.connect()
        con.execute(f"INSERT INTO variants VALUES ({COLS})",
                    _row("snapF", "Maria", "3", 700200, "A", "G", "missense_variant", af=None,
                         sig="Uncertain_significance"))
        con.close()
        provenance.record_provenance("snapF")
        con = load.connect()
        rows = con.execute("SELECT value_summary, used_in_acmg FROM annotation_provenance "
                           "WHERE snapshot_id='snapF' AND source='spliceai'").fetchall()
        con.close()
        check("provenance records the SpliceAI evidence", len(rows) == 1)
        check("marked as used by ACMG", bool(rows) and rows[0][1] is True)


def test_label_needs_alleles() -> None:
    print("#44 enformer.label needs the alleles of a scored variant:")
    with Sandbox():
        em._calibration = {"snv": [i / 100 for i in range(1, 401)],
                           "indel": [float(i) for i in range(1, 401)]}
        en = (3.99, 0.0, 0, "CAGE:invented")
        raised = False
        try:
            em.label(en)
        except ValueError:
            raised = True
        check("a scored variant without alleles is refused", raised)
        check("no score still renders a dash", em.label(None) == "—")
        check("an SNV grades against the SNV null", "p99.8 high" in em.label(en, "A", "G"))


TESTS = (test_splice_label_lookups, test_error_code_forms, test_calibration_never_partial,
         test_local_calibration_never_partial, test_strand_pair_validation,
         test_targets_download_atomic, test_missense_spliceai, test_label_needs_alleles)


if __name__ == "__main__":
    for test in TESTS:
        try:
            test()
        except Exception as exc:  # noqa: BLE001 — one crash must not hide the other tests
            traceback.print_exc()
            check(f"{test.__name__} ran to completion ({type(exc).__name__})", False)
    print()
    if test_units._failures:
        print(f"{len(test_units._failures)} check(s) FAILED: {test_units._failures}")
        sys.exit(1)
    print("ALL TESTS PASSED")
