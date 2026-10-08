"""Stage 5: annotate a normalized sample against a pinned snapshot with VEP.

Gene / consequence / gnomAD AF / SIFT / PolyPhen come from the offline VEP cache
(its version is recorded in the snapshot manifest). Clinical significance comes
from the snapshot's dated ClinVar VCF via `--custom`, so it is version-pinned by
us and therefore diffable across snapshots.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

from . import config
from .normalize import normalized_path
from .snapshot import load_manifest, snapshot_dir
from .util import append_tsv, log, run, sha256, today


def _vep_env() -> dict:
    """VEP needs its own root on PERL5LIB to find the source-installed
    Bio::DB::HTS (used for --custom VCFs and --fasta)."""
    env = dict(os.environ)
    env["PERL5LIB"] = f"{config.VEP_DIR}:{env.get('PERL5LIB', '')}".rstrip(":")
    return env

RUN_FIELDS = ["snapshot_id", "sample_id", "run_date",
              "vep_cache_version", "clinvar_file", "annotated_path"]

# VEP --tab output columns we request. ZYG (zygosity) comes from --individual all.
VEP_FIELDS = ",".join([
    "Uploaded_variation", "Location", "Allele", "SYMBOL", "Consequence",
    "gnomADe_AF", "gnomADg_AF", "SIFT", "PolyPhen", "ZYG",
    "CLINVAR_CLNSIG", "CLINVAR_CLNREVSTAT", "CLINVAR_CLNDN", "CLINVAR_ALLELEID",
])


def annotated_path(snapshot_id: str, sample: str) -> Path:
    return config.ANNOTATED_DIR / snapshot_id / f"{sample}.vep.tsv.gz"


def clinvar_vcf(snapshot_id: str) -> Path:
    return snapshot_dir(snapshot_id) / "clinvar.vcf.gz"


def _run_vep(input_vcf: Path, clinvar: Path, out_tsv: Path, cache_version: str) -> None:
    """The one VEP invocation, shared by the full `annotate` and the incremental subset
    re-annotation. ClinVar is attached `type=exact` (coordinate+allele exact match), so the
    same variant always gets the same ClinVar columns regardless of how many are processed.

    Writes to a per-process .part file and renames on success, so the final name never holds
    a partial gzip stream. VEP used to write straight to it, and the skip check in `annotate`
    tests only existence — so a run killed mid-VEP (or the machine losing power) left a
    truncated TSV under exactly the name the next attempt would skip to and load. That became
    materially more likely once a failed scan started retrying itself rather than staying
    stopped; the retry would consume the corpse of the run that died.
    """
    custom = (
        f"file={clinvar},short_name=CLINVAR,format=vcf,type=exact,coords=0,"
        f"fields=CLNSIG%CLNREVSTAT%CLNDN%ALLELEID"
    )
    part = out_tsv.with_name(f"{out_tsv.name}.{os.getpid()}.part")
    part.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(config.VEP_BIN), "--offline", "--cache", "--dir_cache", str(config.VEP_CACHE_DIR),
        "--species", "homo_sapiens", "--assembly", config.BUILD,
        "--cache_version", cache_version,
        "--fasta", str(config.REF_FASTA),
        "-i", str(input_vcf), "--format", "vcf",
        "--tab", "--fields", VEP_FIELDS,
        "--symbol", "--sift", "b", "--polyphen", "b",
        "--af_gnomade", "--af_gnomadg",
        "--individual", "all",
        "--pick",  # one canonical consequence per variant (clean 1-row-per-variant)
        "--fork", str(config.VEP_FORKS),
        "--custom", custom,
        "--compress_output", "gzip",
        "-o", str(part), "--force_overwrite", "--no_stats",
    ]
    try:
        run(cmd, env=_vep_env())
        os.replace(part, out_tsv)          # atomic: the final name appears only when complete
    except BaseException:
        part.unlink(missing_ok=True)       # covers SystemExit/KeyboardInterrupt too
        raise


def _pinned_cache_version(snapshot_id: str) -> str:
    """The VEP cache version this snapshot was pinned to, not whatever config says today.

    Every run used config.VEP_CACHE_VERSION, so re-annotating an older snapshot after a cache
    bump produced rows that disagreed with the version recorded in that snapshot's own
    manifest and in runs.tsv — the provenance said one thing and the data was another. A scan
    of the CURRENT snapshot is unaffected (create_snapshot pins from the same config), so
    this only ever mattered where it mattered most: reproducing history.

    Refuses rather than silently substituting when the pinned cache is not installed — a
    reproduction that quietly uses a different reference is worse than one that stops.
    """
    pinned = str(load_manifest(snapshot_id).get("vep_cache_version") or "").strip()
    if not pinned:
        log.warning("Snapshot %s records no vep_cache_version — using the configured %s.",
                    snapshot_id, config.VEP_CACHE_VERSION)
        return config.VEP_CACHE_VERSION
    if pinned != config.VEP_CACHE_VERSION:
        installed = config.VEP_CACHE_DIR / "homo_sapiens" / f"{pinned}_{config.BUILD}"
        if not installed.exists():
            raise SystemExit(
                f"Snapshot {snapshot_id} is pinned to VEP cache r{pinned}, which is not "
                f"installed ({installed} missing). The configured cache is "
                f"r{config.VEP_CACHE_VERSION}; annotating with it would contradict this "
                f"snapshot's recorded provenance. Install r{pinned} or re-create the "
                f"snapshot against the current cache.")
        log.info("Snapshot %s is pinned to VEP cache r%s (configured: r%s) — using the "
                 "pinned one.", snapshot_id, pinned, config.VEP_CACHE_VERSION)
    return pinned


def _stamp_path(out: Path) -> Path:
    return out.with_name(f"{out.name}.inputs.json")


def _input_stamp(sample: str, snapshot_id: str) -> dict:
    """What this TSV's contents depend on: the normalized VCF, the snapshot's pinned ClinVar
    and VEP cache, and the field list we ask VEP for. Mirrors incremental.source_digests
    (computed here rather than imported — incremental imports this module).

    The VCF is identified by content, as it is there. This stamp used (size, whole-second
    mtime) while the scan ledger hashed, so the two could disagree about one file: a VCF
    replaced by different bytes of the same size with its timestamp preserved (`rsync -a`,
    `cp -p`, a restore) made the ledger ask for a full load while this check reused the old
    TSV — and the scan then recorded the new hash beside the old annotations."""
    meta = load_manifest(snapshot_id)
    return {
        "normalized_vcf": f"sha256:{sha256(normalized_path(sample))}",
        "clinvar": meta.get("clinvar_date", "?"),
        "vep_cache": meta.get("vep_cache_version", "?"),
        "vep_fields": VEP_FIELDS,
    }


def _legacy_fingerprint(sample: str) -> str:
    """The (size, mtime) form stamps carried before they hashed the VCF."""
    st = normalized_path(sample).stat()
    return f"{st.st_size}:{int(st.st_mtime)}"


def _read_stamp(out: Path) -> dict | None:
    try:
        return json.loads(_stamp_path(out).read_text())
    except (OSError, ValueError):
        return None


def annotate(sample: str, snapshot_id: str, force: bool = False) -> Path:
    if not config.VEP_BIN.exists():
        raise SystemExit(f"VEP not installed at {config.VEP_BIN}. Run `setup` first.")
    norm = normalized_path(sample)
    if not norm.exists():
        raise SystemExit(f"{sample} not normalized. Run `normalize {sample}` first.")

    clinvar = clinvar_vcf(snapshot_id)
    if not clinvar.exists():
        raise SystemExit(f"Snapshot {snapshot_id} has no ClinVar VCF. Create it first.")

    out = annotated_path(snapshot_id, sample)
    out.parent.mkdir(parents=True, exist_ok=True)
    want = _input_stamp(sample, snapshot_id)
    if out.exists() and not force:
        # Existence alone is not a reason to reuse: it says nothing about whether the inputs
        # that produced the file are the inputs we have now. A re-normalized sample against
        # an unchanged same-day snapshot id would otherwise keep the previous run's
        # annotations forever, and `scan --full` would quietly not re-annotate.
        have = _read_stamp(out)
        if have and not str(have.get("normalized_vcf", "")).startswith("sha256:") \
                and have.get("normalized_vcf") == _legacy_fingerprint(sample):
            # Written before the stamp hashed its input. Its fingerprint still matches, so it
            # is reused on the same terms as an unstamped TSV below, and re-stamped by content
            # so the weaker test is applied to it once rather than forever.
            have = {**have, "normalized_vcf": want["normalized_vcf"]}
            if have == want:
                _stamp_path(out).write_text(json.dumps(want, indent=2))
        if have == want:
            log.info("%s already annotated against %s (inputs unchanged) — skipping.",
                     sample, snapshot_id)
            return out
        if have is None:
            log.warning("%s/%s has no input stamp (written before stamping existed) — "
                        "reusing it. Use --force to re-annotate if it may be stale.",
                        sample, snapshot_id)
            return out
        changed = [k for k in want if have.get(k) != want[k]]
        log.info("%s/%s: re-annotating, %s changed since that TSV was written.",
                 sample, snapshot_id, ", ".join(changed))

    log.info("Running VEP: %s against snapshot %s …", sample, snapshot_id)
    cache_version = _pinned_cache_version(snapshot_id)
    _run_vep(norm, clinvar, out, cache_version)
    _stamp_path(out).write_text(json.dumps(want, indent=2))

    append_tsv(config.RUNS_TSV, {
        "snapshot_id": snapshot_id, "sample_id": sample, "run_date": today(),
        # The cache VEP was actually run with. Recording config.VEP_CACHE_VERSION here made
        # the ledger contradict both the snapshot and the run whenever an older snapshot was
        # re-annotated after a cache bump — the one case the ledger is consulted for.
        "vep_cache_version": cache_version,
        "clinvar_file": str(clinvar), "annotated_path": str(out),
    }, RUN_FIELDS)
    log.info("Wrote %s", out)
    return out


def annotate_sites(sample: str, snapshot_id: str, sites_vcf: Path) -> Path:
    """VEP-annotate a *subset* VCF (already extracted from the sample's normalized VCF) against
    a snapshot's ClinVar. Returns a TSV path that `load.load_subset` merges in. The VEP path is
    identical to the full run, so the rows are byte-identical to what a full annotate produces."""
    clinvar = clinvar_vcf(snapshot_id)
    out = sites_vcf.with_suffix(".vep.tsv.gz")
    log.info("Running VEP on %d-byte subset for %s against %s …",
             sites_vcf.stat().st_size, sample, snapshot_id)
    _run_vep(sites_vcf, clinvar, out, _pinned_cache_version(snapshot_id))
    return out
