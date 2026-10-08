"""Stage 2: ingest a raw per-sample VCF as immutable source of truth."""

from __future__ import annotations

import gzip
import os
import re
import shutil
from pathlib import Path

from . import config
from .util import append_tsv, find_row, log, read_tsv, run, sha256, today

SAMPLE_FIELDS = [
    "sample_id", "relationship", "sex", "build", "source",
    "sha256", "raw_path", "ingest_date",
]

# Off-machine archive mount (optional: a NAS share, an external volume). Set RAW_ARCHIVE to
# have every ingested raw file mirrored there.
RAW_ARCHIVE = os.environ.get("RAW_ARCHIVE")


def _header_lines(vcf: Path) -> list[str]:
    opener = gzip.open if vcf.suffix == ".gz" else open
    out = []
    with opener(vcf, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            out.append(line.rstrip("\n"))
    return out


# chr1's length is the assembly's fingerprint — it differs in every build we might be handed.
_CHR1_LENGTH = {"GRCh37": 249250621, "GRCh38": 248956422}
_ASSEMBLY_BY_CHR1 = {249250621: "GRCh37/hg19", 248956422: "GRCh38/hg38",
                     248387328: "CHM13v2.0 (T2T)"}


def _chr1_length(header: list[str]) -> int | None:
    """The declared length of chromosome 1, from the ##contig header lines."""
    for line in header:
        if not line.startswith("##contig="):
            continue
        # ID=1 or ID=chr1, and nothing else: the trailing [,>] keeps ID=10, ID=11 and
        # ID=chr1_KI270706v1_random out.
        if not re.search(r"ID=(?:chr)?1[,>]", line):
            continue
        m = re.search(r"length=(\d+)", line)
        if m:
            return int(m.group(1))
    return None


def _looks_like_build(header: list[str]) -> tuple[bool, str]:
    """(ok, why-not) — guard against ingesting a mismatched reference build.

    A declared chr1 length is authoritative and VETOES: it is a fact about the assembly the
    caller used, where a build name in the header is only a claim, and the two disagree more
    often than one would like. The old test was a flat OR over the whole header text, so any
    mention of "GRCh38" anywhere — a ##bcftools_viewCommand line quoting a path, a ##source
    naming a pipeline — outvoted a contig length that said GRCh37 outright. Only when the
    header declares no usable chr1 contig do we fall back to the name.
    """
    length = _chr1_length(header)
    expected = _CHR1_LENGTH.get(config.BUILD)
    if length is not None and expected is not None:
        if length == expected:
            return True, ""
        return False, (
            f"its chromosome 1 is declared as {length:,} bp — "
            f"{_ASSEMBLY_BY_CHR1.get(length, 'an assembly we do not recognise')}, not "
            f"{config.BUILD} ({expected:,} bp)")
    names = ("GRCh37", "hg19") if config.BUILD == "GRCh37" else ("GRCh38", "hg38")
    text = "\n".join(header)
    if any(n in text for n in names):
        return True, ""
    return False, ("it declares no chromosome-1 contig length and no build name in its header")


def _sample_names(header: list[str]) -> list[str]:
    """Sample columns from the #CHROM line (fields 10 onwards)."""
    for line in header:
        if line.startswith("#CHROM"):
            return line.split("\t")[9:]
    return []


def _copy_readonly(src: Path, dest: Path, what: str) -> None:
    """Copy src → dest and make dest read-only, restartably.

    A read-only destination left behind by an interrupted ingest is not an error to abort on:
    if it already matches the source byte-for-byte the copy is simply already done. Otherwise
    the copy goes to a temp file in the same directory and is renamed over the destination,
    so an interrupted run never leaves a partial file under the final name.
    """
    if dest.exists():
        if dest.stat().st_size == src.stat().st_size and sha256(dest) == sha256(src):
            log.info("%s already present and identical (%s) — reusing.", what, dest)
            return
        raise SystemExit(
            f"{dest} already exists and differs from {src}. Refusing to overwrite an "
            f"ingested {what}; move it aside deliberately if the replacement is intended.")
    tmp = dest.with_name(dest.name + ".part")
    tmp.unlink(missing_ok=True)
    log.info("Copying %s → %s", src, dest)
    try:
        shutil.copy2(src, tmp)
        os.chmod(tmp, 0o444)      # read-only: immutable source of truth
        os.replace(tmp, dest)     # atomic: dest appears complete or not at all
    except BaseException:
        try:
            os.chmod(tmp, 0o644)
        except OSError:
            pass
        tmp.unlink(missing_ok=True)
        raise


def ingest(vcf_path: str, sample: str, relation: str, sex: str,
           warn_unmirrored: bool = True) -> None:
    src = Path(vcf_path).expanduser().resolve()
    if not src.exists():
        raise SystemExit(f"Input VCF not found: {src}")

    if find_row(config.SAMPLES_TSV, sample_id=sample):
        log.info("Sample %s already ingested — skipping (idempotent).", sample)
        return

    header = _header_lines(src)
    ok, why = _looks_like_build(header)
    if not ok:
        raise SystemExit(
            f"{src.name} does not look like {config.BUILD}: {why}. Refusing to ingest a "
            f"mismatched build — every downstream coordinate would be wrong. Lift it over "
            f"first, or ingest it under a suffixed id once it is on {config.BUILD}.")

    # One VCF, one sample. `annotate` runs VEP with --individual all, which emits a row per
    # sample, and `load` files every row under the single sample_id recorded here — so a
    # multisample VCF silently attributes everyone's genotypes to one person.
    names = _sample_names(header)
    if not names:
        raise SystemExit(
            f"{src.name} has no sample columns (a sites-only VCF). There are no genotypes "
            f"to ingest — zygosity, panel status and every report downstream need them.")
    if len(names) > 1:
        shown = ", ".join(names[:5]) + (", …" if len(names) > 5 else "")
        raise SystemExit(
            f"{src.name} carries {len(names)} sample columns ({shown}) but ingest records "
            f"one sample id ({sample!r}). VEP's --individual all would emit a row per "
            f"sample and all of them would be filed under {sample!r}. Split it first: "
            f"bcftools view -s <name> -Oz -o <name>.vcf.gz {src.name}")
    if names[0] != sample:
        log.info("%s: VCF sample column is %r, ingesting as %r.", src.name, names[0], sample)

    dest_dir = config.RAW_DIR / sample
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / src.name
    # The idempotency check above is the manifest row, which is written LAST — after a
    # multi-GB copy and a sha256 over it. A crash in that window used to leave a 0444
    # destination and no manifest row, so the retry took the copy path again and died with
    # PermissionError on its own read-only file, with no way forward but a manual chmod.
    # Copy to a temp name and rename into place: an interrupted attempt leaves only the temp
    # file, and a complete-but-unrecorded one is verified and reused.
    _copy_readonly(src, dest, "raw file")
    digest = sha256(dest)

    if RAW_ARCHIVE:
        arch = Path(RAW_ARCHIVE) / "raw" / sample / src.name
        arch.parent.mkdir(parents=True, exist_ok=True)
        _copy_readonly(dest, arch, "archive mirror")
        log.info("Mirrored to archive: %s", arch)
    elif warn_unmirrored:
        log.warning("RAW_ARCHIVE not set — raw file is NOT mirrored off-machine.")

    append_tsv(config.SAMPLES_TSV, {
        "sample_id": sample, "relationship": relation, "sex": sex,
        "build": config.BUILD, "source": str(src), "sha256": digest,
        "raw_path": str(dest), "ingest_date": today(),
    }, SAMPLE_FIELDS)

    _upsert_family(sample, sex)
    log.info("Ingested %s (sha256=%s…). Next: `normalize %s`.", sample, digest[:12], sample)


_FAMILY_HEADER = ["id", "display_name", "sex", "father", "mother", "partner"]


def _upsert_family(sample: str, sex: str) -> None:
    """Add the PERSON to the family graph (family.tsv), then regenerate pedigree.ped from it.

    Relationships (parents, partner) are set by editing family.tsv afterwards — the ingest
    only seeds the person with their sex. This replaces the old stub-PED writer: pedigree.ped
    is now a *generated* artifact (pipeline/pedigree.py).

    An id that resolves to someone already in the graph is an alternate CALLSET of theirs
    (`Adam_t2t` -> `Adam`) and adds no row: callsets live in samples.tsv only. This is
    how `Jan_t2t` and `Maria_t2t` became parentless people in the first place, giving the
    family a duplicate tree — see pipeline/pedigree.py."""
    from . import pedigree
    fp = config.FAMILY_FILE
    existing = read_tsv(fp) if fp.exists() else []
    person = pedigree.load(fp).person_of(sample) if fp.exists() else sample
    if person != sample:
        log.info("%s is another callset of %s — family graph unchanged.", sample, person)
    elif not any(r.get("id") == sample for r in existing):
        fp.parent.mkdir(parents=True, exist_ok=True)
        if not fp.exists():
            fp.write_text("\t".join(_FAMILY_HEADER) + "\n")
        with open(fp, "a") as fh:
            fh.write("\t".join([sample, sample, sex.lower(), "", "", ""]) + "\n")
        log.info("Added %s to the family graph — edit %s to set parents/partner.", sample, fp)
    pedigree.generate_pedigree_ped()
