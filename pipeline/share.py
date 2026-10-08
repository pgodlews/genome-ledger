"""`share`: a copy of a person's VCF or CRAM that is fit to hand to someone else.

The working files keep a record of how they were made, and that record is personal. A CRAM's
`@PG` lines quote every command that produced it with full paths, the home folder included;
its `@SQ UR:` fields point at the reference on this machine; its read group names the callset
id, which can be an old or internal one. A VCF carries `##GATKCommandLine`,
`##bcftools_*Command` and `##reference` lines with the same paths and run dates, and the sample
column of an ingested callset is often the provider's kit id. None of that belongs in a file
sent to a doctor, a researcher or an upload service.

The working files are left alone — that provenance is how the pipeline explains its own
results. This writes a separate copy under SHARE_DIR whose header is rebuilt from an allowlist:

  VCF   ##fileformat, ##INFO, ##FORMAT, ##FILTER, ##ALT and ##contig lines are kept, a
        `##reference=<build>` line replaces the path, and the sample column becomes the
        person's id (or --name).
  CRAM  @HD is kept; @SQ keeps SN, LN, M5 and AS; @RG keeps ID and PL with SM and LB set to
        the name; every @PG and @CO line is dropped.

What it cannot remove, and what the person sharing should know:
  * The variants and reads ARE the genome. Nothing here makes a file anonymous.
  * CRAM read names carry the sequencer's instrument, run and flowcell ids.
  * The CRAM read-group ID is left as it is, because every read refers to it; it is the
    sample id this pipeline used.
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path

from . import config
from .call import aligned_cram, sample_vcf
from .util import log, require_tools, run

_VCF_KEEP = re.compile(r"^##(fileformat|INFO|FORMAT|FILTER|ALT|contig)=")
_SQ_KEEP = {"SN", "LN", "M5", "AS"}
_RG_KEEP = {"ID", "PL"}


def _person(sample: str) -> str:
    try:
        from . import pedigree
        return pedigree.load().person_of(sample)
    except Exception:  # noqa: BLE001 — no family graph: the id is the best name there is
        return sample


def vcf_header(lines: list[str], name: str) -> list[str]:
    """The shareable header: allowlisted meta lines, the build, one renamed sample column."""
    meta = [ln for ln in lines if _VCF_KEEP.match(ln)]
    fileformat = [ln for ln in meta if ln.startswith("##fileformat=")]
    rest = [ln for ln in meta if not ln.startswith("##fileformat=")]
    chrom = next((ln for ln in lines if ln.startswith("#CHROM")), None)
    if chrom is None:
        raise SystemExit("VCF header has no #CHROM line.")
    cols = chrom.split("\t")
    if len(cols) != 10:
        raise SystemExit(f"Expected one sample column, found {max(len(cols) - 9, 0)}; "
                         "`share` writes single-person files.")
    return [*fileformat, f"##reference={config.BUILD}", *rest, "\t".join(cols[:9] + [name])]


def cram_header(lines: list[str], name: str) -> list[str]:
    """The shareable SAM header: allowlisted @HD/@SQ/@RG fields, no @PG or @CO."""
    out = []
    for ln in lines:
        tag, *fields = ln.split("\t")
        if tag == "@HD":
            out.append(ln)
        elif tag == "@SQ":
            out.append("\t".join([tag, *[f for f in fields if f[:2] in _SQ_KEEP]]))
        elif tag == "@RG":
            kept = [f for f in fields if f[:2] in _RG_KEEP]
            out.append("\t".join([tag, *kept, f"SM:{name}", f"LB:{name}"]))
        # @PG, @CO and anything else: dropped
    return out


def share(sample: str, kind: str = "vcf", name: str | None = None) -> Path:
    name = name or _person(sample)
    if not re.fullmatch(r"[A-Za-z0-9._-]+", name):
        raise SystemExit(f"Name '{name}' must be letters, digits, '.', '_' or '-'.")
    out_dir = config.SHARE_DIR / name
    out_dir.mkdir(parents=True, exist_ok=True)

    if kind == "vcf":
        require_tools("bcftools", "tabix")
        src = sample_vcf(sample)
        if not src.exists():
            raise SystemExit(f"No VCF for {sample} ({src}). Run `call` or `normalize` first.")
        header = run(["bcftools", "view", "-h", str(src)], capture_output=True,
                     text=True).stdout.splitlines()
        new = vcf_header(header, name)
        out = out_dir / f"{name}.{config.BUILD}.vcf.gz"
        with tempfile.NamedTemporaryFile("w", suffix=".hdr", delete=False) as fh:
            fh.write("\n".join(new) + "\n")
        try:
            run(["bcftools", "reheader", "-h", fh.name, "-o", str(out), str(src)])
        finally:
            Path(fh.name).unlink(missing_ok=True)
        run(["tabix", "-f", "-p", "vcf", str(out)])
    elif kind == "cram":
        require_tools("samtools")
        src = aligned_cram(sample)
        if not src.exists():
            raise SystemExit(f"No CRAM for {sample} ({src}) — only samples aligned by this "
                             "pipeline have one.")
        header = run(["samtools", "view", "-H", "--no-PG", str(src)], capture_output=True,
                     text=True).stdout.splitlines()
        new = cram_header(header, name)
        out = out_dir / f"{name}.cram"
        with tempfile.NamedTemporaryFile("w", suffix=".sam", delete=False) as fh:
            fh.write("\n".join(new) + "\n")
        try:
            with open(out, "wb") as dest:
                run(["samtools", "reheader", "--no-PG", fh.name, str(src)], stdout=dest)
        except BaseException:
            out.unlink(missing_ok=True)
            raise
        finally:
            Path(fh.name).unlink(missing_ok=True)
        run(["samtools", "index", str(out)])
    else:
        raise SystemExit(f"Unknown kind '{kind}' (vcf or cram).")

    log.info("Wrote %s — header rebuilt from an allowlist (no command lines, paths or "
             "internal ids). The contents are still %s's genome; read names in a CRAM still "
             "carry the sequencer's run and flowcell ids.", out, name)
    return out
