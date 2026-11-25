"""Stage: full classical HLA typing via T1K.

Supersedes the tag-SNP proxies in `pipeline/hla.py` with classical alleles
(HLA-A/B/C/DPA1/DPB1/DQA1/DQB1/DRB1) typed from the CRAM. samtools pulls the read pairs that
touch the extended MHC, plus every pair with an unmapped mate, and T1K genotypes them against
genomic IPD-IMGT/HLA allele sequences. T1K ships via bioconda, so (like Manta) it runs in a
pinned container, here with networking off. The allele index is built locally from a pinned
IPD-IMGT/HLA release, because that licence needs permission to redistribute modified data.

Re-derives the same two risk calls `hla.py` makes from tag SNPs (HLA-B27, celiac
DQ2.5/DQ8) directly from the typed alleles, so the two can be cross-checked.
"""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

from . import reportpaths
from . import config, hla
from .call import aligned_cram
from .util import log, require_tools, run, sha256, write_report


def _docker(mounts: list[str], *cmd: str) -> list[str]:
    """A T1K container command: no network, and files owned by the calling user."""
    args = ["docker", "run", "--rm", "--network", "none",
            "--user", f"{os.getuid()}:{os.getgid()}"]
    for mount in mounts:
        args += ["-v", mount]
    return [*args, config.T1K_IMAGE, *cmd]


def ensure_t1k_index() -> Path:
    """The T1K allele index, built once from the pinned, checksummed IPD-IMGT/HLA release."""
    index = config.T1K_INDEX_DIR / "hlaidx_dna_seq.fa"
    if index.exists():
        return index
    require_tools("docker")
    config.T1K_DIR.mkdir(parents=True, exist_ok=True)
    source = config.T1K_DIR / "hla.dat.zip"
    if not source.exists():
        from .setup import _download
        _download(config.T1K_IMGT_URL, source, f"IPD-IMGT/HLA {config.T1K_IMGT_VERSION}")
    if sha256(source) != config.T1K_IMGT_SHA256:
        raise SystemExit(f"{source} does not match the pinned IPD-IMGT/HLA "
                         f"{config.T1K_IMGT_VERSION} checksum. Delete it and retry.")
    with tempfile.TemporaryDirectory(prefix=".build-", dir=config.T1K_DIR) as tmp:
        work = Path(tmp)
        with zipfile.ZipFile(source) as archive, archive.open("hla.dat") as src, \
                open(work / "hla.dat", "wb") as dst:
            shutil.copyfileobj(src, dst)
        (work / "index").mkdir()
        run(_docker([f"{work}:/build"], "t1k-build.pl", "-d", "/build/hla.dat",
                    "-o", "/build/index", "--prefix", "hlaidx"),
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if not (work / "index" / index.name).is_file():
            raise SystemExit("T1K did not build its HLA allele index.")
        os.replace(work / "index", config.T1K_INDEX_DIR)
    return index


def _mhc_contig(cram: Path) -> str:
    """The CRAM's chromosome 6 name. The MHC coordinates below are GRCh38, so another build
    must stop here rather than extract the wrong region and type it without complaint."""
    out = subprocess.run(["samtools", "view", "-H", str(cram)],
                         capture_output=True, check=True, text=True)
    for line in out.stdout.splitlines():
        if not line.startswith("@SQ"):
            continue
        tags = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
        if tags.get("SN") in ("6", "chr6"):
            if tags.get("LN") != str(config.HLA_CHR6_LENGTH):
                raise SystemExit(f"{cram}: chromosome 6 is not GRCh38 length "
                                 f"({tags.get('LN')} vs {config.HLA_CHR6_LENGTH}).")
            return tags["SN"]
    raise SystemExit(f"No chromosome 6 (6 or chr6) in the {cram} header.")


def _read_filter(contig: str) -> str:
    """Pairs touching the extended MHC, and pairs with an unmapped mate. HLA reads too
    divergent for the primary assembly usually end up unmapped."""
    start, end = config.HLA_MHC_START, config.HLA_MHC_END
    return (f'flag.unmap || flag.munmap '
            f'|| (rname=="{contig}" && endpos>={start} && pos<={end}) '
            f'|| (rnext=="{contig}" && pnext>={start} && pnext<={end})')


def _extract_reads(cram: Path, contig: str, work: Path, threads: int, log_file) -> None:
    helpers = min(4, threads)
    steps = [
        ["samtools", "view", "-@", threads, "-T", config.REF_FASTA, "-u", "-F", "2304",
         "-e", _read_filter(contig), cram],
        ["samtools", "collate", "-@", helpers, "-u", "-O", "-T", work / "collate", "-"],
        ["samtools", "fastq", "-@", helpers, "-n", "-1", work / "r1.fq.gz",
         "-2", work / "r2.fq.gz", "-0", "/dev/null", "-s", "/dev/null"],
    ]
    script = "set -o pipefail; " + " | ".join(shlex.join(map(str, s)) for s in steps)
    run(["bash", "-c", script], stdout=log_file, stderr=subprocess.STDOUT)


def _run_t1k(cram: Path, out_dir: Path, sample: str) -> Path:
    """Extract reads and type them; keep T1K's call files and none of the temporary reads.
    `HLA_TYPE_TMP` can move those reads to a larger disk or a RAM disk."""
    require_tools("docker", "samtools")
    index = ensure_t1k_index()
    contig = _mhc_contig(cram)
    threads = os.cpu_count() or 4
    scratch = Path(config.HLA_TYPE_TMP) if config.HLA_TYPE_TMP else out_dir
    with tempfile.TemporaryDirectory(prefix=f".t1k-{sample}-", dir=scratch) as tmp, \
            open(config.LOGS_DIR / f"hla_typing_{sample}.log", "wb") as log_file:
        work = Path(tmp)
        (work / "out").mkdir()
        _extract_reads(cram, contig, work, threads, log_file)
        run(_docker([f"{work}:/work", f"{index.parent}:/index:ro"], "run-t1k",
                    "-1", "/work/r1.fq.gz", "-2", "/work/r2.fq.gz", "--preset", "hla-wgs",
                    "-f", f"/index/{index.name}", "-t", str(threads),
                    "-o", sample, "--od", "/work/out"),
            stdout=log_file, stderr=subprocess.STDOUT)
        # The genotype file goes last: its presence marks a finished run.
        for suffix in ("allele.tsv", "allele.vcf", "genotype.tsv"):
            dest = out_dir / f"{sample}_{suffix}"
            part = dest.with_name(dest.name + ".part")
            shutil.copyfile(work / "out" / dest.name, part)
            os.replace(part, dest)
    return out_dir / f"{sample}_genotype.tsv"


def _four_digit(allele: str) -> str:
    """'A*01:01:01' -> 'A*01:01' (fields beyond the 4-digit resolution are synonymous
    substitutions / noncoding differences T1K still reports; risk calls key off the first
    two fields)."""
    gene, _, fields = allele.partition("*")
    parts = fields.split(":")
    return f"{gene}*{':'.join(parts[:2])}"


def _parse_t1k(path: Path) -> dict[str, dict]:
    """{locus: {"alleles": [{"names", "quality"}], "low": bool}} from a T1K genotype file.
    One called allele can list several tied alleles. It is low confidence when it scores
    under T1K_MIN_QUALITY or its ties span more than one 2-field allele. T1K's docs say to
    ignore quality <= 0, and losing an allele that way leaves the locus low confidence."""
    calls: dict[str, dict] = {}
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        gene = fields[0].removeprefix("HLA-")
        if gene not in config.HLA_TYPE_GENES or len(fields) < 8:
            continue
        called, rejected = [], False
        for allele, quality in ((fields[2], fields[4]), (fields[5], fields[7])):
            if allele == ".":
                continue
            if float(quality) <= 0:
                rejected = True
                continue
            names = list(dict.fromkeys(_four_digit(a.removeprefix("HLA-"))
                                       for a in allele.split(",")))
            called.append({"names": names, "quality": float(quality)})
        if called:
            low = rejected or any(a["quality"] < config.T1K_MIN_QUALITY or len(a["names"]) > 1
                                  for a in called)
            calls[gene] = {"alleles": called, "low": low}
    return calls


# HLA-B*27 subtypes NOT associated with ankylosing spondylitis (Khan 2013; the classic
# exceptions to the B27 association).
_B27_NON_ASSOCIATED = frozenset({"B*27:06", "B*27:09"})
# DQA1 alleles that form the DQ2.5 heterodimer with a DQB1*02 beta chain — *05:05 encodes the
# same mature alpha chain as *05:01.
_DQA1_DQ25 = frozenset({"DQA1*05:01", "DQA1*05:05"})
# *02:02 differs from *02:01 by one residue outside the peptide-binding groove. It is the beta
# chain of DR7-DQ2.2, so it is the one a DR5-DQ7 / DR7-DQ2.2 person pairs with DQA1*05:05.
_DQB1_DQ25 = frozenset({"DQB1*02:01", "DQB1*02:02"})


def _risk_calls(genotype: dict[str, list[str]],
                tags: dict[str, bool | None] | None = None,
                low: set[str] | frozenset[str] = frozenset()) -> list[str]:
    """Re-derive the same two tag-SNP-based calls in hla.py from typed alleles. `tags` is
    `hla.tag_status`: where it disagrees with the typing, the call is unresolved, because a
    thinly covered locus can type to the wrong allele and that must not read as a negative.
    `low` names low-confidence loci, whose calls stand only where the tag SNP agrees."""
    tags = tags or {}
    b = [_four_digit(a) for a in genotype.get("B", [])]
    dqa1 = [_four_digit(a) for a in genotype.get("DQA1", [])]
    dqb1 = [_four_digit(a) for a in genotype.get("DQB1", [])]

    calls = []
    # B*27:06 and B*27:09 are the well-documented non-associated subtypes: they differ from
    # B*27:05 in the B pocket / alpha-2 domain and are NOT linked to ankylosing spondylitis
    # (B*27:06 in Southeast Asia, B*27:09 in Sardinia). A bare "B*27" prefix test reported
    # them as positive for a risk they don't carry.
    b27_alleles = [a for a in b if a.startswith("B*27")]
    assoc = [a for a in b27_alleles if a not in _B27_NON_ASSOCIATED]
    b27 = bool(assoc)
    if not b:
        # No call means too few reads, which says nothing about whether B*27 is there.
        calls.append("**HLA-B27**: not determined (no B calls: HLA-B could not be typed)")
    elif tags.get("B27") is True and not b27_alleles:
        calls.append(f"**HLA-B27**: unresolved (tag SNP positive, typed alleles show no B*27: "
                     f"{', '.join(b)})")
    elif tags.get("B27") is False and b27_alleles:
        calls.append(f"**HLA-B27**: unresolved (typed {', '.join(b27_alleles)}, "
                     f"tag SNP negative)")
    elif "B" in low:
        if tags.get("B27") is False:
            verdict = "absent (low-confidence typing, tag SNP agrees"
        elif tags.get("B27") is True and assoc:
            verdict = "present (low-confidence typing, tag SNP agrees"
        else:
            verdict = "not determined (low-confidence HLA-B typing"
        calls.append(f"**HLA-B27**: {verdict}: {', '.join(b)})")
    elif b27_alleles and not assoc:
        calls.append(f"**HLA-B27**: present but **non-associated subtype** "
                     f"({', '.join(b27_alleles)}) — not linked to ankylosing spondylitis")
    else:
        calls.append(f"**HLA-B27**: {'present' if b27 else 'absent'} ({', '.join(b)})")

    # The heterodimer forms whether its two chains come from one haplotype (DR3: DQA1*05:01 +
    # DQB1*02:01) or from both (DR5-DQ7 / DR7-DQ2.2: DQA1*05:05 + DQB1*02:02), so this needs
    # no phasing. Requiring *02:01 reported every trans carrier as "neither".
    dq25 = bool(set(dqa1) & _DQA1_DQ25) and bool(set(dqb1) & _DQB1_DQ25)
    dq8 = "DQB1*03:02" in dqb1
    # "neither" needs DQB1, and DQA1 too wherever a DQ2 beta chain could complete DQ2.5.
    excluded = bool(dqb1) and (bool(dqa1) or not set(dqb1) & _DQB1_DQ25)
    typed = [h for h, p in [("DQ2.5", dq25), ("DQ8", dq8)] if p]
    # A typed DQ2.5 with a negative tag stands: the tag cannot see DQ2.5 formed in trans.
    missed = [h for h in ("DQ2.5", "DQ8") if tags.get(h) is True and h not in typed]
    unbacked = dq8 and tags.get("DQ8") is False
    if missed or unbacked:
        parts = [f"typed {h}" for h in typed if not (h == "DQ8" and unbacked)]
        parts += [f"tag SNP shows {h}, typing did not confirm it" for h in missed]
        parts += ["typing shows DQ8, tag SNP does not"] if unbacked else []
        hap = f"unresolved ({'; '.join(parts)})"
    elif set(low) & {"DQA1", "DQB1"}:
        # A negative tag cannot back "neither", because it does not see DQ2.5 formed in trans.
        if typed and all(tags.get(h) is True for h in typed):
            hap = f"{' + '.join(typed)} (low-confidence typing, tag SNP agrees)"
        else:
            hap = "not determined (low-confidence DQ typing)"
    elif typed:
        hap = " + ".join(typed)
    elif any(a.startswith("DQA1*05") for a in dqa1) and any(a.startswith("DQB1*02") for a in dqb1):
        # One wrong second field away from DQ2.5, so this cannot count as an exclusion.
        hap = "not determined (DQA1*05 with a DQB1*02 allele outside the DQ2.5 set)"
    else:
        hap = "neither" if excluded else "not determined"
    calls.append(f"**Celiac haplotype** (DQ2.5/DQ8): {hap} "
                 f"(DQA1 {', '.join(dqa1) or '?'}; DQB1 {', '.join(dqb1) or '?'})")
    return calls


def _locus_row(gene: str, call: dict | None) -> str:
    if not call:
        return f"| HLA-{gene} | no call | – | – |"
    alleles = ", ".join("/".join(a["names"]) for a in call["alleles"])
    quality = ", ".join(f"{a['quality']:.0f}" for a in call["alleles"])
    return f"| HLA-{gene} | {alleles} | {quality} | {'low' if call['low'] else 'ok'} |"


def hla_typing_report(sample: str, force: bool = False) -> Path:
    cram = aligned_cram(sample)
    if not cram.exists():
        raise SystemExit(f"{cram} missing. Run `align {sample}` first.")

    out_dir = config.HLA_TYPE_DIR / sample
    summary = reportpaths.genome_report("hla_typing", sample)
    if summary.exists() and not force:
        log.info("%s HLA typing already reported (%s) — skipping (use --force).",
                 sample, summary)
        return summary
    out_dir.mkdir(parents=True, exist_ok=True)

    genotype_tsv = out_dir / f"{sample}_genotype.tsv"
    if force or not genotype_tsv.exists():
        _run_t1k(cram, out_dir, sample)
    if not genotype_tsv.exists():
        raise SystemExit(f"T1K did not produce {genotype_tsv} — check "
                         f"{config.LOGS_DIR}/hla_typing_{sample}.log")

    calls = _parse_t1k(genotype_tsv)
    genotype = {g: [n for a in c["alleles"] for n in a["names"]] for g, c in calls.items()}
    low = {g for g, c in calls.items() if c["low"]}
    for gene in config.HLA_TYPE_GENES:
        log.info("  HLA-%-6s %s", gene, ", ".join(genotype.get(gene, ["no call"])))
    tags = hla.tag_status(sample)

    lines = [
        f"# HLA typing (classical, T1K) — {sample}",
        "",
        "| Locus | Alleles (2-field) | Quality | Confidence |",
        "|-------|-------------------|---------|------------|",
        *[_locus_row(g, calls.get(g)) for g in config.HLA_TYPE_GENES],
        "",
        "## Risk calls re-derived from typed alleles",
        *[f"- {c}" for c in _risk_calls(genotype, tags, low)],
        "",
        "## Notes",
        f"- Typed with T1K {config.T1K_VERSION} against IPD-IMGT/HLA "
        f"{config.T1K_IMGT_VERSION}, from read pairs in the extended MHC and pairs with an "
        "unmapped mate. One allele at a locus usually means a homozygous call.",
        f"- `low` confidence means an allele scored under quality {config.T1K_MIN_QUALITY}, "
        "an allele was rejected, or T1K could not separate several 2-field alleles (shown "
        "with `/`). A risk call from a `low` locus stands only where the tag SNP agrees.",
        "- `no call` means T1K could not type that locus. It is not evidence that any "
        "allele is absent.",
        ("- Risk calls were cross-checked against the tag SNPs in the `hla` report. "
         "`unresolved` means the two disagree, and neither side counts as a negative. A typed "
         "DQ2.5 with a negative tag can be real, because the tag cannot see DQ2.5 formed in "
         "trans." if tags else
         "- The tag SNPs could not be read for this run, so these risk calls were not "
         "cross-checked against the `hla` report."),
        "- Typing from short reads can still be wrong, especially at `low` loci. Check "
        "alleles against relatives or an orthogonal assay before relying on them.",
        "- Research-grade; clinically actionable HLA results (e.g. B27, celiac) warrant "
        "confirmatory testing.",
        "",
    ]
    write_report(summary, lines)
    log.info("HLA typing report: %s", summary)
    return summary
