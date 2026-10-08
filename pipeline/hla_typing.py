"""Stage: full classical HLA typing via arcasHLA.

Supersedes the tag-SNP proxies in `pipeline/hla.py` with real 4-digit classical alleles
(HLA-A/B/C/DPA1/DPB1/DQA1/DQB1/DRB1) called straight from the CRAM: arcasHLA extracts
chr6 (+ HLA alt/decoy, if present) reads with samtools, then kallisto-pseudoaligns them
against the IMGT/HLA reference to resolve genotypes. arcasHLA ships via bioconda, so
(like Manta) it runs containerized rather than adding conda to this project's toolchain —
wherever the CRAM + Docker live.

Re-derives the same two risk calls `hla.py` makes from tag SNPs (HLA-B27, celiac
DQ2.5/DQ8) directly from the typed alleles, so the two can be cross-checked.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

from . import reportpaths
from . import config
from .call import aligned_cram
from .util import log, require_tools, run, write_report


def _cram_ref_dir(cram: Path) -> str | None:
    """Return the directory of the CRAM header's embedded @SQ UR tag, or None if it
    already matches our own REF_FASTA's directory. The CRAM stores an absolute path
    from the aligning machine baked in at align time — on a different machine that
    path may not exist. Unlike Manta, which takes an explicit --referenceFasta, arcasHLA's
    raw samtools calls expose no such flag, and htslib's REF_CACHE/REF_PATH env vars did
    not take priority over the UR in testing — so instead we bind-mount our own reference
    directory at the CRAM's embedded path, the same trick sv.py uses for GENOMES_ROOT itself."""
    out = subprocess.run(
        ["docker", "run", "--rm", "-v", f"{cram.parent}:{cram.parent}",
         config.ARCASHLA_IMAGE, "samtools", "view", "-H", str(cram)],
        capture_output=True, check=True, text=True)
    for line in out.stdout.splitlines():
        if line.startswith("@SQ") and "\tUR:" in line:
            ur = next(f[3:] for f in line.split("\t") if f.startswith("UR:"))
            ref_dir = str(Path(ur).parent)
            return None if ref_dir == str(config.REF_FASTA.parent) else ref_dir
    raise SystemExit(f"No @SQ UR tag found in {cram} header")


def _run_arcashla(cram: Path, out_dir: Path, sample: str) -> Path:
    """Run arcasHLA extract + genotype in Docker; return the genotype JSON path.
    Mounts GENOMES_ROOT at the same path inside the container (matches sv.py's Manta
    pattern) so the CRAM's embedded reference URI still resolves and no CRAM->BAM
    conversion is needed."""
    require_tools("docker")
    root = config.GENOMES_ROOT
    uid, gid = os.getuid(), os.getgid()
    threads = str(os.cpu_count() or 4)
    ref_dir = _cram_ref_dir(cram)

    # arcasHLA's index_bam() only ever checks for a `.bai` sidecar, even on CRAM input
    # (where samtools produces `.crai`) — so indexing silently succeeds but the script
    # still exits "unable to index bam file." Nothing downstream reads this sentinel's
    # contents (only its existence gates the check), so a symlink to the real .crai
    # satisfies it without ever converting the CRAM to BAM.
    bai_sentinel = Path(f"{cram}.bai")
    # lexists, not exists(): a DANGLING symlink (target .crai deleted/regenerated) reports
    # exists()==False, and symlink_to would then raise FileExistsError on every re-run.
    if not os.path.lexists(bai_sentinel):
        bai_sentinel.symlink_to(f"{cram.name}.crai")
    genes = ",".join(config.HLA_TYPE_GENES)
    inner = (
        f"arcasHLA extract {cram} -o {out_dir} -t {threads} --unmapped -v && "
        f"arcasHLA genotype {out_dir}/{sample}.extracted.1.fq.gz "
        f"{out_dir}/{sample}.extracted.2.fq.gz -g {genes} -o {out_dir} -t {threads} -v"
    )
    cmd = ["docker", "run", "--rm",
           "-v", f"{root}:{root}",
           "-w", str(out_dir),
           "-e", f"HOME={out_dir}"]
    if ref_dir is not None:
        cmd += ["-v", f"{config.REF_FASTA.parent}:{ref_dir}:ro"]
    cmd += ["--user", f"{uid}:{gid}",
            config.ARCASHLA_IMAGE, "bash", "-c", inner]
    run(cmd,
        stdout=open(config.LOGS_DIR / f"hla_typing_{sample}.log", "wb"),
        stderr=subprocess.STDOUT)
    return out_dir / f"{sample}.genotype.json"


def _four_digit(allele: str) -> str:
    """'A*01:01:01' -> 'A*01:01' (fields beyond the 4-digit resolution are synonymous
    substitutions / noncoding differences arcasHLA still reports; risk calls key off
    the first two fields)."""
    gene, _, fields = allele.partition("*")
    parts = fields.split(":")
    return f"{gene}*{':'.join(parts[:2])}"


# HLA-B*27 subtypes NOT associated with ankylosing spondylitis (Khan 2013; the classic
# exceptions to the B27 association).
_B27_NON_ASSOCIATED = frozenset({"B*27:06", "B*27:09"})
# DQA1 alleles that form the DQ2.5 heterodimer with DQB1*02:01 — *05:05 encodes the same
# mature alpha chain as *05:01.
_DQA1_DQ25 = frozenset({"DQA1*05:01", "DQA1*05:05"})


def _risk_calls(genotype: dict[str, list[str]]) -> list[str]:
    """Re-derive the same two tag-SNP-based calls in hla.py from typed alleles."""
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
    if b27_alleles and not assoc:
        calls.append(f"**HLA-B27**: present but **non-associated subtype** "
                     f"({', '.join(b27_alleles)}) — not linked to ankylosing spondylitis")
    else:
        calls.append(f"**HLA-B27**: {'present' if b27 else 'absent'} "
                     f"({', '.join(b) or 'no B calls'})")

    # DQA1*05:05 encodes the same mature DQ alpha chain as *05:01 and pairs with DQB1*02:01
    # to form the identical DQ2.5 heterodimer (it is the allele carried on the common
    # DR7-DQ2.2 / DR5-DQ7 trans haplotype). Testing only *05:01 missed those carriers.
    dq25 = bool(set(dqa1) & _DQA1_DQ25) and "DQB1*02:01" in dqb1
    dq8 = "DQB1*03:02" in dqb1
    hap = " + ".join(h for h, p in [("DQ2.5", dq25), ("DQ8", dq8)] if p) or "neither"
    calls.append(f"**Celiac haplotype** (DQ2.5/DQ8): {hap} "
                 f"(DQA1 {', '.join(dqa1) or '?'}; DQB1 {', '.join(dqb1) or '?'})")
    return calls


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

    genotype_json = out_dir / f"{sample}.genotype.json"
    if force or not genotype_json.exists():
        _run_arcashla(cram, out_dir, sample)
    if not genotype_json.exists():
        raise SystemExit(f"arcasHLA did not produce {genotype_json} — check "
                         f"{config.LOGS_DIR}/hla_typing_{sample}.log")

    genotype = json.loads(genotype_json.read_text())
    for gene in config.HLA_TYPE_GENES:
        log.info("  HLA-%-6s %s", gene, ", ".join(genotype.get(gene, ["no call"])))

    rows = [(g, ", ".join(_four_digit(a) for a in genotype.get(g, [])) or "no call")
            for g in config.HLA_TYPE_GENES]

    lines = [
        f"# HLA typing (classical, arcasHLA) — {sample}",
        "",
        "| Locus | Alleles (2-field) |",
        "|-------|-------------------|",
        *[f"| HLA-{g} | {a} |" for g, a in rows],
        "",
        "## Risk calls re-derived from typed alleles",
        *[f"- {c}" for c in _risk_calls(genotype)],
        "",
        "## Notes",
        "- arcasHLA was built for RNA-seq but is widely used on WGS/exome BAMs too; it "
        "pseudoaligns chr6-extracted reads against the IMGT/HLA reference with kallisto. "
        "Homozygous loci report a single allele (no true second call to distinguish from "
        "a het that collapsed to one), the tool's known limitation, not a data gap.",
        "- Compare against `hla` (tag-SNP report) for the same sample — the risk calls "
        "here should agree; a disagreement means the tag SNP tags an uncommon haplotype "
        "for this genome and the typed allele here is the more trustworthy call.",
        "- Research-grade; clinically actionable HLA results (e.g. B27, celiac) warrant "
        "confirmatory testing.",
        "",
    ]
    write_report(summary, lines)
    log.info("HLA typing report: %s", summary)
    return summary
