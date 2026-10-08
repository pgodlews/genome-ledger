"""Stage: trio de-novo mutation calling.

A de-novo variant is present in the child and in neither parent — the one class of finding
that *needs* the whole family.

The naive implementation — intersect the child's VCF against the parents' and call whatever
is left de-novo — does not work, and fails in a direction that looks like success. Both
parents' VCFs are **variant-only**: a site absent from them is *either* genuinely
homozygous-reference *or* simply never called there. Against an expected true rate of
~60-100 de-novo SNVs genome-wide, the parental false-negatives dominate by one to two
orders of magnitude, and they dominate hardest when child and parents were called by
different callers or lifted from different reference builds (for example, when a child
is called via DeepVariant lifted from GRCh37 and parents are native GATK GRCh38).

So the stage is two passes:

    child normalized VCF ──filter──▶ het, GQ/DP/VAF-plausible, autosomal
        │  bcftools isec -C  (drop anything either parent's VCF already carries)
        ▼
    candidates  ──▶ sites VCF ──▶ GATK force-call BOTH PARENTS off their CRAMs
        │                          (--alleles, EMIT_ALL_ACTIVE_SITES: hom-ref is *evidenced*)
        ▼
    confirmed: both parents called, at depth, with no alt-read support

Pass 1 is a cheap recall-oriented filter; pass 2 is what makes the output mean anything.
A candidate the parents' reads *can't* resolve (low depth, no-call) is reported as
**indeterminate**, never as de-novo — the same honesty rule as `callable` and
`force-call`: absence of evidence is reported as such, not as evidence of absence.

Runs off the PARENTS' CRAMs, so it runs where those live (e.g. on the local machine, or
a Linux compute host where GATK's native vector PairHMM runs).
"""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

from . import reportpaths
from . import config, pedigree
from .call import aligned_cram
from .forcecall import _fai_lengths, call_at_sites
from .normalize import normalized_path
from .util import file_token, log, require_tools, run, write_report

def _file_stamp(path: Path) -> str:
    """Content identity for a de-novo cache key — util.file_token, not a second copy of it.

    This was its own `size:mtime` implementation, described as "the same idiom
    incremental.source_digests uses". That stopped being true when both of those moved to
    content hashes, leaving this the only place still keying expensive cached artifacts on a
    fingerprint that `rsync -a` and `cp -p` preserve across different bytes. The duplication
    was the reason it got left behind, so it now calls the shared helper.
    """
    return file_token(path)


def _key(*parts: str) -> str:
    h = hashlib.sha256()
    for part in parts:
        h.update(part.encode() + b"\0")
    return h.hexdigest()[:12]


def _keyed_path(wd: Path, stem: str, key: str, suffix: str) -> Path:
    """`<wd>/<stem>.<key><suffix>`, pruning siblings carrying a different key.

    Putting the key in the *filename* is what makes serving a stale result impossible — the
    convention forcecall.build_sites already uses (`sites_<hash>`). Pruning is what stops
    that guarantee from filling the disk: a carrier view is ~100 MB and there is one per
    parental callset. Also removes the pre-keying unsuffixed name, so an existing working
    directory upgrades cleanly instead of leaving a file nothing will ever read again."""
    out = wd / f"{stem}.{key}{suffix}"
    stale = [q for q in wd.glob(f"{stem}.*{suffix}") if q != out]
    stale += [q for q in (wd / f"{stem}{suffix}",) if q.exists()]
    for old_file in stale:
        old_file.unlink(missing_ok=True)
        Path(f"{old_file}.tbi").unlink(missing_ok=True)
        log.info("de-novo: pruned stale %s", old_file.name)
    return out


def _ad_counts(ad: str) -> list[int]:
    """AD ('ref,alt[,alt2…]') as per-allele integers, POSITIONALLY.

    A missing per-allele value ('.') becomes 0 in place rather than being dropped: filtering
    non-numeric fields out before indexing would shift every later allele down a slot, so
    'AD=10,.,5' would credit allele 2's reads to allele 1."""
    return [int(x) if x.lstrip("-").isdigit() else 0 for x in ad.split(",")]


def _chrom_key(chrom: str) -> tuple[int, str]:
    """Natural contig order (1…22, then X/Y/MT) rather than lexical ('10' after '9')."""
    return (int(chrom), "") if chrom.isdigit() else (10 ** 6, chrom)


# Status of one candidate after the parental force-call pass.
CONFIRMED = "de_novo"
INHERITED = "parent_carries"
INDETERMINATE = "indeterminate"
CLUSTERED = "clustered"


# --- paths -------------------------------------------------------------------
def report_path(child: str) -> Path:
    return reportpaths.genome_report("denovo", child)


def _work_dir(child: str) -> Path:
    return config.DENOVO_DIR / child


# --- pass 1: child-side candidates -------------------------------------------
def _autosomes() -> list[str]:
    """Autosome contig names as the normalized VCFs spell them (Ensembl: '1'…'22')."""
    return [str(i) for i in range(1, 23)]


def is_homopolymer_indel(ref: str, alt: str) -> bool:
    """Is this indel a single-base slip in a homopolymer run? Pure, so it is testable.

    Polymerase slippage in a homopolymer inserts or deletes copies of ONE base (CT>C, T>TA,
    G>GTT). Real de-novo indels are not restricted that way. Only the changed bases are
    examined, after stripping the shared prefix bcftools leaves on a normalized record.
    """
    if len(ref) == len(alt) or not ref or not alt:
        return False
    longer, shorter = (ref, alt) if len(ref) > len(alt) else (alt, ref)
    delta = longer[len(shorter):] if longer.startswith(shorter) else longer
    return bool(delta) and len(set(delta)) == 1


def child_passes(dp: int, gq: int, vaf: float, max_dp: int,
                 ref: str = "", alt: str = "") -> tuple[bool, str]:
    """Is one child het call a plausible de-novo candidate? Returns (ok, reason-if-not).

    Pure, so the thresholds are testable without a VCF. Each rejection maps to a specific
    false-positive mode rather than a generic quality score.

    `ref`/`alt` are optional so the depth/quality thresholds stay callable on their own; when
    they are supplied the indel-specific cuts apply as well.
    """
    if gq < config.DN_MIN_GQ:
        return False, "low_gq"
    if dp < config.DN_MIN_DP:
        return False, "low_depth"
    if max_dp and dp > max_dp:
        return False, "excess_depth"        # collapsed paralog / mismapping pile-up
    is_indel = bool(ref) and bool(alt) and (len(ref) != 1 or len(alt) != 1)
    if is_indel:
        if config.DN_DROP_HOMOPOLYMER_INDELS and is_homopolymer_indel(ref, alt):
            return False, "homopolymer_indel"   # polymerase slippage, not mutation
        if vaf < config.DN_MIN_VAF_INDEL:
            return False, "low_vaf_indel"       # indel reference bias / mismapping
    elif vaf < config.DN_MIN_VAF:
        return False, "low_vaf"             # mosaic, or reference-bias artefact
    if vaf > config.DN_MAX_VAF:
        return False, "high_vaf"            # not a balanced het
    return True, ""


def _vaf(ad: str) -> float:
    """Alt fraction from an AD field ('ref,alt[,alt2…]'). 0.0 when unparseable."""
    counts = _ad_counts(ad)
    if len(counts) < 2:
        return 0.0
    total = sum(counts)
    return (sum(counts[1:]) / total) if total else 0.0


def _median_het_depth(vcf: Path) -> int:
    """Median DP over the child's heterozygous calls — the scale for the excess-depth cut.

    Taken from the child's own genome rather than a constant: depth is a property of this
    library, and a fixed ceiling would be meaningless at 15x or at 100x."""
    out = subprocess.run(
        ["bcftools", "query", "-i", 'GT="het"', "-f", "[%DP]\n", str(vcf)],
        capture_output=True, text=True, check=True).stdout
    dps = sorted(int(x) for x in out.split() if x.isdigit())
    return dps[len(dps) // 2] if dps else 0


def parent_callsets(parent: str) -> list[str]:
    """Every ingested callset belonging to one parent, by this repo's id convention:
    the bare id plus any `<parent>_<suffix>` variant (e.g. `parent` + `parent_t2t`).

    The same person called twice by different pipelines is two rows in samples.tsv but one
    genome, and each callset's *false negatives* are largely caller-specific. Excluding
    against all of them removes exactly the candidates a single caller happened to miss —
    which is the dominant false-positive source in a cross-caller trio."""
    from .util import read_tsv
    ids = [r["sample_id"] for r in read_tsv(config.SAMPLES_TSV)]
    return [i for i in ids
            if (i == parent or i.startswith(f"{parent}_")) and normalized_path(i).exists()]


def _carrier_view(sample: str, wd: Path, force: bool = False) -> Path:
    """A parent callset filtered to records where the parent actually CARRIES an alt allele.

    Not cosmetic. `bcftools isec` matches on position/alleles and ignores genotypes, and a
    callset lifted from T2T-CHM13 can carry a large share of explicit `0/0` records (liftover residue —
    CHM13 differs from GRCh38 at millions of bases, recovered by RECOVER_SWAPPED_REF_ALT).
    Excluding on record *presence* would therefore delete real de-novo candidates at every
    such site — thousands of them in a genome."""
    out = _keyed_path(wd, f"{sample}.carrier", _key(_file_stamp(normalized_path(sample))),
                      ".vcf.gz")
    if out.exists() and out.with_suffix(".gz.tbi").exists() and not force:
        return out
    run(["bcftools", "view", "-i", 'GT="alt"', str(normalized_path(sample)),
         "-Oz", "-o", str(out)])
    run(["tabix", "-f", "-p", "vcf", str(out)])
    return out


def _private_to_child(child_vcf: Path, parent_vcfs: list[Path], out: Path) -> Path:
    """Records in the child's VCF that no parental callset carries (exact REF/ALT match).

    `bcftools isec -C -w1` — allele-aware by default, so a *different* ALT at a shared
    position is correctly treated as not-transmitted rather than as a match."""
    run(["bcftools", "isec", "-C", "-w1", str(child_vcf), *[str(p) for p in parent_vcfs],
         "-Oz", "-o", str(out)])
    run(["tabix", "-f", "-p", "vcf", str(out)])
    return out


def candidates(child: str, parents: list[str], force: bool = False) -> tuple[list[dict], dict]:
    """Pass 1: plausible de-novo candidates from the child's VCF alone. (candidates, stats)."""
    require_tools("bcftools", "tabix")
    child_vcf = normalized_path(child)
    if not child_vcf.exists():
        raise SystemExit(f"{child_vcf} missing. Run `normalize {child}` first.")
    for p in parents:
        if not normalized_path(p).exists():
            raise SystemExit(f"{normalized_path(p)} missing. Run `normalize {p}` first.")

    wd = _work_dir(child)
    wd.mkdir(parents=True, exist_ok=True)
    med = _median_het_depth(child_vcf)
    max_dp = int(config.DN_MAX_DP_FACTOR * med) if med else 0
    log.info("de-novo: child median het depth %dx → excess-depth cut at %dx.", med, max_dp)

    # Cheap prefilter in bcftools (het + GQ + DP) so isec runs over a much smaller file;
    # the VAF/excess-depth cuts need per-record arithmetic and happen in child_passes below.
    regions = ",".join(_autosomes()) if config.DN_AUTOSOMES_ONLY else ""
    # Keyed by the inputs that determine the contents: the child's VCF and every threshold
    # that goes into the bcftools filter. Tuning DN_MIN_GQ via the environment and re-running
    # must not silently report numbers computed against the previous filter.
    filt_key = _key(_file_stamp(child_vcf), str(config.DN_MIN_GQ), str(config.DN_MIN_DP),
                    str(config.DN_AUTOSOMES_ONLY))
    filt = _keyed_path(wd, f"{child}.dn_filtered", filt_key, ".vcf.gz")
    if not filt.exists() or force:
        cmd = ["bcftools", "view", "-i",
               f'GT="het" && FMT/GQ>={config.DN_MIN_GQ} && FMT/DP>={config.DN_MIN_DP}']
        if regions:
            cmd += ["-r", regions]
        cmd += [str(child_vcf), "-Oz", "-o", str(filt)]
        run(cmd)
        run(["tabix", "-f", "-p", "vcf", str(filt)])

    callsets = [c for p in parents for c in parent_callsets(p)]
    log.info("de-novo: excluding against %d parental callset(s): %s",
             len(callsets), ", ".join(callsets))
    carriers = [_carrier_view(c, wd, force) for c in callsets]
    # Keyed by the filter *and* every parental callset that went into the exclusion, so
    # ingesting another callset for a parent invalidates it rather than being ignored.
    priv_key = _key(filt_key, *(f"{c}:{_file_stamp(normalized_path(c))}" for c in callsets))
    priv = _keyed_path(wd, f"{child}.dn_private", priv_key, ".vcf.gz")
    if not priv.exists() or force:
        _private_to_child(filt, carriers, priv)

    n_filtered = int(subprocess.run(["bcftools", "index", "-n", str(filt)],
                                    capture_output=True, text=True, check=True).stdout)
    rows = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%DP\t%GQ\t%AD]\n",
         str(priv)], capture_output=True, text=True, check=True).stdout

    cands, rejected = [], {}
    for line in rows.splitlines():
        f = line.split("\t")
        if len(f) < 8:
            continue
        chrom, pos, ref, alt, _gt, dp, gq, ad = f[:8]
        dp_i = int(dp) if dp.isdigit() else 0
        gq_i = int(gq) if gq.isdigit() else 0
        vaf = _vaf(ad)
        ok, why = child_passes(dp_i, gq_i, vaf, max_dp, ref, alt)
        if not ok:
            rejected[why] = rejected.get(why, 0) + 1
            continue
        cands.append({"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt,
                      "dp": dp_i, "gq": gq_i, "vaf": vaf,
                      "kind": "SNV" if len(ref) == 1 and len(alt) == 1 else "indel"})
    cands.sort(key=lambda c: (_chrom_key(c["chrom"]), c["pos"]))
    stats = {"child_het_prefiltered": n_filtered, "private_to_child": len(rows.splitlines()),
             "callsets": callsets,
             "median_het_dp": med, "max_dp": max_dp, "rejected": rejected,
             "candidates": len(cands)}
    log.info("de-novo: %d child hets → %d private to child → %d candidates.",
             n_filtered, stats["private_to_child"], len(cands))
    return cands, stats


# --- pass 2: parental confirmation off the CRAMs ------------------------------
def _write_sites(wd: Path, cands: list[dict]) -> tuple[Path, Path, str]:
    """Sites VCF + interval list for a set of positions. Returns (vcf, intervals, key).

    Takes the working directory rather than a sample so the segregation stage can build its
    own site list with the same content-keyed guarantees."""
    lengths = _fai_lengths()
    order = {c: i for i, c in enumerate(lengths)}
    recs = sorted(((c["chrom"], c["pos"], c["ref"], c["alt"]) for c in cands),
                  key=lambda r: (order.get(r[0], len(order)), r[1]))
    unknown = sorted({c for c, *_ in recs if c not in lengths}, key=_chrom_key)
    if unknown:
        # Emitting body records for a contig with no ##contig header produces a sites VCF
        # GATK rejects — or, worse, a tool somewhere silently tolerates. The cause is always
        # a contig-naming mismatch upstream (an un-renamed 'chr1'), which is worth stopping
        # for rather than quietly dropping part of the candidate set.
        n = sum(1 for c, *_ in recs if c in set(unknown))
        raise SystemExit(
            f"de-novo: {n} candidate(s) sit on contig(s) absent from "
            f"{config.REF_FASTA.name}.fai: {', '.join(unknown)}. The child's normalized VCF "
            f"is not on the same contig naming as the reference — re-run `normalize`.")

    key = _key(*(f"{c}\t{p}\t{ref}\t{alt}" for c, p, ref, alt in recs))
    used = list(dict.fromkeys(c for c, *_ in recs))
    header = ["##fileformat=VCFv4.2",
              *[f"##contig=<ID={c},length={lengths[c]}>" for c in used],
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    body = [f"{c}\t{p}\t.\t{ref}\t{alt}\t.\t.\t." for c, p, ref, alt in recs]
    for stale in wd.glob("sites_*"):
        if not stale.name.startswith(f"sites_{key}"):
            stale.unlink(missing_ok=True)
            log.info("de-novo: pruned stale %s", stale.name)
    raw = wd / f"sites_{key}.vcf"
    vcf = wd / f"sites_{key}.vcf.gz"
    intervals = wd / f"sites_{key}.list"
    if not vcf.exists():
        raw.write_text("\n".join(header + body) + "\n")
        run(["bgzip", "-f", str(raw)])
        run(["tabix", "-f", "-p", "vcf", str(vcf)])
    intervals.write_text("".join(f"{c}:{p}-{p}\n" for c, p, *_ in recs))
    return vcf, intervals, key


def parent_verdict(gt: str, dp: int, alt_ad: int) -> tuple[str, str]:
    """What one parent's force-called genotype says about a candidate. (status, detail).

    Pure. The three outcomes are deliberately asymmetric: only an *evidenced* hom-ref at
    adequate depth with no alt-read support supports a de-novo call; anything the reads
    leave unresolved is INDETERMINATE, never a silent pass."""
    if not gt or gt.replace("|", "/").split("/")[0] in (".", ""):
        return INDETERMINATE, "no-call"
    if dp < config.DN_MIN_PARENT_DP:
        return INDETERMINATE, f"DP {dp}"
    alleles = gt.replace("|", "/").split("/")
    if any(a not in ("0", ".") for a in alleles):
        return INHERITED, f"GT {gt}"
    if alt_ad > config.DN_MAX_PARENT_ALT:
        return INHERITED, f"{alt_ad} alt reads"
    return CONFIRMED, f"0/0 at {dp}x"


def _force_call_parent(parent: str, vcf_sites: Path, intervals: Path, key: str,
                       wd: Path, force: bool = False) -> Path:
    out = _keyed_path(wd, f"{parent}.parental", key, ".vcf.gz")
    if out.exists() and out.with_suffix(".gz.tbi").exists() and not force:
        log.info("de-novo: %s already force-called at these sites — reusing.", parent)
        return out
    cram = aligned_cram(parent)
    if not cram.exists():
        raise SystemExit(
            f"{cram} missing — de-novo confirmation reads the PARENTS' CRAMs. Run this "
            f"stage where they live (or `align {parent}` first).")
    if not config.GATK_JAR.exists():
        raise SystemExit(f"GATK not installed at {config.GATK_JAR}. Run `setup`.")
    return call_at_sites(parent, cram, vcf_sites, intervals, out, key, tag="denovo")


def _read_parent(vcf: Path) -> dict[tuple[str, int], tuple[str, int, int]]:
    """(chrom, pos) → (GT, DP, alt-supporting AD) from a force-called parental VCF."""
    out = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS[\t%GT\t%DP\t%AD]\n", str(vcf)],
        capture_output=True, text=True, check=True).stdout
    geno: dict[tuple[str, int], tuple[str, int, int]] = {}
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        chrom, pos, gt, dp, ad = f[:5]
        counts = _ad_counts(ad)
        geno[(chrom, int(pos))] = (gt, int(dp) if dp.isdigit() else 0,
                                   sum(counts[1:]) if len(counts) > 1 else 0)
    return geno


def _pileup_parent(parent: str, regions: Path, key: str, wd: Path,
                   force: bool = False) -> dict[tuple[str, int], int]:
    """(chrom, pos) → non-reference read count in the parent's CRAM, from a RAW pileup.

    Deliberately caller-independent. The force-call pass asks HaplotypeCaller whether the
    parent carries the allele — but at exactly the sites that reach this stage (found by the
    child's caller, absent from the parents' HC callsets) HC's local assembly is the thing
    that failed, and it answers `0/0` at full depth with complete confidence. A pileup has no
    assembly and no genotype model: it reports the bases the reads actually carry.

    Counts *any* non-reference support (bcftools' `<*>` allele included), so a parent showing
    a different alt at the position also blocks the de-novo call. Conservative in the only
    direction that is safe here.

    Cached under the candidate-set `key`, exactly like the force-call output beside it. That
    is not tidiness: `confirm` treats a site missing from this dict as "no pileup opinion" and
    falls back to HaplotypeCaller's verdict alone, so a pileup keyed only by parent name would
    — after any change to the candidate set — silently stop covering the new candidates while
    appearing to work. The check would go missing precisely where it was needed."""
    out = _keyed_path(wd, f"{parent}.pileup", key, ".tsv")
    if not out.exists() or force:
        cmd = ["bcftools", "mpileup", "-f", str(config.REF_FASTA), "-R", str(regions),
               "-a", "AD", "-d", "200", "-q", str(config.DN_PILEUP_MIN_MQ),
               "-Q", str(config.DN_PILEUP_MIN_BQ), "--no-BAQ", str(aligned_cram(parent))]
        # Write to a temp file and rename only on success: `out.exists()` is the cache gate,
        # so a half-written pileup left behind by a failure would be served as a complete one
        # on the next run — the same silent-degradation hazard the key above prevents.
        tmp = out.with_suffix(".partial")
        errlog = config.LOGS_DIR / f"denovo_pileup_{parent}.log"
        try:
            with open(tmp, "wb") as fh, open(errlog, "wb") as lf:
                mp = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=lf)
                try:
                    subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\t[%AD]\n"],
                                   stdin=mp.stdout, stdout=fh, check=True)
                finally:
                    # Always reached, so mpileup is reaped even when query raises (it takes
                    # SIGPIPE and exits); without this the process was left dangling and its
                    # exit status never inspected.
                    mp.stdout.close()
                    rc = mp.wait()
            if rc:
                raise SystemExit(
                    f"de-novo: mpileup failed for {parent} (exit {rc}) — see {errlog}")
            tmp.replace(out)
        finally:
            tmp.unlink(missing_ok=True)
    geno: dict[tuple[str, int], int] = {}
    for line in out.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 3:
            continue
        counts = _ad_counts(f[2])
        geno[(f[0], int(f[1]))] = sum(counts[1:]) if len(counts) > 1 else 0
    return geno


def confirm(child: str, cands: list[dict], parents: list[str],
            force: bool = False) -> list[dict]:
    """Pass 2: force-call each parent at the candidate sites and classify every candidate."""
    require_tools("bcftools", "tabix", "bgzip", "samtools")
    wd = _work_dir(child)
    vcf_sites, intervals, key = _write_sites(wd, cands)
    parental = {p: _read_parent(_force_call_parent(p, vcf_sites, intervals, key, wd, force))
                for p in parents}

    # Second pass over the SAME sites with a raw pileup — see _pileup_parent.
    regions = wd / f"sites_{key}.regions.txt"     # .txt, not .bed: 2-column chrom/pos form
    regions.write_text("".join(f"{c['chrom']}\t{c['pos']}\n" for c in cands))
    pileups = {p: _pileup_parent(p, regions, key, wd, force) for p in parents}

    for c in cands:
        verdicts = {}
        for p in parents:
            gt, dp, alt_ad = parental[p].get((c["chrom"], c["pos"]), ("", 0, 0))
            verdicts[p] = parent_verdict(gt, dp, alt_ad)
            nonref = pileups[p].get((c["chrom"], c["pos"]))
            if nonref is not None and nonref >= config.DN_PILEUP_MIN_ALT:
                verdicts[p] = (INHERITED, f"{nonref} non-ref reads (pileup)")
        c["parents"] = verdicts
        statuses = {v[0] for v in verdicts.values()}
        # Any parent carrying it settles the question; otherwise every parent must have
        # given positive hom-ref evidence before we call it de-novo.
        c["status"] = (INHERITED if INHERITED in statuses
                       else INDETERMINATE if INDETERMINATE in statuses
                       else CONFIRMED)
    _flag_clusters(cands)
    return cands


def _flag_clusters(cands: list[dict]) -> int:
    """Demote confirmed calls that sit within DN_CLUSTER_WINDOW of another confirmed call.

    Independent point mutations do not arrive in pairs a few hundred bases apart; a cluster
    is one unreliable locus firing repeatedly. Applied to the confirmed set only — a
    candidate already explained as inherited needs no further explanation."""
    conf = sorted((c for c in cands if c["status"] == CONFIRMED),
                  key=lambda c: (_chrom_key(c["chrom"]), c["pos"]))
    flagged = set()
    for i, c in enumerate(conf):
        for other in conf[i + 1:]:
            if other["chrom"] != c["chrom"] or other["pos"] - c["pos"] > config.DN_CLUSTER_WINDOW:
                break
            flagged.add(id(c))
            flagged.add(id(other))
    for c in conf:
        if id(c) in flagged:
            c["status"] = CLUSTERED
    return len(flagged)


# --- annotation (best-effort) -------------------------------------------------
def _annotate(child: str, cands: list[dict]) -> None:
    """Attach gene/consequence/ClinVar from DuckDB where the child's variants are loaded.

    Best-effort: the stage runs off VCFs and CRAMs, so a child who hasn't been annotated
    yet still gets a full de-novo report — just without gene names."""
    if not cands:
        return
    try:
        from .load import TABLE, connect
        from .snapshot import latest_snapshot
        snap = latest_snapshot()
        if not snap:
            return
        from .panel import sql_in_list
        con = connect()
        # Filter by position IN SQL. Without it this grouped every row the child has in the
        # snapshot (millions) and threw all but the candidate positions away in Python — ~1.8 GB
        # of peak RSS to retrieve a few thousand annotations. Positions are ints parsed from
        # our own VCF, so inlining them is safe; chroms go through the shared quoting helper.
        chroms = sql_in_list({c["chrom"] for c in cands})
        positions = ", ".join(str(p) for p in sorted({c["pos"] for c in cands}))
        rows = con.execute(
            f"SELECT chrom, pos, any_value(gene), any_value(consequence), "
            f"any_value(clinvar_sig) FROM {TABLE} WHERE snapshot_id = ? AND sample_id = ? "
            f"AND chrom IN ({chroms}) AND pos IN ({positions}) "
            f"GROUP BY chrom, pos", [snap, child]).fetchall()
        loaded = con.execute(
            f"SELECT count(*) FROM {TABLE} WHERE snapshot_id = ? AND sample_id = ?",
            [snap, child]).fetchone()[0]
        con.close()
        if not loaded:
            # Silent before: the query simply matched nothing and every row lost its gene
            # name with no indication why. The usual cause is a child annotated against an
            # older snapshot than the latest.
            log.warning("de-novo: %s has no variants loaded for the latest snapshot %s — "
                        "reporting positions without gene annotation. Run `annotate %s %s` "
                        "and `load %s %s` to annotate the report.",
                        child, snap, child, snap, child, snap)
            return
        ann = {(r[0], r[1]): r[2:] for r in rows}
        for c in cands:
            gene, csq, sig = ann.get((c["chrom"], c["pos"]), (None, None, None))
            c["gene"], c["csq"], c["clinvar"] = gene, csq, sig
    except Exception as e:  # noqa: BLE001 — annotation is a nicety, never the point
        log.warning("de-novo: gene annotation unavailable (%s) — reporting positions only.", e)


# --- report -------------------------------------------------------------------
def _row(c: str, parents: list[str], cand: dict) -> str:
    pv = cand.get("parents", {})
    detail = " · ".join(f"{p}: {pv[p][1]}" for p in parents if p in pv)
    gene = cand.get("gene") or "—"
    csq = (cand.get("csq") or "—").split(",")[0]
    return (f"| {cand['chrom']}:{cand['pos']} | {cand['ref']}>{cand['alt']} | {cand['kind']} | "
            f"{cand['dp']}× | {cand['vaf']:.2f} | {gene} | {csq} | "
            f"{cand.get('clinvar') or '—'} | {detail} |")


def _filter_lines(stats: dict) -> list[str]:
    """The "how the candidates were filtered" bullets — shared by the confirmed and the
    `--no-confirm` report, which must describe the same pass-1 filter."""
    rej = stats.get("rejected", {})
    return [
        "",
        f"- Child heterozygous calls passing GQ≥{config.DN_MIN_GQ}, DP≥{config.DN_MIN_DP}"
        f"{', autosomes only' if config.DN_AUTOSOMES_ONLY else ''}: "
        f"**{stats['child_het_prefiltered']:,}**",
        f"- …of which carried by none of the {len(stats.get('callsets', []))} parental "
        f"callset(s) `{'`, `'.join(stats.get('callsets', []))}` (`bcftools isec -C`, "
        f"carrier records only): **{stats['private_to_child']:,}**",
        f"- …surviving the per-record cuts: **{stats['candidates']:,}** "
        f"(rejected: " + (", ".join(f"{k} {v:,}" for k, v in sorted(rej.items())) or "none")
        + ")",
        f"- Child median het depth **{stats['median_het_dp']}×** → excess-depth cut at "
        f"**{stats['max_dp']}×** (collapsed-paralog guard)",
        f"- Balanced-het window: VAF {config.DN_MIN_VAF}–{config.DN_MAX_VAF} "
        f"(indels {config.DN_MIN_VAF_INDEL}–{config.DN_MAX_VAF}"
        + ("; single-base homopolymer slips dropped)" if config.DN_DROP_HOMOPOLYMER_INDELS
           else ")"),
    ]


def _caveat_lines(confirmed: bool = True) -> list[str]:
    """Method & caveats — shared, so the `--no-confirm` report carries the same limits.

    `confirmed` switches the bullets that assert the parental CRAM pass actually ran; stating
    it did in a `--no-confirm` report would misdescribe exactly the check that is missing."""
    parental = ([
        "- **Absence from a parent's VCF is not evidence.** Both parents' VCFs are "
        "variant-only, so a site missing from them is either hom-ref or never called. "
        "Every candidate here was re-genotyped **off the parents' CRAMs** with GATK "
        "`--alleles … EMIT_ALL_ACTIVE_SITES`, which emits an *evidenced* `0/0` — that "
        "confirmation, not the VCF intersection, is what the de-novo count rests on.",
        f"- A parent needs DP≥{config.DN_MIN_PARENT_DP} and ≤{config.DN_MAX_PARENT_ALT} "
        "alt-supporting read(s) to count as hom-ref. A single alt read at high depth is "
        "indistinguishable from sequencing error, so it is not treated as transmission.",
        "- **The parental check runs twice, on purpose.** The force-call asks GATK "
        "HaplotypeCaller; but at sites the child's caller found and HC systematically "
        "cannot, HC returns a confident *wrong* hom-ref — it is not an independent test "
        "there. So every candidate is also checked against a **raw pileup** (no assembly, "
        "no genotype model): ≥"
        f"{config.DN_PILEUP_MIN_ALT} non-reference reads in a parent marks it inherited "
        "regardless of what any caller said.",
    ] if confirmed else [
        "- **Absence from a parent's VCF is not evidence, and that is all this report has.** "
        "Both parents' VCFs are variant-only, so a site missing from them is either hom-ref "
        "or never called — and the pass that tells those apart (re-genotyping the parents' "
        "CRAMs, plus a caller-independent pileup) was skipped. Run without `--no-confirm` "
        "before reading any row as a finding.",
    ])
    return parental + [
        f"- **Expected scale:** roughly {config.DN_EXPECTED_MIN}–{config.DN_EXPECTED_MAX} "
        "de-novo SNVs genome-wide per generation. A count far above that means the filters "
        "are leaking (mismapping, parental dropout, cross-caller artefacts); far below means "
        "they are too strict.",
        "- **Hom-alt calls are excluded by construction** — two independent de-novo hits at "
        "one position is vanishingly rare, so `1/1` indicates parental miscalling instead.",
        "- **Post-zygotic mosaicism** presents below the VAF floor and is filtered out here; "
        "this stage measures germline de-novo events only.",
        "- Cross-caller/cross-build trios (child DeepVariant-lifted vs parents native GATK) "
        "carry extra artefact load — the indeterminate column is where that surfaces. Each "
        "parent is excluded against **every** callset ingested for them, so a candidate that "
        "only one caller missed is filtered before the (expensive) CRAM pass.",
        "- Research-grade, not a clinical diagnosis. A confirmed de-novo variant in a "
        "disease gene warrants orthogonal validation (Sanger) before it means anything.",
    ]


def _report(child: str, parents: list[str], cands: list[dict], stats: dict,
            confirmed_pass: bool) -> Path:
    conf = [c for c in cands if c.get("status") == CONFIRMED]
    indet = [c for c in cands if c.get("status") == INDETERMINATE]
    inh = [c for c in cands if c.get("status") == INHERITED]
    clust = [c for c in cands if c.get("status") == CLUSTERED]
    snv = sum(1 for c in conf if c["kind"] == "SNV")

    hdr = ["| Position | Change | Type | Child DP | Child VAF | Gene | Consequence | "
           "ClinVar | Parental evidence |",
           "|----------|--------|------|---------:|----------:|------|-------------|"
           "---------|-------------------|"]
    head = [
        f"# Trio de-novo variants — {child}",
        "",
        f"- **Trio:** child `{child}`  ·  parents {', '.join(f'`{p}`' for p in parents)}",
    ]
    if confirmed_pass:
        lines = head + [
            f"- **Confirmed de-novo:** **{len(conf)}** ({snv} SNV, {len(conf) - snv} indel)",
            f"- **Indeterminate** (parental reads couldn't resolve): {len(indet)}",
            f"- **Inherited after all** (parent carries it, missed by their VCF): {len(inh)}",
            f"- **Dropped as clustered** (within {config.DN_CLUSTER_WINDOW:,} bp of another "
            f"call — one bad locus, not independent mutations): {len(clust)}",
            "",
        ]
    else:
        # Never headline a confirmed/inherited/clustered tally here: without the parental
        # pass every candidate is unclassified, so those counts are all zero and reporting
        # them reads as "nothing found" rather than "nothing checked".
        lines = head + [
            f"- **Candidates (unconfirmed):** **{len(cands):,}** "
            f"({sum(1 for c in cands if c['kind'] == 'SNV'):,} SNV, "
            f"{sum(1 for c in cands if c['kind'] == 'indel'):,} indel)",
            "",
        ]
    # Plausibility gate. The germline de-novo rate is one of the few numbers in this whole
    # pipeline with a firm published prior, which makes it a usable self-check: a count far
    # above it is not a discovery, it is a leak. Saying so in the headline — rather than
    # burying it in the caveats under a confident-looking number — is the point.
    # Compare LIKE WITH LIKE: the 40-150 prior is a de-novo *SNV* rate, so the gate keys on the
    # confirmed SNV count. Testing the SNV+indel total against it overstated the excess: an SNV
    # count just above the range read as a large breach once the indels were added.
    n_snv = sum(1 for c in conf if c["kind"] == "SNV")
    n_indel = len(conf) - n_snv
    if confirmed_pass and n_snv > config.DN_EXPECTED_MAX:
        ratio = n_snv / config.DN_EXPECTED_MAX
        lines += [
            f"> **This count is not a de-novo call set.** {n_snv:,} confirmed SNVs against an "
            f"expected {config.DN_EXPECTED_MIN}–{config.DN_EXPECTED_MAX} germline de-novo "
            f"SNVs per generation — **{ratio:.1f}× the top of that range** "
            f"(plus {n_indel:,} indels, which have no comparable published prior here). The "
            "parental side has been checked two independent ways (force-call + raw pileup), so "
            "the residual is on the **child's** side: artefacts of the child's own callset that "
            "cannot be separated from real de-novo events without piling up the child's own "
            "**reads**. Treat the rows below as a candidate list, not as findings.",
            "",
        ]
    elif confirmed_pass and n_snv < config.DN_EXPECTED_MIN:
        lines += [
            f"> **Fewer than expected.** {n_snv} confirmed SNVs against an expected "
            f"{config.DN_EXPECTED_MIN}–{config.DN_EXPECTED_MAX}; the filters may be too "
            "strict, or parental coverage may be masking real events (see indeterminate).",
            "",
        ]
    if not confirmed_pass:
        lines += [
            "> **UNCONFIRMED — candidates only.** The parental pass (`--no-confirm`) did not "
            "run, so nothing below has been checked against the parents' reads. Unconfirmed "
            "candidates are typically dominated by parental false-negatives rather than real "
            "de-novo events; most are reclassified once the parents' reads are checked. Do not "
            "read the count as a de-novo rate.",
            "",
            "## Candidates (unconfirmed)",
            "",
            *hdr,
            *([_row(child, parents, c) for c in cands[:200]] or
              ["| _none_ | | | | | | | | |"]),
            *([f"", f"_…and {len(cands) - 200:,} more (full list is the candidate VCF in "
                    f"`{_work_dir(child)}`)._"] if len(cands) > 200 else []),
            "",
            "## How the candidates were filtered",
        ]
        lines += _filter_lines(stats)
        lines += ["", "## Method & caveats", ""] + _caveat_lines(confirmed=False)
        out = report_path(child)
        out.parent.mkdir(parents=True, exist_ok=True)
        write_report(out, lines)
        return out

    lines += [
        "## Confirmed de-novo",
        "",
        *hdr,
        *([_row(child, parents, c) for c in conf] or
          ["| _none_ | | | | | | | | |"]),
        "",
        "## Indeterminate — parental reads inadequate",
        "",
        "Not de-novo and not inherited: the parents' reads could not resolve the site "
        "(no-call, or below the depth floor). Reported rather than dropped — the same rule "
        "the callability report applies to the genome at large.",
        "",
        *hdr,
        *([_row(child, parents, c) for c in indet] or
          ["| _none_ | | | | | | | | |"]),
        "",
        "## Dropped as clustered",
        "",
        f"Confirmed against both parents, but within {config.DN_CLUSTER_WINDOW:,} bp of "
        "another such call. Independent de-novo point mutations do not arrive in bursts; a "
        "cluster means one locus is mismapping or collapsed, so all of its calls are "
        "suspect. Listed rather than deleted — a genuine multi-nucleotide event would also "
        "land here.",
        "",
        *hdr,
        *([_row(child, parents, c) for c in clust[:100]] or
          ["| _none_ | | | | | | | | |"]),
        *([f"", f"_…and {len(clust) - 100:,} more._"] if len(clust) > 100 else []),
        "",
        "## How the candidates were filtered",
        *_filter_lines(stats),
        "",
        "## Method & caveats",
        "",
        *_caveat_lines(),
        "",
    ]
    out = report_path(child)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_report(out, lines)
    return out


# --- entry point --------------------------------------------------------------
def denovo(child: str, force: bool = False, confirm_parents: bool = True) -> Path:
    """Full trio de-novo pass for `child`, resolving the trio from the family graph."""
    fam = pedigree.load(sequenced=pedigree.sequenced_ids())
    trio = fam.trio(child)
    if not trio:
        raise SystemExit(
            f"No trio for {child}: both parents must be set in {config.FAMILY_FILE} and "
            f"sequenced. Edit the family graph, then re-run.")
    _child, father, mother = trio
    parents = [father, mother]
    log.info("de-novo: trio %s ← (%s, %s)", child, father, mother)

    out = report_path(child)
    if out.exists() and not force:
        log.info("%s already has a de-novo report (%s) — skipping (use --force).", child, out)
        return out

    cands, stats = candidates(child, parents, force=force)
    if confirm_parents and cands:
        cands = confirm(child, cands, parents, force=force)
    _annotate(child, cands)
    # Whether the parental pass RAN — not whether it found anything. Gating on `cands` made
    # a zero-candidate run print "the parental force-call pass did not run", which is false
    # and reads like a failure.
    report = _report(child, parents, cands, stats, confirmed_pass=confirm_parents)
    n = sum(1 for c in cands if c.get("status") == CONFIRMED)
    log.info("de-novo report: %s (%d candidates, %d confirmed)", report, len(cands), n)
    return report
