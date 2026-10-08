"""Stage: polygenic risk scores.

Computes a raw weighted dosage score from PGS Catalog **harmonized** (GRCh38) scoring
files: score = Σ (dosage of effect allele) × weight, over the score's variants. The sample's
genotype at each position is read from the called VCF; positions absent from the
variant-only VCF are taken as homozygous reference (their effect-allele dosage inferred
from the reference base) — the same high-coverage-WGS assumption the other engines make.

Percentile calibration places the raw score on a population distribution: 1000 Genomes
(NYGC, GRCh38) carries per-superpopulation allele frequencies in its INFO, so we remote-query
just the score's positions and model the population score as Normal(Σ 2·f·w, Σ 2·f(1−f)·w²)
under HWE + CLT — an analytic stand-in for scoring all 3202 reference individuals. The
percentile is reported ancestry-matched (EUR) and overall, and cached per PGS.

Honest limits (documented in the report):
  - The percentile is *relative genetic load*, not absolute risk; most PGS are
    European-derived and miscalibrate across ancestries (hence the EUR-matched column).
  - Per-variant lookup suits small scores; genome-wide scores (millions of variants) use
    **plink2 --score**, with an *empirical* percentile from scoring the same PGS over the
    1000G reference panel (restricted to the sample's covered variants so the comparison is
    fair). Rigorous absolute dosage still needs the force-call step (P3/P4 prereq).
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from . import reportpaths
from . import config, consent, forcecall
from .util import (PASS_ONLY, faidx_bases, file_token, log, require_tools,
                   stamp_is_current,
                   write_report, write_stamp)


def _called_vcf(sample: str) -> Path:
    """This sample's VCF — the pipeline's own call, or the normalized ingested one."""
    from .call import sample_vcf

    return sample_vcf(sample)


def _scoring_file(pid: str) -> Path:
    """The harmonized PGS scoring file, downloaded once and cached.

    Atomic and bounded. urlretrieve wrote straight to the final name with no timeout, so an
    interrupted or stalled download left a truncated .txt.gz that every later run treated as
    already present — and a truncated gzip either raises deep inside _load_variants or, worse,
    parses to a short variant list that scores as a real result.
    """
    dest = config.PGS_DIR / f"{pid}_hmPOS_GRCh38.txt.gz"
    if dest.exists():
        return dest
    config.PGS_DIR.mkdir(parents=True, exist_ok=True)
    log.info("Downloading PGS scoring file %s …", pid)
    tmp = dest.with_name(f"{dest.name}.{os.getpid()}.part")
    try:
        with urllib.request.urlopen(config.PGS_URL_TMPL.format(pid=pid),
                                    timeout=config.PGS_DOWNLOAD_TIMEOUT) as r, \
                open(tmp, "wb") as fh:
            shutil.copyfileobj(r, fh)
        # Verify it is a readable gzip before publishing it under the cached name: a
        # truncated body that still transferred "successfully" is the case this guards.
        with gzip.open(tmp, "rb") as fh:
            fh.read(1 << 16)
        os.replace(tmp, dest)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    return dest


def _load_variants(path: Path) -> list[tuple[str, str, str, str, float]]:
    """(chrom, pos, effect_allele, other_allele, weight) from a harmonized file."""
    out, header = [], None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = {c: i for i, c in enumerate(f)}
                continue
            try:
                chrom, pos = f[header["hm_chr"]], f[header["hm_pos"]]
                ea = f[header["effect_allele"]]
                # Not every harmonized file carries `other_allele`; many (PGS000024,
                # PGS000040, PGS000149, …) publish only the harmonizer's inferred
                # `hm_inferOtherAllele`. Requiring the former dropped *every* row of such a
                # file, and the score then reported a clean-looking 0.0 over 0 variants —
                # silent and wrong. Fall back, and let the 0-variant guard below catch the
                # genuinely unusable ones.
                oa_col = "other_allele" if "other_allele" in header else "hm_inferOtherAllele"
                oa = f[header[oa_col]]
                w = float(f[header["effect_weight"]])
            except (KeyError, ValueError, IndexError):
                continue
            if chrom and pos:
                out.append((chrom, pos, ea, oa, w))
    return out


_PASS_ONLY = PASS_ONLY   # shared with traits/HLA — see util.PASS_ONLY


def _present_genotypes(vcf: Path, variants) -> dict[tuple[str, str], tuple[str, str, str]]:
    """(chrom,pos) → (REF, ALT, GT) for the score's PASS-filtered positions present in the VCF."""
    regions = "".join(f"{c}\t{p}\n" for c, p, *_ in variants)
    # Per-process name: PGS_DIR is shared across samples, so a fixed name would let two
    # concurrent runs clobber each other's region file.
    # PGS_DIR is not in config.ALL_DIRS, so `setup` never creates it — the other writers here
    # each mkdir on demand and these two were missed. In a full prs run _scoring_file happens
    # to create it first, but a caller that reaches _score_one another way would die on the
    # region file on a machine that has never downloaded a scoring file.
    config.PGS_DIR.mkdir(parents=True, exist_ok=True)
    reg = config.PGS_DIR / f".regions.{os.getpid()}.tmp"
    reg.write_text(regions)
    try:
        proc = subprocess.run(
            ["bcftools", "query", "-R", str(reg), *_PASS_ONLY,
             "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
             str(vcf)], capture_output=True, text=True)
    finally:
        reg.unlink(missing_ok=True)
    # A failed query used to be indistinguishable from "the sample carries none of these":
    # empty stdout -> every position absent -> _ref_bases supplies the reference base ->
    # scored as homozygous reference throughout. A missing index or an unreadable VCF then
    # produced a confident low percentile instead of an error, and _score_entry discards
    # n_carry so nothing downstream noticed. The scoring-file branch already refuses to
    # score what it could not parse ("the worst possible outcome for a health report");
    # the genotype query needs the same treatment. Raising here is caught by prs_report's
    # per-score isolation, so the trait lands in the "Scores that failed" table.
    if proc.returncode:
        raise RuntimeError(
            f"bcftools query failed on {vcf.name} (rc={proc.returncode}): "
            f"{' '.join(proc.stderr.split())[:300]}")
    out = proc.stdout
    present = {}
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 5:   # malformed/truncated bcftools output line — skip, don't crash
            continue
        c, p, ref, alt, gt = f[:5]
        present[(c, p)] = (ref, alt, gt)
    return present


def _ref_bases(absent) -> dict[tuple[str, str], str]:
    """Reference base at each absent position (→ the sample is hom-ref there). Batched
    (see util.faidx_bases / util._FAIDX_BATCH) — a large marker panel blows past the OS
    exec() argument-list limit in one unbatched shot."""
    return faidx_bases(absent)


def _dosage(effect: str, alleles: list[str]) -> int:
    return sum(1 for a in alleles if a == effect)


def _score_one(vcfs, variants) -> tuple[float, int, int, dict]:
    """Returns (raw score, n_used, n_carrying, contrib) where contrib maps
    (chrom,pos,effect,other) → dosage×weight for every scored variant — used to recompute the score over exactly the
    subset that has a 1000G frequency, so the percentile compares like with like.

    `vcfs` is a fallback chain consulted in order. The first VCF holding a record at a
    position wins — including a no-call, which is evidence that the site could not be
    genotyped and must stay out of the score. Only a position no VCF carries at all falls
    through to the reference base. Ordering the force-called VCF ahead of the called VCF
    (rather than substituting it) is what keeps curated sites evidenced without turning
    every marker outside the curated union into an assumed hom-ref.

    Contributions are keyed by the scoring row (chrom, pos, effect, other), not by position:
    a scoring file may carry two rows at one coordinate (two alleles of a multiallelic site),
    and a position key let the second overwrite the first in `contrib` while `total` summed
    both — so the headline raw score and the percentile's recomputed score disagreed.
    """
    present: dict[tuple[str, str], tuple[str, str, str]] = {}
    remaining = list(variants)
    for vcf in vcfs:
        if not remaining:
            break
        found = _present_genotypes(vcf, remaining)
        present.update(found)
        remaining = [v for v in remaining if (v[0], v[1]) not in found]
    absent = [(c, p) for c, p, *_ in variants if (c, p) not in present]
    refbases = _ref_bases(absent)
    total, n_carry, n_used = 0.0, 0, 0
    contrib: dict[tuple[str, str, str, str], float] = {}
    for c, p, ea, oa, w in variants:
        if (c, p) in present:
            ref, alt, gt = present[(c, p)]
            coding = [ref, *alt.split(",")]
            fields = gt.replace("|", "/").split("/")
            called = [g for g in fields if g not in (".", "")]
            if not called:
                continue                      # "./." — nothing was called here
            if len(called) < len(fields):
                # A PARTIAL no-call ("1/."): one allele is known, the other was not called.
                # The doubling below is for haploid genotypes, and applying it here invented
                # a homozygote — "1/." scored two effect alleles, indistinguishable from a
                # confirmed "1/1". The site is unresolved, so it stays out of the score
                # exactly as "./." does rather than being guessed at in the risk-raising
                # direction.
                continue
            try:
                alleles = [coding[int(g)] for g in called]
            except (ValueError, IndexError):
                continue
            if len(alleles) == 1:
                # A genuinely haploid call — ONE GT field, not a half-missing diploid one:
                # chrX/chrY in a male, or chrM. Doubling puts it on the same 0-2 dosage
                # scale as an autosomal genotype.
                alleles *= 2
            d = _dosage(ea, alleles)
        else:
            rb = refbases.get((c, p))
            if rb is None:
                continue
            d = _dosage(ea, [rb, rb])
        # Counted for BOTH branches. n_carry used to be incremented only for sites present
        # in the VCF, while `total` also summed the reference-inferred ones — so wherever the
        # effect allele IS the reference allele (absent ⇒ two copies), the headline count and
        # the "sites with at least one copy" figure disagreed with each other.
        if d:
            n_carry += 1
        total += d * w
        contrib[(c, p, ea, oa)] = d * w
        n_used += 1
    return total, n_used, n_carry, contrib


# --- percentile calibration against 1000 Genomes ----------------------------
def _kg_cache(pid: str) -> Path:
    return config.PGS_DIR / f"{pid}_1kg_af.tsv"


def _reference_freqs(pid: str, variants) -> dict[tuple[str, str], tuple]:
    """(chrom,pos) → (ref, alt, af_eur, af_all) for the score's positions, from 1000G.
    Cached per PGS; first call remote-queries the per-chromosome VCFs at just these
    positions. Returns {} if the lookup is unavailable (calibration is then skipped)."""
    cache = _kg_cache(pid)
    # An empty cache is treated as absent, not as "this score has no reference frequencies".
    # A zero-byte file here is always a bug's residue (see the write guard below), and read
    # back as authoritative it silently suppresses calibration on every future run with no
    # invalidation path short of deleting the file by hand.
    if cache.exists() and cache.stat().st_size > 0:
        out = {}
        for line in cache.read_text().splitlines():
            f = line.split("\t")
            if len(f) != 6:
                continue
            out[(f[0], f[1])] = (f[2], f[3], f[4], f[5])
        return out

    by_chrom: dict[str, list[str]] = {}
    for c, p, *_ in variants:
        by_chrom.setdefault(c, []).append(p)
    freqs: dict[tuple[str, str], tuple] = {}
    incomplete: list[str] = []  # chromosomes queried but that yielded zero rows — a
    # transient remote-1000G hiccup, not "this score has no variants on that contig".
    failed: list[str] = []      # chromosomes whose query REPORTED an error. Separate from
    # `incomplete` on purpose: bcftools can stream rows and then fail (a truncated remote
    # read, a timeout partway), so a non-zero exit that still produced output looked complete
    # to the zero-rows check and got cached. Every later run then read that truncated cache
    # as authoritative, with no invalidation path short of deleting the file by hand.
    config.PGS_DIR.mkdir(parents=True, exist_ok=True)
    for chrom, positions in sorted(by_chrom.items()):
        reg = config.PGS_DIR / f".kg_reg_{chrom}.{os.getpid()}.tmp"
        reg.write_text("".join(f"chr{chrom}\t{p}\n" for p in positions))
        url = config.kg_vcf_source(chrom)
        log.info("  1000G AF query chr%s (%d positions) …", chrom, len(positions))
        try:
            proc = subprocess.run(
                ["bcftools", "query", "-R", str(reg),
                 "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF_EUR\t%INFO/AF\n", url],
                capture_output=True, text=True, timeout=900)
            out = proc.stdout
            # bcftools reports real failures (missing INFO tag, unreadable index, bad region
            # file) on stderr while still exiting 0 with empty stdout. Discarding it turned
            # every such failure into an unexplained "no rows for chrN" warning.
            if proc.returncode or (not out.strip() and proc.stderr.strip()):
                failed.append(chrom)
                log.warning("  bcftools chr%s rc=%s: %s", chrom, proc.returncode,
                            " ".join(proc.stderr.split())[:300])
        except Exception as e:  # noqa: BLE001 — calibration is best-effort
            log.warning("  1000G query failed for chr%s (%s) — calibration skipped.", chrom, e)
            return {}
        finally:
            reg.unlink(missing_ok=True)
        n_before = len(freqs)
        for line in out.splitlines():
            f = line.split("\t")
            if len(f) < 6:
                continue
            ch = f[0][3:] if f[0].startswith("chr") else f[0]
            freqs[(ch, f[1])] = (f[2], f[3], f[4], f[5])
        if len(freqs) == n_before:
            incomplete.append(chrom)
    if failed:
        # Rows from a failed query are still usable for THIS run (the percentile is computed
        # over whatever subset has frequencies) but must never be persisted as the score's
        # reference set.
        log.warning("  1000G AF query reported an error for chr%s — cache NOT written "
                    "(a non-zero exit can still have emitted rows).", ", chr".join(failed))
        return freqs
    if incomplete:
        # Don't cache a partial result as if complete — every later run would read the
        # truncated cache and silently calibrate against a subset of the score's variants,
        # with no invalidation path short of deleting the file by hand.
        log.warning("  1000G AF lookup returned no rows for chr%s — cache NOT written "
                    "(would silently stick a partial calibration).", ", chr".join(incomplete))
        return freqs
    if not freqs:
        # Nothing was found — including the case where the caller passed no variants at all
        # (a scoring file that parsed to zero rows). The per-chromosome `incomplete` check
        # above cannot catch that: with no chromosomes to query the loop never runs, so
        # without this guard we would write a zero-byte cache and permanently disable
        # calibration for this score.
        log.warning("  no 1000G frequencies resolved for %s — cache NOT written.", pid)
        return {}
    tmp = cache.with_name(f"{cache.name}.{os.getpid()}.part")   # atomic, as for scores_all
    tmp.write_text("".join(f"{c}\t{p}\t{r}\t{a}\t{e}\t{l}\n"
                           for (c, p), (r, a, e, l) in freqs.items()))
    os.replace(tmp, cache)
    return freqs


def _effect_freq(ea: str, oa: str, rec: tuple) -> tuple[float, float] | None:
    """(eur_freq, all_freq) of the effect allele, oriented against the 1000G REF/ALT."""
    ref, alt, afe, afa = rec
    alts, fe, fa = alt.split(","), afe.split(","), afa.split(",")
    try:
        if ea in alts:
            i = alts.index(ea)
            return float(fe[i]), float(fa[i])
        if ea == ref:                      # effect allele is REF → 1 − Σ(alt freqs)
            se = sum(float(x) for x in fe if x not in (".", ""))
            sa = sum(float(x) for x in fa if x not in (".", ""))
            return max(0.0, 1 - se), max(0.0, 1 - sa)
    except (ValueError, IndexError):
        return None
    return None


def _norm_cdf(z: float) -> float:
    return 0.5 * (1 + math.erf(z / math.sqrt(2)))


def _calibrate(contrib, variants, freqs) -> tuple[float | None, float | None, int]:
    """Percentile of the sample's score in the EUR and overall population distributions,
    each modelled as Normal(Σ 2·f·w, Σ 2·f(1−f)·w²) under Hardy-Weinberg + CLT. The sample
    score and the distribution moments are summed over the *same* variant subset — those
    present in `contrib` (the sample was scored) AND in `freqs` (a 1000G frequency exists)."""
    raw = me = ve = ma = va = 0.0
    n = 0
    for c, p, ea, oa, w in variants:
        rec = freqs.get((c, p))
        key = (c, p, ea, oa)          # per scoring ROW — see _score_one on duplicate coords
        if not rec or key not in contrib:
            continue
        ef = _effect_freq(ea, oa, rec)
        if ef is None:
            continue
        fe, fa = ef
        raw += contrib[key]
        me += 2 * fe * w
        ve += 2 * fe * (1 - fe) * w * w
        ma += 2 * fa * w
        va += 2 * fa * (1 - fa) * w * w
        n += 1

    def pct(mean, var):
        return _norm_cdf((raw - mean) / math.sqrt(var)) * 100 if var > 0 else None
    return pct(me, ve), pct(ma, va), n


# --- genome-wide scoring via plink2 ------------------------------------------
def _plink2_bin() -> str:
    hits = list(config.PLINK2_DIR.glob("plink2"))
    b = str(hits[0]) if hits else shutil.which("plink2")
    if not b:
        raise SystemExit(f"plink2 not found (under {config.PLINK2_DIR} or PATH). Run `setup`.")
    return b


# Bump when the pgen build recipe changes — the stamp then invalidates every cached pgen.
# v2: added --var-filter (non-PASS records were being imported as real genotypes).
_PGEN_RECIPE = "v2-varfilter"


def _pgen_token(sample: str) -> str:
    return f"{_PGEN_RECIPE}|{file_token(_called_vcf(sample))}"


def _pgen(sample: str) -> Path:
    """Build (cached) a plink2 pgen from the called VCF; variant IDs are chr:pos, duplicate
    positions dropped. Reused across scores and runs (the import dominates plink2 runtime).

    The cache is keyed on the source VCF's identity and the build recipe, not just the sample
    name: re-calling or re-ingesting a sample used to leave PRS silently scoring
    the previous genome's genotypes."""
    prefix = config.PLINK2_DIR / sample
    token = _pgen_token(sample)
    if prefix.with_suffix(".pgen").exists() and stamp_is_current(prefix, token):
        return prefix
    if prefix.with_suffix(".pgen").exists():
        log.info("%s pgen is stale (VCF or build recipe changed) — rebuilding.", sample)
    config.PLINK2_DIR.mkdir(parents=True, exist_ok=True)
    log.info("Building plink2 pgen for %s (one-time) …", sample)
    # --vcf-allow-no-nonvar: our VCF is variant-sites-only (every GT carries an ALT), which
    # plink2 otherwise rejects as malformed; REF/ALT come correctly from the caller.
    # --autosome: restrict to chr1-22 (chrX needs sex/ploidy handling; PGS are autosomal).
    # --var-filter (no argument) keeps only records whose FILTER is PASS or '.' — the same
    # gate _PASS_ONLY applies to the direct genotype reads. Without it the pgen (and every
    # consumer of it) carries the caller's rejected calls.
    subprocess.run([_plink2_bin(), "--vcf", str(_called_vcf(sample)), "--vcf-allow-no-nonvar",
                    "--var-filter",
                    "--autosome", "--set-all-var-ids", "@:#", "--rm-dup", "exclude-all",
                    "--make-pgen", "--out", str(prefix)],
                   check=True, stdout=open(config.LOGS_DIR / f"prs_pgen_{sample}.log", "wb"),
                   stderr=subprocess.STDOUT)
    # The VCF's own sample column is whatever name was baked in at call time (e.g. a
    # pre-family-migrate "father"/"mother") — not necessarily our pipeline's `sample` id.
    # Single-individual pgen, so just force the .psam IID to match; callers elsewhere
    # match reference-vs-target rows by this id.
    psam = prefix.with_suffix(".psam")
    lines = psam.read_text().splitlines()
    psam.write_text(lines[0] + "\n" + "\t".join([sample] + lines[1].split("\t")[1:]) + "\n")
    write_stamp(prefix, token)   # last: a crash mid-build leaves the cache marked stale
    return prefix


def _score_plink2(sample: str, pid: str, variants, out_dir: Path) -> tuple[float, int, int]:
    """Genome-wide score via `plink2 --score`. Returns (Σ dosage·weight, variants used,
    total). Scores over variants present in the WGS calls; absent (hom-ref) sites contribute
    0 — exact for alt-effect alleles, an undercount for ref-effect alleles pending force-call."""
    prefix = _pgen(sample)
    sf = out_dir / f"{pid}.score.txt"
    # Via _write_score_file, which drops duplicated coordinates. This path used to write the
    # file inline and kept them, so plink2 aborted the whole score with
    #   Error: --score: ALT1 allele for variant '8:141738107' appears multiple times
    # whenever the sample happened to carry a position the scoring file lists twice with the
    # allele plink2 knows as ALT1. That is sample-dependent, so the same score could fail for
    # one person and succeed for the next.
    _write_score_file(sf, variants)
    op = out_dir / f"{sample}.{pid}"
    log.info("plink2 --score %s (%d variants) for %s …", pid, len(variants), sample)
    logf = config.LOGS_DIR / f"prs_{pid}_{sample}.log"
    # no-mean-imputation: a single-sample pgen has no allele frequencies to impute from, and
    # we want absent/missing to contribute 0 (not a population mean) — our scoring convention.
    subprocess.run([_plink2_bin(), "--pfile", str(prefix), "--score", str(sf),
                    "1", "2", "3", "header", "cols=+scoresums", "no-mean-imputation",
                    "--out", str(op)],
                   check=True, stdout=open(logf, "wb"), stderr=subprocess.STDOUT)
    rows = (op.parent / f"{op.name}.sscore").read_text().splitlines()
    header = rows[0].lstrip("#").split("\t")
    cols = {c: i for i, c in enumerate(header)}
    score = float(rows[1].split("\t")[cols["SCORE1_SUM"]])
    m = re.search(r"(\d[\d,]*) variants? processed", logf.read_text())
    used = int(m.group(1).replace(",", "")) if m else 0
    return score, used, len(variants)


# --- genome-wide percentile: empirical, against the 1000G reference panel ----
# The analytic AF model (small scores) doesn't scale to millions of variants, so we score the
# SAME PGS over the 1000G individuals with plink2 and take an empirical percentile. PRS is
# additive, so we score per-chromosome and sum. The per-chromosome reference pgen (1000G
# genotypes at the score's positions) is sample-independent and cached, so the heavy network
# extraction happens once per PGS.
#
# This used to restrict the reference scoring to the sample's own plink2-covered sites, on the
# reasoning that scoring both sides over the same variant set made the comparison fair. It did
# not. The sample's VCF is **variant-only**, so its covered sites are precisely the sites where
# the sample is non-reference — choosing the comparison set from the sample's own genotypes is
# textbook ascertainment bias. The sample carried alt alleles at essentially every scored site
# by construction while reference individuals carried them at population frequency, so every
# sample landed in the top percentile.
#
# Both sides are now scored over the FULL scoring-file variant set. That alone would trade an
# upward bias for a downward one, because plink2 cannot see the sample's hom-ref sites at all:
# where a score variant's effect allele is the REFERENCE allele, absent-from-VCF means the
# sample carries TWO copies of the effect allele, not zero. `_ref_effect_supplement` adds that
# dosage back, reading REF from the reference .pvar we already build. The result is an exact
# dosage under the pipeline-wide convention that a site absent from a PASS-filtered
# variant-only VCF is homozygous reference — the same assumption `traits`/`hla` make, and the
# one `callable` (C3) exists to bound.
def _kg_population_map() -> dict[str, str]:
    """1000G sample ID → superpopulation (EUR/EAS/AMR/SAS/AFR), from the cached panel file."""
    if not config.KG_PANEL_FILE.exists():
        config.PGS_DIR.mkdir(parents=True, exist_ok=True)
        log.info("Downloading 1000G population map …")
        try:
            urllib.request.urlretrieve(config.KG_PANEL_URL, config.KG_PANEL_FILE)
        except Exception as e:  # noqa: BLE001 — best-effort, percentile degrades to none
            log.warning("1000G population map download failed (%s)", e)
            return {}
    lines = config.KG_PANEL_FILE.read_text().splitlines()
    if not lines:
        return {}
    header = lines[0].split()
    name_pairs = (("SampleID", "Superpopulation"),   # original 3202-sample ped_population file
                  ("sample", "super_pop"))            # phase-3 2504-sample .panel fallback
    for sample_col, pop_col in name_pairs:
        if sample_col in header and pop_col in header:
            si, pi = header.index(sample_col), header.index(pop_col)
            break
    else:                                  # no recognized header — documented column order
        si, pi = 0, 5
    out = {}
    for line in lines[1:]:
        f = line.split()
        if len(f) > max(si, pi):
            out[f[si]] = f[pi]
    return out


def _kg_ref_pgen(chrom: str, chrom_vars, workdir: Path) -> Path | None:
    """Cached per-chromosome reference pgen of 1000G genotypes at the score's positions.
    Sample-independent → built once per PGS and reused for the whole family."""
    prefix = workdir / f"chr{chrom}"
    if prefix.with_suffix(".pgen").exists():
        return prefix
    workdir.mkdir(parents=True, exist_ok=True)
    # Per-process temp names. The {pid}_1kg workdir is deliberately shared across the family
    # (the reference panel is sample-independent), so fixed temp names let two members scored
    # concurrently unlink and overwrite each other's intermediates mid-read — the same race
    # already fixed in _present_genotypes, missed here. The CACHED outputs below stay shared.
    tag = os.getpid()
    reg = workdir / f".reg_chr{chrom}.{tag}.tmp"
    reg.write_text("".join(f"chr{chrom}\t{p}\n" for _c, p, *_ in chrom_vars))
    subset = workdir / f".subset_chr{chrom}.{tag}.vcf.gz"
    url = config.kg_vcf_source(chrom)
    log.info("  1000G genotype extract chr%s (%d positions) …", chrom, len(chrom_vars))
    try:
        # -T, not -R. Both select by position file, but -R seeks through the tabix index per
        # region while -T streams the file once. At this scale -R is pathological: 86,160
        # scattered positions on chr22 (0.5 GB) had not finished after 13 minutes, where -T
        # took 30 seconds — and chr1 is five times the size. That single flag is most of the
        # ~38 h the serial panel build was projected to take.
        #
        # It is also the correct semantic here, not merely the fast one. -R additionally
        # returns records that merely OVERLAP a target (a deletion starting earlier and
        # spanning it); -T matches on start position only. Measured on 1,500 positions, -R
        # returned 16 such extra records and -T returned nothing -R did not. Those 16 can
        # never match a score entry anyway: the panel is built with --set-all-var-ids @:#, so
        # a record is keyed by its own start, while the score file is keyed by the marker
        # position. They would only add variants for --rm-dup to trip over.
        subprocess.run(["bcftools", "view", "-T", str(reg),
                        "--threads", str(config.PRS_BCFTOOLS_THREADS),
                        "-Oz", "-o", str(subset), url],
                       check=True, capture_output=True, timeout=7200)
    except Exception as e:  # noqa: BLE001 — best-effort network extraction
        log.warning("  1000G extract failed chr%s (%s) — genome-wide percentile skipped.",
                    chrom, e)
        reg.unlink(missing_ok=True)
        subset.unlink(missing_ok=True)
        return None
    reg.unlink(missing_ok=True)
    # --max-alleles 2: keep clean biallelic records so the score's A1 matches unambiguously.
    # plink2 strips any "chr" prefix on output regardless of source contig naming (verified
    # empirically), so @:# var IDs come out bare ("1:pos") — matched by the score file.
    # Build under a per-process prefix, then rename each part into place: a second run that
    # finds .pgen present must not find a half-written .pvar/.psam beside it.
    build = workdir / f".build_chr{chrom}.{tag}"
    # APPEND the extension, never Path.with_suffix: this prefix carries a per-process tag as
    # its own dot-segment, so with_suffix(".pgen") on ".build_chr1.68067" yields
    # ".build_chr1.pgen" — a name plink2 never writes. The existence check below therefore
    # reported "build failed" on every successful build, every genome-wide percentile came
    # back n/a, and the cleanup missed the real files too, leaking gigabytes of orphans. repeats.py
    # carries the same warning; this call site did not get the memo.
    ext_of = lambda e: workdir / f".build_chr{chrom}.{tag}{e}"
    pp = subprocess.run([_plink2_bin(), "--vcf", str(subset), "--max-alleles", "2",
                         "--set-all-var-ids", "@:#", "--rm-dup", "exclude-all",
                         "--make-pgen", "--out", str(build)],
                        capture_output=True, text=True)
    subset.unlink(missing_ok=True)
    if not ext_of(".pgen").exists():
        log.warning("  plink2 reference pgen build failed chr%s (rc=%s): %s", chrom,
                    pp.returncode, " ".join(pp.stderr.split())[:200] or "no .pgen produced")
        for e in (".pgen", ".pvar", ".psam", ".log"):
            ext_of(e).unlink(missing_ok=True)
        return None
    for ext in (".pvar", ".psam", ".pgen"):   # .pgen last: it is the existence check above
        src = ext_of(ext)
        if src.exists():
            os.replace(src, prefix.with_suffix(ext))
    ext_of(".log").unlink(missing_ok=True)
    return prefix


def _write_score_file(path: Path, cov) -> None:
    """plink2 --score file (ID/A1/W) for the covered variants. IDs must be the BARE
    'chrom:pos' — plink2 strips any 'chr' prefix on pgen output (verified
    empirically); a 'chr' prefix here matches zero variants and silently kills the percentile.

    Positions carrying more than one scoring row are DROPPED, which is what the reference
    panel already does to itself via `--rm-dup exclude-all`. plink2 keys variants as
    chrom:pos, so two rows at one coordinate are indistinguishable to it and it aborts the
    entire chromosome:

        Error: --score: ALT1 allele for variant '19:39893126' appears multiple times

    One such position is enough to void a score's genome-wide percentile — PGS003446 has
    exactly one, in 538,001 variants, and its percentile has been n/a because of it.
    """
    counts = Counter((c, p) for c, p, *_ in cov)
    dups = {k for k, n in counts.items() if n > 1}
    with open(path, "w") as fh:
        fh.write("ID\tA1\tW\n")
        for c, p, ea, _oa, w in cov:
            if (c, p) in dups:
                continue
            fh.write(f"{c}:{p}\t{ea}\t{w}\n")
    if dups:
        log.info("  %s: dropped %d row(s) at %d duplicated coordinate(s) — plink2 cannot "
                 "tell them apart by chrom:pos and the panel excludes them too.",
                 path.name, sum(counts[k] for k in dups), len(dups))


_AUTOSOMES = frozenset(str(i) for i in range(1, 23))


def _autosomal(variants):
    """Score variants on the canonical autosomes 1-22.

    The sample's pgen is built with `--autosome`, and the 1000G reference panel only exists
    per canonical chromosome — a PGS scoring file that includes ALT/patch contigs (PGS000018
    carries one on chr17_KI270857v1_alt) has no counterpart on either side. Previously these
    were filtered out only as a side effect of restricting to the sample's covered sites;
    scoring the full variant set surfaced them, and a single failed 1000G fetch aborted the
    entire percentile. Filtering explicitly keeps both sides on the same variant set for a
    reason that has nothing to do with the sample's own genotypes."""
    return [v for v in variants if v[0] in _AUTOSOMES]


def _ref_score_chrom(chrom: str, chrom_vars, workdir: Path):
    """Build (or reuse) the 1000G panel for one chromosome and score every reference
    individual over this PGS's variants there. Returns (chrom, {iid: score}) or (chrom, None).

    Split out of _reference_scores so the chromosomes can run concurrently. Safe to do so:
    every temp path it touches carries the chromosome in its name, so two workers never
    collide on one file, and it returns its totals instead of accumulating into shared state.
    """
    prefix = _kg_ref_pgen(chrom, chrom_vars, workdir)
    if prefix is None:
        return chrom, None
    sf = workdir / f".score_all_chr{chrom}.{os.getpid()}.txt"
    _write_score_file(sf, chrom_vars)
    op = workdir / f".sc_all_chr{chrom}.{os.getpid()}"
    sscore = Path(f"{op}.sscore")
    try:
        pp = subprocess.run([_plink2_bin(), "--pfile", str(prefix), "--score", str(sf),
                             "1", "2", "3", "header", "cols=+scoresums", "no-mean-imputation",
                             "--out", str(op)], capture_output=True, text=True)
        if pp.returncode:
            # Same reasoning as the bcftools guard in _reference_freqs: plink2 can write a
            # partial .sscore and then fail, and a partial reference panel does not degrade
            # the percentile, it inflates it. Treat the chromosome as missing.
            log.warning("  plink2 reference scoring chr%s rc=%s: %s", chrom, pp.returncode,
                        " ".join(pp.stderr.split())[:200])
            return chrom, None
        rows = sscore.read_text().splitlines() if sscore.exists() else []
        hdr = {c: i for i, c in enumerate(rows[0].lstrip("#").split("\t"))} if rows else {}
        ii, si = hdr.get("IID"), hdr.get("SCORE1_SUM")
        if ii is None or si is None:
            # plink2 can exit 0 having produced nothing usable; such a chromosome used to
            # just vanish from the sums.
            return chrom, None
        out = {}
        for r in rows[1:]:
            f = r.split("\t")
            out[f[ii]] = float(f[si])
        return chrom, out
    finally:
        sscore.unlink(missing_ok=True)
        Path(f"{op}.log").unlink(missing_ok=True)
        sf.unlink(missing_ok=True)


def _reference_scores(pid: str, variants) -> list[tuple[str, float]] | None:
    """Per-individual 1000G scores for this PGS over the full scoring-file variant set.

    Returns [(superpopulation, score)] or None if the reference scoring can't be completed.
    Cached per **PGS** — no longer per sample, because the variant set no longer depends on
    the sample (that dependence was the bias; see the section comment above). One consequence
    worth having: the reference distribution is now computed once and reused family-wide."""
    workdir = config.PGS_DIR / f"{pid}_1kg"
    cache = workdir / "scores_all.tsv"
    if cache.exists():
        rows = (ln.split("\t") for ln in cache.read_text().splitlines())
        return [(f[1], float(f[2])) for f in rows if len(f) == 3]
    popmap = _kg_population_map()
    if not popmap:
        return None
    # The cache write at the bottom needs this to exist. It only ever did because
    # _kg_ref_pgen happened to create it first, which is not this function's to rely on.
    workdir.mkdir(parents=True, exist_ok=True)
    by_chrom: dict[str, list] = defaultdict(list)
    for v in variants:
        by_chrom[v[0]].append(v)
    items = [(c, v) for c, v in sorted(by_chrom.items()) if v]
    workers = max(1, min(config.PRS_REF_WORKERS, len(items)))
    log.info("  1000G reference panel for %s: %d chromosome(s), %d in parallel "
             "(%d bcftools thread(s) each).", pid, len(items), workers,
             config.PRS_BCFTOOLS_THREADS)
    totals: dict[str, float] = defaultdict(float)
    scored_any = False
    incomplete: list[str] = []   # chromosomes whose plink2 --score produced nothing
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(_ref_score_chrom, c, v, workdir): c for c, v in items}
        for fut in as_completed(futures):
            chrom, got = fut.result()
            if got is None:
                incomplete.append(chrom)
            else:
                for iid, sc in got.items():       # summed here, not in the workers, so the
                    totals[iid] += sc             # accumulator needs no lock
                scored_any = True
    if not scored_any:
        return None
    if incomplete:
        # A reference distribution missing a chromosome is not a degraded percentile, it is a
        # WRONG one: the sample is still scored genome-wide, so every reference individual
        # sits artificially low and the sample's percentile is inflated. `_reference_freqs`
        # refuses to cache a partial result for the same reason; this must also refuse to
        # *use* one, and the cache is now per-PGS, so writing it would stick family-wide.
        log.warning("  1000G scoring produced no output for chr%s — genome-wide percentile "
                    "skipped and NOT cached (a partial reference panel inflates it).",
                    ", chr".join(incomplete))
        return None
    result = [(popmap[iid], s) for iid, s in totals.items() if iid in popmap]
    # Atomic: a concurrent scorer reading this cache must see the previous complete file or
    # the new complete file, never the half of it that had been flushed so far.
    ctmp = cache.with_name(f"{cache.name}.{os.getpid()}.part")
    ctmp.write_text("".join(f"{iid}\t{popmap[iid]}\t{s}\n"
                            for iid, s in totals.items() if iid in popmap))
    os.replace(ctmp, cache)
    return result


def _ref_alleles(variants, workdir: Path) -> dict[tuple[str, str], str]:
    """(chrom, pos) → REF allele for the supplement below.

    When the 1000G reference pvars exist they are the ONLY source used, and positions they
    do not carry are deliberately left unresolved: the reference panel was not scored on
    those variants either, so crediting the sample for them would make the two sides cover
    different variant sets and inflate the percentile — the asymmetry this whole supplement
    exists to remove.

    When no pvars exist at all the reference run never happened, so there is no percentile to
    stay symmetric with, and REF comes from the local FASTA instead. That is what keeps the
    reported raw score a property of the genome: previously a missing reference panel left
    the supplement silently at zero, so the same sample scored differently depending on
    whether the network was up when it ran.
    """
    want = {(v[0], v[1]) for v in variants}
    pvars = sorted(workdir.glob("chr*.pvar"))
    if pvars:
        refs: dict[tuple[str, str], str] = {}
        for pvar in pvars:
            for line in pvar.read_text().splitlines():
                if line.startswith("#"):
                    continue
                f = line.split("\t")
                if len(f) >= 4 and (f[0], f[1]) in want:
                    refs[(f[0], f[1])] = f[3].upper()      # #CHROM POS ID REF ALT
        return refs
    if not want:
        return {}
    # Bounded at the same threshold that separates a small score from a genome-wide one.
    # Below it the FASTA pass is a handful of faidx batches and buys a raw score that does
    # not move with network luck. Above it, it is millions of lookups to adjust a number
    # that has no reference distribution to be comparable with — the supplement exists only
    # for that comparison. Measured unbounded: 29.3M lookups per sample across 22 scores,
    # 3.5 h for one sample, and a correction larger than the score itself applied to a score
    # reported with no percentile at all. The caller reports the plain plink2 sum instead.
    if len(want) > config.PRS_PLINK_MIN:
        log.info("  no cached 1000G pvars under %s and %d positions to resolve — skipping the "
                 "hom-ref supplement; the raw score is the plink2 sum over the sample's own "
                 "called sites.", workdir.name, len(want))
        return {}
    log.info("  no cached 1000G pvars under %s — resolving %d reference allele(s) from %s.",
             workdir.name, len(want), config.REF_FASTA.name)
    return {k: b.upper() for k, b in faidx_bases(sorted(want)).items() if b}


def _ref_effect_supplement(variants, covered, workdir: Path) -> tuple[float, int]:
    """(dosage owed, n variants) at score sites absent from the sample's variant-only VCF.

    A variant-only VCF lists only sites where the sample is non-reference, so a score variant
    missing from it means homozygous reference. Where the effect allele is the ALT that
    correctly contributes zero — but where the **effect allele is the REFERENCE allele**,
    hom-ref means the sample carries two copies of it, and plink2 scores it as zero because
    the site simply is not in the pgen. Without this the sample's score is biased *downward*
    by 2·w for every such variant, which is exactly the error introduced by scoring the
    reference over the full variant set and the sample over only its own sites.

    Variants whose REF we cannot establish (absent from the 1000G reference pvar) are skipped:
    they are equally absent from the reference scoring, so the comparison stays symmetric."""
    # Resolve REF only where it is actually needed — a covered site already has an explicit
    # genotype in the pgen. On the FASTA path above this is the difference between a handful
    # of lookups and one per score variant.
    uncovered = [v for v in variants if (v[0], v[1]) not in covered]
    refs = _ref_alleles(uncovered, workdir)
    total, n = 0.0, 0
    for chrom, pos, ea, _oa, w in uncovered:
        ref = refs.get((chrom, pos))
        if ref is not None and ea.upper() == ref:
            total += 2.0 * w
            n += 1
    return total, n


def _sample_covered(prefix: Path) -> set[tuple[str, str]]:
    """(chrom,pos) the sample has an explicit genotype for, from its plink2 .pvar."""
    covered: set[tuple[str, str]] = set()
    for line in prefix.with_suffix(".pvar").read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) >= 2:
            covered.add((f[0], f[1]))
    return covered


def _empirical_percentile(score: float, ref: list[tuple[str, float]]):
    """(EUR percentile, overall percentile) — share of reference individuals scoring below."""
    def pct(vals):
        return 100.0 * sum(1 for v in vals if v < score) / len(vals) if vals else None
    eur = [s for p, s in ref if p == config.KG_POP]
    return pct(eur), pct([s for _p, s in ref])


def _score_entry(sample: str, entry: dict, gt_vcfs, out_dir: Path) -> dict:
    """Score one panel row. Raises on failure; the caller isolates it."""
    pid, trait = entry["pgs_id"], entry["trait"]
    variants = _load_variants(_scoring_file(pid))
    if not variants:
        # A score nothing could be parsed from must fail loudly. Scored as-is it yields a
        # raw 0.0 over 0 variants, which renders as an ordinary row and reads as a real
        # result — the worst possible outcome for a health report.
        raise ValueError(f"no usable variants parsed from the {pid} scoring file")
    if len(variants) > config.PRS_PLINK_MIN:
        # Genome-wide: score with plink2, then take an *empirical* percentile by scoring the
        # same PGS over the 1000G reference panel (restricted to the sample's covered
        # variants). Best-effort — the raw score still stands if the reference run can't run.
        score, used, total = _score_plink2(sample, pid, variants, out_dir)
        pe = pa = None
        n_supp = 0
        scored = _autosomal(variants)
        ref = None
        if config.PRS_GENOMEWIDE_PERCENTILE:
            try:
                ref = _reference_scores(pid, scored)
            except Exception as e:  # noqa: BLE001 — percentile is best-effort
                log.warning("  genome-wide reference panel unavailable for %s (%s)", pid, e)
        # The hom-ref supplement is a property of the SAMPLE, not of the reference run: a
        # variant-only VCF omits sites where the sample is hom-ref, and where the effect
        # allele IS the reference the sample carries two copies plink2 cannot see. Adding it
        # only when the reference scoring happened to succeed made the *reported raw score*
        # depend on network availability — the same genome, a different headline number on a
        # bad day. It is computed after the reference run so the cached pvars are used when
        # they exist (see _ref_alleles on why that matters for symmetry).
        try:
            covered = _sample_covered(_pgen(sample))
            supp, n_supp = _ref_effect_supplement(
                scored, covered, config.PGS_DIR / f"{pid}_1kg")
            score += supp
            if n_supp:
                log.info("  %s: +%.4f from %d reference-effect allele(s) the "
                         "variant-only VCF omits (hom-ref = 2 copies).", pid, supp, n_supp)
        except Exception as e:  # noqa: BLE001 — the plink2 sum still stands
            log.warning("  %s: hom-ref supplement could not be computed (%s) — the raw score "
                        "covers only the sample's own called sites.", pid, e)
        if ref:
            try:
                pe, pa = _empirical_percentile(score, ref)
            except Exception as e:  # noqa: BLE001 — percentile is best-effort
                log.warning("  genome-wide percentile unavailable for %s (%s)", pid, e)
        log.info("  %-32s %s  plink2 raw=%.4f  EUR %%ile=%s (%d/%d variants used)",
                 trait, pid, score, f"{pe:.1f}" if pe is not None else "n/a",
                 used + n_supp, total)
        return {**entry, "trait": trait, "pid": pid, "method": "plink2", "raw": score,
                "used": used + n_supp, "total": total, "pe": pe, "pa": pa}
    else:
        score, n_used, _n_carry, contrib = _score_one(gt_vcfs, variants)
        try:
            # _autosomal, matching the plink2 path above: a scoring file with an
            # alt/patch-contig variant (PGS000018 carries one on chr17_KI270857v1_alt)
            # has no 1000G counterpart, so that chromosome always returns zero rows, the
            # partial-result guard correctly refuses to cache — and every run then
            # re-queried every chromosome over the network, forever.
            freqs = _reference_freqs(pid, _autosomal(variants))
        except Exception as e:  # noqa: BLE001 — calibration is best-effort
            log.warning("  calibration unavailable for %s (%s)", pid, e)
            freqs = {}
        pe, pa, _nc = _calibrate(contrib, variants, freqs) if freqs else (None, None, 0)
        log.info("  %-32s %s  raw=%.4f  EUR %%ile=%s (%d/%d variants)", trait, pid, score,
                 f"{pe:.1f}" if pe is not None else "n/a", n_used, len(variants))
        return {**entry, "trait": trait, "pid": pid, "method": "exact", "raw": score,
                "used": n_used, "total": len(variants), "pe": pe, "pa": pa}


def _panel_scores() -> list[dict]:
    """Panel rows to score, newest-format panel preferred, PRS_SCORES as the fallback.

    Each row carries the metadata the report needs alongside the score itself — trait,
    actionability, licence, published discrimination — so a reader can see *why* a trait is
    in the report and how much weight the score deserves, rather than a bare percentile.
    """
    path = config.COMPLEX_PANEL
    if not path.exists():
        log.info("No trait panel at %s — falling back to config.PRS_SCORES.", path)
        return [{"pgs_id": pid, "trait": trait, "actionability": "review",
                 "s1_gate": "no", "licence": "", "metric": "", "metric_kind": "",
                 "n_variants": ""} for pid, trait in config.PRS_SCORES]
    rows = []
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not row.get("pgs_id"):
                continue
            rows.append(row)
    if not rows:
        raise SystemExit(f"{path} has no usable rows.")
    return rows


def _panel_filters(rows: list[dict], sample: str) -> tuple[list[dict], list[str]]:
    """Apply the S1 consent gate and the variant-count tier cap. Returns (kept, notes).

    The gate is per PERSON. It used to consult only config.PRS_S1_GATE, a single global flag,
    so one setting decided disclosure for the whole family: turning it off to see your own
    severe non-actionable traits released everyone else's in the same run, and turning it on
    withheld them from people who had consented. consent.tsv already records this per
    (sample, category) and repeats.py already honours it; PRS did not consult it at all.

    The global flag is kept as the policy switch. It is on by default, and a trait is then
    withheld only from a sample that has not consented to that category; with it off
    (PRS_S1_GATE=0 / --no-s1-gate) nothing is gated for anyone.
    """
    notes, kept = [], []
    consented = consent.has_consented(sample, consent.CATEGORY_ADULT_ONSET_UNTREATABLE)
    gate_on = config.PRS_S1_GATE and not consented
    gated = [r for r in rows if gate_on and r.get("s1_gate") == "yes"]
    for row in rows:
        if row in gated:
            continue
        cap = config.PRS_MAX_VARIANTS
        n = int(row["n_variants"]) if str(row.get("n_variants") or "").isdigit() else 0
        if cap is not None and n > cap:
            continue
        kept.append(row)
    if gated:
        notes.append(f"**{len(gated)} trait(s) withheld** by the S1 incidental-findings "
                     f"consent filter (severe and non-actionable). {sample} has not "
                     f"consented to the `{consent.CATEGORY_ADULT_ONSET_UNTREATABLE}` "
                     "category — record consent for this person "
                     "(`run.py consent ...`), or set `PRS_S1_GATE=0` to disable the "
                     "gate for everyone.")
    skipped = len(rows) - len(gated) - len(kept)
    if skipped:
        notes.append(f"**{skipped} score(s) skipped** — over the "
                     f"`PRS_MAX_VARIANTS = {config.PRS_MAX_VARIANTS:,}` tier cap. Re-run "
                     "without the cap to score them (scoring files and reference panels "
                     "cache, and are reused family-wide).")
    return kept, notes


def _band(pct: float | None) -> str:
    if pct is None:
        return "—"
    if pct >= 90:
        return "high (top 10%)"
    if pct >= 75:
        return "above average"
    if pct >= 25:
        return "average"
    if pct >= 10:
        return "below average"
    return "low (bottom 10%)"


# Report prominence follows actionability, not score strength: a 60th-percentile result for
# something with a screening pathway is worth more of the reader's attention than a 95th for
# something nothing can be done about. Ordering the table by percentile alone would invert that.
_ACTION_SECTIONS = [
    ("treat", "Preventive treatment exists"),
    ("screen", "Screening or surveillance pathway exists"),
    ("life", "Lifestyle-modifiable"),
    ("none", "No established intervention"),
    ("review", "Actionability not yet reviewed"),
]


def _pct(r: dict) -> str:
    if r["pe"] is not None:
        return f"{r['pe']:.1f}%"
    return "— (genome-wide)" if r["method"] == "plink2" else "—"


def _result_tables(rows: list[dict]) -> list[str]:
    """Result tables grouped by actionability, each sorted by percentile (highest first)."""
    def pct_of(r):
        return r["pe"] if r["pe"] is not None else -1.0

    out: list[str] = []
    for key, heading in _ACTION_SECTIONS:
        group = sorted((r for r in rows if (r.get("actionability") or "review") == key),
                       key=pct_of, reverse=True)
        if not group:
            continue
        out += [
            f"### {heading}",
            "",
            "| Trait | PGS ID | Raw score | EUR percentile | Band (EUR) | "
            "Variants (used / total) | Score AUROC/C | Engine | Licence |",
            "|-------|--------|----------:|---------------:|------------|"
            "-------------------------|--------------:|--------|---------|",
        ]
        for r in group:
            metric = f"{r['metric']} ({r['metric_kind']})" if r.get("metric") else "—"
            out.append(
                f"| {r['trait']} | {r['pid']} | {r['raw']:.4f} | {_pct(r)} | {_band(r['pe']) if r['pe'] is not None else '—'} | "
                f"{r['used']} / {r['total']} | {metric} | {r['method']} | "
                f"{r.get('licence') or '—'} |")
        out.append("")
    return out


def _failure_lines(failures: list[tuple[str, str, str]]) -> list[str]:
    """Name what did not score. A silently short table reads as a clean run."""
    if not failures:
        return []
    return ["### Scores that failed", "",
            f"{len(failures)} score(s) could not be computed this run; every other trait "
            "above is unaffected. Re-run to retry — scoring files and reference panels cache.",
            "",
            "| Trait | PGS ID | Error |", "|-------|--------|-------|",
            *[f"| {t} | {p} | {e[:120]} |" for t, p, e in failures], ""]


def s1_token(sample: str) -> str:
    """The disclosure state a PRS report for `sample` is rendered under right now — the
    same test `_panel_filters` applies."""
    withheld = config.PRS_S1_GATE and not consent.has_consented(
        sample, consent.CATEGORY_ADULT_ONSET_UNTREATABLE)
    return f"s1={'withheld' if withheld else 'shown'}"


def consent_stamp(sample: str) -> Path:
    """Records which `s1_token` the report on disk was rendered under (cf. repeats.py)."""
    return config.PRS_DIR / sample / f"prs_{sample}.consent-stamp"


def _consent_current(sample: str) -> bool:
    """Was the report on disk rendered under today's disclosure state?

    The report's content depends on it, so existence alone is not enough: the cached report
    was returned before `_panel_filters` ever ran, and a report scored while consent stood
    kept its gated traits after consent was withdrawn. A report with no stamp predates
    stamping; it cannot disclose more than "shown" allows, so it is reused in that state
    and regenerated in the other."""
    token = s1_token(sample)
    try:
        return consent_stamp(sample).read_text().strip() == token
    except OSError:
        return token == "s1=shown"


def panel_token() -> str:
    """What a report was scored against: the trait panel's content, and the variant cap.

    A digest of every panel row, so registering a score — or correcting a trait name, an
    actionability class or a published AUROC, all of which are printed — makes the reports
    scored before it stale. The cap is part of it because a capped run is a partial report:
    it must not pass for a full one when the next uncapped run comes round."""
    rows = sorted(({str(k): v for k, v in r.items()} for r in _panel_scores()),
                  key=lambda r: r.get("pgs_id") or "")
    digest = hashlib.sha256(json.dumps(rows, sort_keys=True).encode()).hexdigest()[:16]
    cap = config.PRS_MAX_VARIANTS
    return f"panel=sha256:{digest};cap={'none' if cap is None else cap}"


def panel_stamp(sample: str) -> Path:
    """Records the `panel_token` the report on disk was scored against."""
    return config.PRS_DIR / sample / f"prs_{sample}.panel-stamp"


def _panel_current(sample: str) -> bool:
    """A report with no stamp predates stamping, so what it was scored against is unknown —
    and unknown is not current."""
    try:
        return panel_stamp(sample).read_text().strip() == panel_token()
    except OSError:
        return False


def _report_is_current(sample: str) -> bool:
    """May the existing report be reused? Only while both things it depends on still hold:
    the disclosure state it was rendered under and the panel it was scored against."""
    return _consent_current(sample) and _panel_current(sample)


def stale_reports(samples: list[str]) -> list[str]:
    """Callsets that HAVE a PRS report which no longer matches the panel or consent state.

    What `check-releases` rescoring acts on. Deliberately not "every callset without a
    current report": a callset that was never scored stays unscored until someone runs
    `prs` for it, so the scheduled job never starts the reference-panel build (hours, tens
    of GB) on a machine where nobody asked for polygenic scores."""
    return [s for s in samples
            if reportpaths.genome_report("prs", s).exists() and not _report_is_current(s)]


def prs_report(sample: str, force: bool = False) -> Path:
    require_tools("bcftools", "samtools")
    vcf = _called_vcf(sample)
    if not vcf.exists():
        raise SystemExit(f"{vcf} missing. Run `call {sample}` first.")
    out_dir = config.PRS_DIR / sample
    summary = reportpaths.genome_report("prs", sample)
    if summary.exists() and not force:
        if _report_is_current(sample):
            log.info("%s PRS already reported (%s) — skipping (use --force).",
                     sample, summary)
            return summary
        log.info("%s PRS: %s since that report was written — regenerating it.", sample,
                 " and ".join(
                     why for why, ok in (
                         (f"the S1 consent state changed (now {s1_token(sample)})",
                          _consent_current(sample)),
                         ("the trait panel changed", _panel_current(sample))) if not ok))
    out_dir.mkdir(parents=True, exist_ok=True)

    # Small (per-variant) scores prefer the force-called VCF: an evidenced 0/0 replaces an
    # assumed one, and a no-call (./.) is dropped from the score rather than silently
    # counted as hom-ref. It is a *preference*, not a substitution — the called VCF stays
    # behind it in the chain. Substituting it meant every marker outside the curated
    # force-call union (a score added since the sites were last built, or one whose scoring
    # file failed to download at build time) fell straight through to "assume reference",
    # discarding the real genotype the called VCF holds. Adding force-called data could
    # therefore make a score *less* accurate than not having it.
    # Genome-wide scores stay on the plink2 path regardless.
    forced = forcecall.available(sample)
    gt_vcfs = [forcecall.forcecalled_vcf(sample), vcf] if forced else [vcf]
    if forced:
        log.info("PRS: small scores prefer the force-called VCF, falling back to %s.",
                 vcf.name)

    panel, notes = _panel_filters(_panel_scores(), sample)
    log.info("PRS: %d trait(s) to score.", len(panel))

    rows, failures = [], []
    for entry in panel:
        pid, trait = entry["pgs_id"], entry["trait"]
        # One score's failure must not cost the rest of the run. At three hardcoded scores a
        # bare exception was survivable; over a 44-trait panel a single 404 or truncated
        # download would discard every score already computed, so each is isolated and the
        # failures are reported together at the end.
        try:
            rows.append(_score_entry(sample, entry, gt_vcfs, out_dir))
        except Exception as e:  # noqa: BLE001 — one bad score must not sink the panel
            log.warning("  %-32s %s  FAILED (%s)", trait, pid, e)
            failures.append((trait, pid, str(e)))
    if failures and not rows:
        raise SystemExit(f"every score in the panel failed ({len(failures)}); "
                         "see the warnings above.")
    return write_prs(sample, rows, notes, failures, forced)


def write_prs(sample: str, rows: list[dict], notes: list[str],
              failures: list[tuple[str, str, str]], forced: bool) -> Path:
    """Write the PRS report from scored rows (the shape `_score_entry` returns)."""
    out_dir = config.PRS_DIR / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = reportpaths.genome_report("prs", sample)

    lines = [
        f"# Polygenic risk scores — {sample}",
        "",
        f"- Small scores are calibrated to a **{config.KG_POP}** 1000 Genomes percentile "
        "(HWE + CLT analytic model); genome-wide scores are computed with **plink2 --score** "
        "and placed on an **empirical** 1000 Genomes percentile (the same PGS scored over the "
        "reference panel).",
        "",
        *_result_tables(rows),
        "",
        *([f"> {n}" for n in notes] + [""] if notes else []),
        *_failure_lines(failures),
        "## How to read this",
        "- Traits are grouped by **actionability**, not by score: a moderate percentile for "
        "something with a screening pathway matters more than a high one for something with "
        "no established intervention.",
        "- **Score AUROC/C** is the published discrimination of the *score itself*, evaluated "
        "PGS-alone (no age/sex/covariates) in an ancestry-matched cohort — it bounds how much "
        "any percentile here can mean. A 0.55 score separates cases from controls barely "
        "better than a coin flip, whatever percentile it returns.",
        "- The **EUR percentile** (small scores) is where this person's raw score falls in the "
        "European 1000 Genomes distribution — e.g. 90% means a higher genetic load than 90% of "
        "that reference. A percentile is *relative genetic load*, not absolute risk.",
        "- **Genome-wide scores** (engine `plink2`) report the raw weighted dosage over the "
        "~millions of score variants. Their **EUR percentile** is *empirical*: the same PGS is "
        "scored over the 1000 Genomes individuals and the sample is ranked among them. Both "
        "sides are scored over the **full scoring-file variant set** (autosomes 1-22). "
        "Restricting to the sites the sample has an explicit genotype for would be an "
        "ascertainment bias that pushes every sample toward the top percentile: those sites are the sample's variant calls by construction, while the "
        "reference individuals carry them only at population frequency. Where a score "
        "variant's effect allele is the REFERENCE allele, absent-from-VCF means two copies, "
        "not zero; `_ref_effect_supplement` adds that dosage back.",
        "",
        "## Caveats",
        ("- **Genotype source** (small scores): read from the **force-called** VCF — evidenced "
         "hom-ref, with no-call sites dropped from the score." if forced else
         "- **Genotype source** (small scores): absent positions are taken as hom-ref from the "
         "reference base. Run `force-call` first to make this rigorous."),
        "- **Calibration model** (small scores): the population score is approximated as "
        "Normal(Σ 2·f·w, Σ 2·f(1−f)·w²) from 1000G effect-allele frequencies (HWE + CLT).",
        "- **Ancestry:** percentiles are against the EUR reference; most PGS are European-derived and "
        "miscalibrate across ancestries.",
        "- **plink2 scores** sum over the variants present in the WGS calls. A variant-only "
        "VCF has no record at hom-ref sites, so plink2 scores them as zero dosage: correct "
        "where the effect allele is the ALT, wrong where the effect allele is the REFERENCE "
        "allele (hom-ref means two copies of it). That missing dosage is added back "
        "explicitly, reading REF from the 1000G reference panel, so the sample's score covers "
        "the same variant set the reference individuals were scored over. The used / total "
        "column counts both.",
        "- **The empirical percentile scores the reference panel over the full scoring-file "
        "variant set**, not over the sample's own covered sites. Restricting it to the "
        "sample's sites — which, from a variant-only VCF, are the sites where the sample is "
        "non-reference — selected the comparison set using the sample's own genotypes and "
        "pushed every sample toward the top percentile.",
        "- **Licence** is the score's own terms of use, recorded for provenance rather than "
        "enforced: personal, non-commercial analysis satisfies every licence in the panel, and "
        "scoring files are fetched at runtime and never redistributed. It matters only if this "
        "pipeline is redeployed commercially, where non-default terms need clearing first.",
        "- **Trait selection** is generated from PGS Catalog by the criteria in "
        "`docs/polygenic-panel-design.md`, not curated by hand. Both the EFO term and the score ID in "
        "`panels/complex_traits.tsv` are human-reviewed.",
        "- Research-grade, not diagnostic.",
        "",
    ]
    write_report(summary, lines)
    consent_stamp(sample).write_text(s1_token(sample) + "\n")
    panel_stamp(sample).write_text(panel_token() + "\n")
    log.info("PRS report: %s", summary)
    return summary
