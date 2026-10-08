"""Stage: repeat-expansion genotyping via ExpansionHunter.

Short-tandem-repeat (STR) expansions cause major diseases — Huntington (HTT), Fragile X
(FMR1), ALS/FTD (C9orf72), myotonic dystrophy (DMPK), the spinocerebellar ataxias — that
SNV/indel calling is blind to (the panel report flags exactly these as "not assessable").
ExpansionHunter sizes the repeat at each catalog locus straight from the CRAM; this stage
runs it and interprets the sizes against established pathogenic thresholds.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from . import reportpaths
from . import config, consent
from .call import aligned_cram
from .util import log, require_tools, run, write_report

# Every locus here is adult-onset with no disease-modifying treatment — a pure prediction,
# not an actionable finding — so S1 gates it behind explicit consent (pipeline/consent.py).
_INCIDENTAL_CATEGORY = consent.CATEGORY_ADULT_ONSET_UNTREATABLE

# VARID → (disorder, premutation_min, pathogenic_min). Research-grade thresholds for the
# well-established loci; loci absent here are reported with their size but no call.
_DISORDERS = {
    "HTT": ("Huntington disease", 27, 36),
    "FMR1": ("Fragile X syndrome", 55, 200),
    "C9orf72": ("ALS / frontotemporal dementia", 24, 30),
    "DMPK": ("Myotonic dystrophy type 1", 35, 50),
    "FXN": ("Friedreich ataxia (recessive — needs both alleles)", None, 66),
    "AR": ("Spinal-bulbar muscular atrophy (Kennedy)", None, 38),
    "ATXN1": ("Spinocerebellar ataxia 1", None, 39),
    "ATXN2": ("Spinocerebellar ataxia 2", None, 33),
    "ATXN3": ("Spinocerebellar ataxia 3 (Machado-Joseph)", None, 56),
    "ATXN7": ("Spinocerebellar ataxia 7", None, 37),
    "CACNA1A": ("Spinocerebellar ataxia 6", None, 20),
    "TBP": ("Spinocerebellar ataxia 17", None, 49),
    "ATN1": ("Dentatorubral-pallidoluysian atrophy (DRPLA)", None, 48),
}


# ExpansionHunter's GRCh38 catalog does not spell every VARID like the HGNC gene symbol:
# the catalog says "C9ORF72" where this table says "C9orf72". An exact-match lookup left
# that locus — the commonest genetic cause of ALS/FTD — genotyped but never interpreted:
# status "—", never flagged, and (once the S1 gate keyed off the same table) never withheld
# either. Match case-insensitively so no locus can silently fall through.
_BY_ID = {k.casefold(): v for k, v in _DISORDERS.items()}


def _disorder(varid: str):
    """(name, premutation_min, pathogenic_min) for a catalog VARID, or None."""
    return _BY_ID.get((varid or "").casefold())


def is_incidental(varid: str) -> bool:
    """Is this locus in the adult-onset/untreatable category the S1 gate covers?"""
    return (varid or "").casefold() in _BY_ID


def _eh_binary() -> Path | None:
    hits = list(config.EXPANSIONHUNTER_DIR.glob("*/bin/ExpansionHunter"))
    return hits[0] if hits else None


def _eh_catalog() -> Path | None:
    hits = list(config.EXPANSIONHUNTER_DIR.glob(
        "*/variant_catalog/grch38/variant_catalog.json"))
    return hits[0] if hits else None


def _info(field: str, info: str) -> str:
    for kv in info.split(";"):
        if kv.startswith(field + "="):
            return kv.split("=", 1)[1]
    return ""


def _status(varid: str, sizes: list[int]) -> tuple[str, bool]:
    d = _disorder(varid)
    if not d or not sizes:
        return "—", False
    name, premut, path = d
    top = max(sizes)
    if top >= path:
        return f"⚠️ EXPANDED — {name} (pathogenic range, ≥{path})", True
    if premut and top >= premut:
        return f"premutation — {name} ({premut}–{path - 1})", True
    return f"normal ({name})", False


def _parse(vcf: Path) -> list[tuple]:
    """(varid, repeat_unit, sizes[list], status, flagged) per locus."""
    rows = []
    for line in vcf.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 10:   # truncated/malformed ExpansionHunter line — skip, don't crash
            continue
        info, fmt, sample = f[7], f[8].split(":"), f[9].split(":")
        varid = _info("VARID", info) or _info("REPID", info)
        unit = _info("RU", info)
        try:
            repcn = sample[fmt.index("REPCN")]
        except (ValueError, IndexError):
            continue
        sizes = [int(x) for x in repcn.replace("|", "/").split("/") if x.isdigit()]
        status, flagged = _status(varid, sizes)
        rows.append((varid, unit, sizes, status, flagged))
    return rows


def _apply_consent(rows: list, consented: bool) -> tuple[list, list]:
    """(visible, withheld) under the S1 consent gate — withholds EVERY locus in the
    incidental category (`_DISORDERS`), not only the flagged ones.

    Withholding just the flagged rows leaks the finding by elimination. Both the
    ExpansionHunter catalog and the threshold table are fixed and public (the latter is
    right here in this file), so a reader diffs the visible table against them and names
    the withheld locus immediately — and a non-zero withheld *count* already discloses
    "you carry an expansion somewhere", which is the incidental finding itself. Redacting
    the whole category is what makes an unconsented report look identical whether or not
    anything was flagged. Loci outside `_DISORDERS` carry no interpretation, so they stay
    visible."""
    if consented:
        return rows, []
    withheld = [r for r in rows if is_incidental(r[0])]
    return [r for r in rows if not is_incidental(r[0])], withheld


def repeats_report(sample: str, force: bool = False) -> Path:
    require_tools("samtools")
    eh, catalog = _eh_binary(), _eh_catalog()
    if not (eh and catalog):
        raise SystemExit(f"ExpansionHunter not installed under {config.EXPANSIONHUNTER_DIR}. "
                         "Run `setup`.")
    cram = aligned_cram(sample)
    if not cram.exists():
        raise SystemExit(f"{cram} missing. Run `align {sample}` first.")

    out_dir = config.REPEATS_DIR / sample
    summary = reportpaths.genome_report("repeats", sample)
    consented = consent.has_consented(sample, _INCIDENTAL_CATEGORY)
    # The report's CONTENT depends on the consent state, so the skip must too. Reusing it on
    # existence alone meant a consent change — in either direction — left the old disclosure
    # (or the old withholding) in place until someone remembered --force. Stamped beside the
    # report rather than via util.stamp_path, which would rewrite a dotted name's last
    # segment.
    stamp = out_dir / f"repeats_{sample}.consent-stamp"
    token = f"consent={'yes' if consented else 'no'}"
    if summary.exists() and not force:
        try:
            current = stamp.read_text().strip() == token
        except OSError:
            current = False        # never stamped (pre-stamp report) → regenerate once
        if current:
            log.info("%s repeats already reported (%s) — skipping (use --force).",
                     sample, summary)
            return summary
        log.info("%s repeats: consent for %s is now %s — regenerating the report.",
                 sample, _INCIDENTAL_CATEGORY, "granted" if consented else "withdrawn")
    out_dir.mkdir(parents=True, exist_ok=True)

    prefix = out_dir / sample
    # ExpansionHunter writes '<prefix>.vcf' LITERALLY — Path.with_suffix would replace a
    # dotted sample id's suffix instead of appending.
    vcf = prefix.parent / f"{prefix.name}.vcf"
    # Reuse existing genotypes unless forced, like the sibling CRAM stages (`sv` reuses its
    # diploidSV.vcf.gz, `hla-type` its genotype.json). Re-rendering after a REPORT-level fix
    # — a corrected threshold key, a consent-policy change — must not cost another multi-hour
    # pass over a 50GB CRAM: the repeat sizes are already computed and unaffected by it.
    if force or not vcf.exists():
        log.info("Running ExpansionHunter on %s …", sample)
        run([str(eh), "--reads", str(cram), "--reference", str(config.REF_FASTA),
             "--variant-catalog", str(catalog), "--output-prefix", str(prefix)],
            stdout=open(config.LOGS_DIR / f"repeats_{sample}.log", "wb"),
            stderr=subprocess.STDOUT)
    else:
        log.info("Reusing existing ExpansionHunter genotypes (%s).", vcf.name)
    return report_from_genotypes(sample, vcf)


def report_from_genotypes(sample: str, vcf: Path) -> Path:
    """Interpret an ExpansionHunter VCF and write the report, under the consent gate.

    Separate from `repeats_report` so the interpretation can run over genotypes that were
    produced elsewhere — the demo feeds it bundled sample output, since a synthetic family
    has no reads to size repeats from."""
    out_dir = config.REPEATS_DIR / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = reportpaths.genome_report("repeats", sample)
    stamp = out_dir / f"repeats_{sample}.consent-stamp"
    consented = consent.has_consented(sample, _INCIDENTAL_CATEGORY)
    token = f"consent={'yes' if consented else 'no'}"

    rows = _parse(vcf)
    rows.sort(key=lambda r: (not r[4], r[0]))  # flagged first, then by locus
    flagged = [r for r in rows if r[4]]
    visible, withheld = _apply_consent(rows, consented)

    def _row(r):
        varid, unit, sizes, status, _ = r
        return f"| {varid} | {unit} | {'/'.join(map(str, sizes)) or '?'} | {status} |"

    # The flagged COUNT is itself the incidental finding ("you carry an expansion somewhere"),
    # so it is stated only once consent is granted. The unconsented line is deliberately
    # constant — it depends on the size of the category, never on this person's results.
    headline = (f"- **Flagged (premutation / expanded):** {len(flagged)}" if consented else
                f"- **{len(withheld)} disease-associated repeat loci are not shown** "
                "(adult-onset / untreatable category — S1 consent gate, see below)")
    lines = [
        f"# Repeat expansions — {sample}",
        "",
        f"- **Loci genotyped:** {len(rows)} (ExpansionHunter {config.EXPANSIONHUNTER_VERSION})",
        headline,
        "",
        "## Disease-associated repeat loci",
        "",
        "| Locus | Repeat | Sizes (repeat units) | Status |",
        "|-------|--------|----------------------|--------|",
        *[_row(r) for r in visible],
    ]
    if withheld:
        lines += [
            "",
            "## Incidental findings — withheld pending consent",
            f"_The {len(withheld)} repeat loci with an established disease interpretation "
            "(Huntington, Fragile X, C9orf72/ALS, the ataxias, …) are adult-onset with no "
            "disease-modifying treatment, so under the S1 policy their sizes are withheld "
            "until you ask for them. **The whole category is withheld, including normal "
            "results** — showing the normal ones would identify any abnormal one by "
            "elimination, and even a count of abnormal loci would disclose that there is "
            "one. So this section reads exactly the same whether or not anything was "
            "flagged. This is a disclosure choice, not a data gate — every size is already "
            "computed and stored. "
            f"To view: `run.py consent {sample} {_INCIDENTAL_CATEGORY} yes`, then re-run "
            f"`run.py repeats {sample} --force`._",
        ]
    lines += [
        "",
        "## Caveats",
        "- Repeat sizing from short reads is approximate, especially for very large "
        "expansions (e.g. FMR1 full mutations, DMPK) — confirm any flag with a clinical "
        "repeat-sizing assay (PCR / Southern blot).",
        "- Thresholds are research-grade; some loci (FXN, recessive) need both alleles "
        "expanded to cause disease.",
        "- Research-grade, not a clinical diagnosis.",
        "",
    ]
    write_report(summary, lines)
    stamp.write_text(token + "\n")   # what consent state this report was rendered under
    log.info("Repeats report: %s (%d loci, %d flagged)", summary, len(rows), len(flagged))
    return summary
