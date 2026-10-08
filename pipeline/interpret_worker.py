"""Variant-effect worker — the GPU kernel behind the Mac's `interpret` stage.

Supports two interchangeable sequence-to-function models behind one interface (--model):
  • enformer — DeepMind Enformer (196,608 bp context → 896 bins × 5,313 human tracks)
  • borzoi   — Calico Borzoi, Enformer's RNA-seq successor (524,288 bp → 6,144 bins × 7,611)

Deliberately *reference-free and genome-agnostic*: it never sees the reference FASTA, a contig
name, or a sample. The Mac (which owns the genome data) extracts the ref and alt context windows
per variant — at the selected model's context length — and ships them as a compact gzipped job;
this script one-hot-encodes each pair, runs the model on the GPU, and returns the change in
predicted tracks. That keeps all genome-specific logic on the Mac and makes the box a clean,
stateless compute node (per docs/hardware-and-scope.md).

Job  (TSV, gzip, --job):   id <TAB> ref_seq <TAB> alt_seq   (seqs are <seq_len> ACGTN chars)
Out  (TSV, --out):         id  delta_max  l2_center  top_track  top_delta  top_desc

Effect summary per variant: alt minus ref over the human head's central bins (the bins flanking
the edit; the receptive-field edges are unreliable), normalised to (batch, bins, tracks) so both
models share one code path. We report the single most-changed (track, bin) magnitude `delta_max`,
the L2 of the whole track vector at the variant's own centre bin `l2_center`, and the top track's
index/description — mirroring how the AlphaMissense/SpliceAI lookups expose one magnitude + label.
"""
import argparse
import csv
import gzip
import sys
import time

import torch
import torch.nn.functional as F

# Each job field is a full context sequence (up to 524,288 chars), past csv's 128 KB default cap.
csv.field_size_limit(10 ** 8)

# Per-model specs. `channels_first` is the one-hot layout the model wants: Enformer takes
# (B, L, 4), Borzoi (a conv1d port) takes (B, 4, L) and emits tracks-first (B, tracks, bins).
MODELS = {
    # `stranded`: whether the human head's tracks come in +/- strand pairs. Enformer's do not
    # (its CAGE targets are strand-agnostic), Borzoi's RNA-seq/CAGE targets do — which decides
    # whether reverse-complement averaging must also permute tracks (see _load_strand_pair).
    "enformer": {"hf": "EleutherAI/enformer-official-rough",
                 "seq_len": 196_608, "n_bins": 896, "channels_first": False, "autocast": False,
                 "stranded": False, "n_tracks": 5_313},
    "borzoi":   {"hf": "johahi/borzoi-replicate-0",
                 "seq_len": 524_288, "n_bins": 6_144, "channels_first": True, "autocast": True,
                 "stranded": True, "n_tracks": 7_611},
}
CENTER_BINS = 5                 # central bins to summarise the effect over (variant sits at centre)
_BASE = {"A": 0, "C": 1, "G": 2, "T": 3}


def _load_targets(path: str) -> dict[int, str]:
    """Optional index -> assay description map (model's targets_human.txt)."""
    out: dict[int, str] = {}
    try:
        with open(path) as fh:
            for i, row in enumerate(csv.DictReader(fh, delimiter="\t")):
                out[i] = (row.get("description") or row.get("identifier") or "").strip()
    except OSError:
        pass
    return out


def _load_strand_pair(path: str, n_tracks: int) -> list[int] | None:
    """index → index of that track's opposite-strand partner, from the targets file's
    `strand_pair` column. None when the file has no such column.

    Needed only for stranded models. Borzoi's RNA-seq/CAGE targets come in +/- pairs, so
    reverse-complementing the input sequence does not merely reverse the bins — it also swaps
    which member of each pair a given track index reports. Averaging forward with a raw RC
    prediction therefore averages each stranded track against its OPPOSITE strand, which
    cancels real signal and biases delta_max toward zero on exactly the assays Borzoi is used
    for. Enformer's human head is unstranded, so it needs no permutation.
    """
    try:
        with open(path) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
    except OSError:
        return None
    if not rows or "strand_pair" not in rows[0]:
        return None
    pair = list(range(n_tracks))
    for i, row in enumerate(rows):
        if i >= n_tracks:
            break
        try:
            j = int(row["strand_pair"])
        except (TypeError, ValueError):
            return None
        if 0 <= j < n_tracks:
            pair[i] = j
    return pair


def _read_job(path: str):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 3:
                yield row[0], row[1].upper(), row[2].upper()


def _encode(seqs: list[str], channels_first: bool, device) -> torch.Tensor:
    """ACGTN strings -> one-hot. N (and any non-ACGT) maps to the all-zero column, as both
    model ports expect. Returns (B, L, 4) or (B, 4, L) per the model's layout."""
    idx = torch.tensor([[_BASE.get(c, -1) for c in s] for s in seqs])   # (B, L), -1 for N
    oh = F.one_hot(idx.clamp(min=0), 4).float()
    oh[idx < 0] = 0.0
    return oh.permute(0, 2, 1).contiguous().to(device) if channels_first else oh.to(device)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", choices=list(MODELS), default="enformer")
    ap.add_argument("--job", required=True, help="gzipped TSV: id, ref_seq, alt_seq")
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch-size", type=int, default=4)
    ap.add_argument("--targets", default="/cache/targets_human.txt")
    ap.add_argument("--no-rc", action="store_true", help="skip reverse-complement averaging")
    args = ap.parse_args()
    spec = MODELS[args.model]

    if not torch.cuda.is_available():
        print("FATAL: CUDA not available in worker", file=sys.stderr)
        return 2
    device = "cuda"
    targets = _load_targets(args.targets)

    t0 = time.time()
    if args.model == "enformer":
        from enformer_pytorch import from_pretrained
        model = from_pretrained(spec["hf"]).to(device).eval()
    else:
        from borzoi_pytorch import Borzoi
        model = Borzoi.from_pretrained(spec["hf"]).to(device).eval()
    print(f"[worker] {args.model} loaded in {time.time()-t0:.1f}s on "
          f"{torch.cuda.get_device_name(0)}", file=sys.stderr)

    n_bins, cf = spec["n_bins"], spec["channels_first"]

    @torch.no_grad()
    def predict(oh: torch.Tensor) -> torch.Tensor:
        """-> (B, bins, tracks), uniform across models."""
        with torch.autocast("cuda", dtype=torch.float16, enabled=spec["autocast"]):
            out = model(oh)
        h = out["human"] if args.model == "enformer" else out.transpose(1, 2)  # borzoi: tracks-first
        return h.float()

    # Resolved just below, before `human` is ever called; declared here so the closure over
    # it is obvious rather than depending on statement order further down.
    strand_pair = None

    def human(seqs: list[str], rc: bool) -> torch.Tensor:
        oh = _encode(seqs, cf, device)
        pred = predict(oh)
        if rc:                            # average with the reverse-complement (robustness)
            oh_rc = (oh.flip(2)[:, [3, 2, 1, 0], :] if cf
                     else oh.flip(1)[..., [3, 2, 1, 0]])
            rev = predict(oh_rc).flip(1)                   # flip bins (dim 1) back
            if strand_pair is not None:
                # ...and swap each stranded track with its partner (dim 2). Reverse-
                # complementing the input also flips which strand a track index reports;
                # without this the average pairs every +-strand track with its own
                # --strand counterpart. Unstranded models get strand_pair=None and skip it.
                rev = rev.index_select(2, strand_pair)
            pred = (pred + rev) / 2
        return pred

    jobs = list(_read_job(args.job))
    lo = n_bins // 2 - CENTER_BINS // 2
    hi, center = lo + CENTER_BINS, n_bins // 2
    rc = not args.no_rc

    # Stranded models need the track permutation above to average forward with RC correctly.
    # If we can't build it, averaging would quietly corrupt the very assays this model is for,
    # so drop RC averaging rather than produce a biased number and say why.
    if rc and spec["stranded"]:
        pairs = _load_strand_pair(args.targets, spec["n_tracks"])
        if pairs is None:
            rc = False
            print(f"[worker] {args.model} is strand-specific but {args.targets} has no "
                  f"strand_pair column — reverse-complement averaging DISABLED (averaging "
                  f"without the track swap would cancel stranded signal). Supply the model's "
                  f"own targets file to re-enable it.", file=sys.stderr)
        elif pairs == list(range(spec["n_tracks"])):
            print(f"[worker] {args.model}: strand_pair column is the identity — treating the "
                  f"head as unstranded.", file=sys.stderr)
        else:
            strand_pair = torch.tensor(pairs, device=device)
            print(f"[worker] {args.model}: RC averaging will swap "
                  f"{sum(1 for i, j in enumerate(pairs) if i != j)} stranded track pair(s).",
                  file=sys.stderr)

    t0 = time.time()
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["id", "delta_max", "l2_center", "top_track", "top_delta", "top_desc"])
        for b in range(0, len(jobs), args.batch_size):
            chunk = jobs[b:b + args.batch_size]
            ids = [c[0] for c in chunk]
            for s in (s for c in chunk for s in (c[1], c[2])):
                if len(s) != spec["seq_len"]:
                    raise SystemExit(f"seq length {len(s)} != {spec['seq_len']} ({args.model})")
            delta = human([c[2] for c in chunk], rc) - human([c[1] for c in chunk], rc)  # (B,bins,tracks)
            win = delta[:, lo:hi, :]
            for i, vid in enumerate(ids):
                flat = win[i].abs()
                pos = torch.argmax(flat)
                track = int(pos % delta.shape[2])
                w.writerow([vid, f"{flat.max().item():.5f}",
                            f"{delta[i, center].pow(2).sum().sqrt().item():.5f}",
                            track, f"{win[i].flatten()[pos].item():.5f}", targets.get(track, "")])
            print(f"[worker] {min(b+args.batch_size, len(jobs))}/{len(jobs)} variants "
                  f"({time.time()-t0:.1f}s)", file=sys.stderr)
    print(f"[worker] done: {len(jobs)} variants in {time.time()-t0:.1f}s", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
