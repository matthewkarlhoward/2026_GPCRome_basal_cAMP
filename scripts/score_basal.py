"""Weighted mean bin (WMB) per barcode per sort day, and overall basal score.

Design: 2 replicate sorts x 4 bins (a-d).
    replicate 1: 20260924, s5-s8    replicate 2: 20260926, s9-s12

Scoring unit is the barcode (one construct: canonical or alt isoform), never the gene.

Per day, each bin's counts are converted to frequency within that bin's library reads
(reads matched to a reference barcode), so sequencing depth cancels out. This assumes
an equal number of cells
was sorted into each bin.
    WMB_day = sum_b(w_b * f_b) / sum_b(f_b),   w = a:1 b:2 c:3 d:4
The reported `score` is that WMB rescaled linearly onto 0-1 by (WMB - 1) / 3, so 0 = all
signal in the lowest bin and 1 = all in the highest. This is a change of units only: it
adds no information and the endpoints are the bin bounds, not a biological zero or max.
NOTE the naming: `<day>_wmb` is the RESCALED 0-1 value, not the raw 1-4 one.
A day counts for a barcode only if its raw reads summed over the 4 bins are >= MIN_READS.
    score = mean of passing days; also SD, n_days.

Reads:  data/counts/s{5..12}_basal_S{n}_counts.csv
Writes: data/processed/basal_scores_per_barcode.csv

Usage: python scripts/score_basal.py   (paths are relative to this file)
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
COUNTS, OUT = ROOT / "data" / "counts", ROOT / "data" / "processed"
DAYS = {"20260924": [5, 6, 7, 8], "20260926": [9, 10, 11, 12]}
assert len(DAYS) == 2, "the screen has two replicate sorts"
BINS, W = "abcd", np.array([1, 2, 3, 4])
WMB_LO, WMB_HI = 1.0, 4.0          # bin-index bounds the raw WMB can take
rescale = lambda v: (v - WMB_LO) / (WMB_HI - WMB_LO)
MIN_READS = 50
KEY = ["symbol", "BarcodeID", "transcript_id", "isoform", "Barcode"]


def load(n):
    return pd.read_csv(COUNTS / f"s{n}_basal_S{n}_counts.csv").set_index("Barcode")


def main():
    out = load(5)[KEY[:-1]].sort_index()

    wmb_cols = []
    for day, samples in DAYS.items():
        raw = np.column_stack([load(n)["Count"].reindex(out.index).to_numpy() for n in samples])
        freq = raw / raw.sum(0, keepdims=True)
        for j, b in enumerate(BINS):
            out[f"{day}_{b}"] = raw[:, j]
        reads = raw.sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            wmb = (freq * W).sum(1) / freq.sum(1)
        out[f"{day}_reads"] = reads
        out[f"{day}_wmb"] = rescale(np.where(reads >= MIN_READS, wmb, np.nan))
        wmb_cols.append(f"{day}_wmb")

    m = out[wmb_cols]
    out["n_days"] = m.notna().sum(1)
    out["score"] = m.mean(1)
    out["score_sd"] = m.std(1, ddof=1)
    out["total_reads"] = out[[f"{d}_reads" for d in DAYS]].sum(1)
    out = out.sort_values(["score", "total_reads"], ascending=[False, False])
    out.reset_index().to_csv(OUT / "basal_scores_per_barcode.csv", index=False, float_format="%.4g")

    print(f"{len(out)} barcodes; scored on >=1 day: {(out.n_days >= 1).sum()}, "
          f"both days: {(out.n_days == 2).sum()}  (days used: {list(DAYS)})")
    print(out.loc[out.n_days == 2, ["symbol", "isoform", "transcript_id", "score", "score_sd", "n_days",
                                     "total_reads"]].head(20).to_string())


if __name__ == "__main__":
    main()
