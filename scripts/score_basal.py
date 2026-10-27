"""Weighted mean bin (WMB) per barcode per sort day, and overall basal score.

Design (from the bench sheet): 3 sort days x 4 bins (a-d), plus s13 = unsorted baseline.
    day 20260922: s1-s4   day 20260924: s5-s8   day 20260926: s9-s12
    (sheet lists the third day as 20240926; assumed to be 2026)
Day 20260922 (s1-s4) and the baseline (s13) are excluded from all analysis: the day-1
bins have near-identical barcode composition (no sort separation; WMB r = -0.2 vs the
other days, no bin relabelling rescues it), and both carry the off-library contaminant.

Scoring unit is the barcode (one construct: canonical or alt isoform), never the gene.

Per day, each bin's counts are converted to frequency within that bin's library reads
(reads matched to a reference barcode), so sequencing depth and the off-library
contaminant barcode in s1-s4 / s13 cancel out. This assumes an equal number of cells
was sorted into each bin. If bin cell counts are known, set CELLS below.
    WMB_day = sum_b(w_b * f_b) / sum_b(f_b),   w = a:1 b:2 c:3 d:4
The reported `score` is that WMB rescaled linearly onto 0-1 by (WMB - 1) / 3, so 0 = all
signal in the lowest bin and 1 = all in the highest. This is a change of units only: it
adds no information and the endpoints are the bin bounds, not a biological zero or max.
The raw 1-4 values are kept as `<day>_wmb_raw14`, so nothing is lost. NOTE the naming:
`<day>_wmb` is the RESCALED 0-1 value, not the raw one.
A day counts for a barcode only if its raw reads summed over the 4 bins are >= MIN_READS.
    score = mean of passing days; also SD, n_days.

A WMB weighted by raw counts rather than within-bin frequencies is kept alongside, as
`score_rawcounts`, for comparison (the two agree at r = 0.9995).

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
BINS, W = "abcd", np.array([1, 2, 3, 4])
WMB_LO, WMB_HI = 1.0, 4.0          # bin-index bounds the raw WMB can take
rescale = lambda v: (v - WMB_LO) / (WMB_HI - WMB_LO)
CELLS = None  # e.g. {"20260922": [n_a, n_b, n_c, n_d], ...}; None = equal cells per bin
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
        if CELLS is not None:
            freq = freq * np.asarray(CELLS[day]) / np.sum(CELLS[day])
        for j, b in enumerate(BINS):
            out[f"{day}_{b}"] = raw[:, j]
        reads = raw.sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            wmb = (freq * W).sum(1) / freq.sum(1)
            wmb_raw = (raw * W).sum(1) / reads
        out[f"{day}_reads"] = reads
        out[f"{day}_wmb_raw14"] = np.where(reads >= MIN_READS, wmb, np.nan)
        out[f"{day}_wmb"] = rescale(out[f"{day}_wmb_raw14"])
        out[f"{day}_wmb_rawcounts"] = rescale(np.where(reads >= MIN_READS, wmb_raw, np.nan))
        wmb_cols.append(f"{day}_wmb")

    m = out[wmb_cols]
    out["n_days"] = m.notna().sum(1)
    out["score"] = m.mean(1)
    out["score_sd"] = m.std(1, ddof=1)
    out["wmb"] = out[[f"{d}_wmb_raw14" for d in DAYS]].mean(1)   # raw 1-4, for reference
    out["score_rawcounts"] = out[[c + "_rawcounts" for c in wmb_cols]].mean(1)
    out["total_reads"] = out[[f"{d}_reads" for d in DAYS]].sum(1)
    out = out.sort_values(["score", "total_reads"], ascending=[False, False])
    out.reset_index().to_csv(OUT / "basal_scores_per_barcode.csv", index=False, float_format="%.4g")

    print(f"{len(out)} barcodes; scored on >=1 day: {(out.n_days >= 1).sum()}, "
          f"both days: {(out.n_days == 2).sum()}  (days used: {list(DAYS)})")
    print(out.loc[out.n_days == 2, ["symbol", "isoform", "transcript_id", "score", "score_sd", "n_days",
                                     "total_reads"]].head(20).to_string())


if __name__ == "__main__":
    main()
