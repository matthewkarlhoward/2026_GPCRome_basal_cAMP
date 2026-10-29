"""QC for the 20261015 GPCRome basal barcode run.

For each trimmed FASTQ, computes every read's Hamming distance to its nearest
reference barcode, then reports:
  - qc_summary.csv          per-sample reads, match rates at 0 / <=1 / <=2 mismatches,
                            barcodes detected, read-length problems
  - top_unmatched.csv       most common reads >1 mismatch from every barcode
  - tol0_vs_tol1_per_barcode.csv  per-barcode counts, exact vs <=1 mismatch

Only 10-bp reads are scored here; other lengths are tallied separately
(count_barcodes.py scores them on a truncated prefix).

Usage: python barcode_qc.py <count_barcodes_tool dir> <out dir>
"""
import glob
import os
import re
import sys
from collections import Counter

import numpy as np
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator

TOOL, OUT = sys.argv[1], sys.argv[2]
SAMPLES = [f"s{i}_basal_S{i}" for i in range(5, 13)]   # the two replicate sorts
LUT = np.zeros(256, dtype=np.uint8)
for i, c in enumerate("ACGTN"):
    LUT[ord(c)] = i


def encode(seqs):
    return LUT[np.frombuffer("".join(seqs).encode(), dtype=np.uint8)].reshape(len(seqs), 10)


ref = pd.read_csv(os.path.join(TOOL, "reference/barcode_mapping.csv"),
                  header=None, names=["Barcode", "BarcodeID"], encoding="utf-8-sig")
annot = pd.read_csv(os.path.join(TOOL, "reference/barcode_annotation.csv"))[["Barcode", "symbol"]]
bc = list(ref["Barcode"])
bc_arr = encode(bc)

# Minimum pairwise distance within the reference: tells us how safe 1-2 mismatches are
pair = (bc_arr[:, None, :] != bc_arr[None, :, :]).sum(2)
np.fill_diagonal(pair, 99)
print(f"Reference: {len(bc)} barcodes, min pairwise Hamming distance = {pair.min()}, "
      f"pairs at distance <=2: {(pair <= 2).sum() // 2}")


def cutadapt_stats(sample):
    txt = open(os.path.join(TOOL, f"logs/{sample}_cutadapt.log")).read()
    get = lambda pat: int(re.search(pat, txt).group(1).replace(",", ""))
    return get(r"Total reads processed:\s+([\d,]+)"), get(r"Reads written \(passing filters\):\s+([\d,]+)")


rows, unmatched_rows, per_bc = [], [], {}
for s in SAMPLES:
    counts = Counter()
    with open(os.path.join(TOOL, f"trimmed/{s}_trimmed.fastq")) as f:
        for _, seq, _ in FastqGeneralIterator(f):
            counts[seq] += 1
    seqs10 = [q for q in counts if len(q) == 10]
    n10 = np.array([counts[q] for q in seqs10])
    other_len = sum(n for q, n in counts.items() if len(q) != 10)

    # nearest-barcode distance, chunked to keep memory small
    dmin = np.empty(len(seqs10), dtype=np.int64)
    best = np.empty(len(seqs10), dtype=np.int64)
    nties = np.empty(len(seqs10), dtype=np.int64)
    arr = encode(seqs10)
    for i in range(0, len(seqs10), 5000):
        d = (arr[i:i + 5000, None, :] != bc_arr[None, :, :]).sum(2)
        dmin[i:i + 5000] = d.min(1)
        best[i:i + 5000] = d.argmin(1)  # first in file order, as count_barcodes.py does
        nties[i:i + 5000] = (d == d.min(1, keepdims=True)).sum(1)

    total, _ = cutadapt_stats(s)
    trimmed = sum(counts.values())
    by_d = {k: int(n10[dmin == k].sum()) for k in (0, 1, 2)}
    tol0 = np.bincount(best[dmin == 0], weights=n10[dmin == 0], minlength=len(bc))
    tol1 = np.bincount(best[dmin <= 1], weights=n10[dmin <= 1], minlength=len(bc))
    per_bc[f"{s}_tol0"], per_bc[f"{s}_tol1"] = tol0.astype(int), tol1.astype(int)
    rows.append({
        "sample": s,
        "raw_pairs": total,
        "trimmed_reads": trimmed,
        "pct_trimmed": 100 * trimmed / total,
        "pct_exact": 100 * by_d[0] / trimmed,
        "pct_le1mm": 100 * (by_d[0] + by_d[1]) / trimmed,
        "pct_le2mm": 100 * (by_d[0] + by_d[1] + by_d[2]) / trimmed,
        "pct_gt2mm": 100 * n10[dmin > 2].sum() / trimmed,
        "pct_not_10bp": 100 * other_len / trimmed,
        "reads_rescued_by_1mm": by_d[1],
        "ambiguous_1mm_reads": int(n10[(dmin == 1) & (nties > 1)].sum()),
        "barcodes_detected_tol0": int((tol0 > 0).sum()),
        "barcodes_detected_tol1": int((tol1 > 0).sum()),
        "barcodes_ge10_reads_tol1": int((tol1 >= 10).sum()),
    })

    um = dmin > 1
    order = np.argsort(-n10[um])[:25]
    um_seqs, um_n, um_d = np.array(seqs10)[um][order], n10[um][order], dmin[um][order]
    for q, n, d in zip(um_seqs, um_n, um_d):
        unmatched_rows.append({"sample": s, "sequence": q, "reads": int(n),
                               "pct_of_trimmed": 100 * n / trimmed, "nearest_bc_distance": int(d)})
    print(f"{s}: done")

summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(OUT, "qc_summary.csv"), index=False, float_format="%.2f")
pd.DataFrame(unmatched_rows).to_csv(os.path.join(OUT, "top_unmatched.csv"), index=False, float_format="%.3f")
pb = pd.concat([ref, pd.DataFrame(per_bc)], axis=1).merge(annot, on="Barcode", how="left")
pb.to_csv(os.path.join(OUT, "tol0_vs_tol1_per_barcode.csv"), index=False)
pd.set_option("display.width", 250)
print(summary.round(1).to_string(index=False))
