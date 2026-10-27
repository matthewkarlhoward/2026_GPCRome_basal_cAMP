# 20261015 GPCRome basal — barcode counting QC

> **Scope.** The tables in `../data/qc/` cover only the two retained replicates (s5–s12).
> The numbers below for s1–s4 and the s13 baseline are kept as the record of why those
> samples were excluded; their count and QC rows are not part of this release.
> See `METHODS.md`.

**Pipeline:** `pipeline/process_one.sh` (fastp merge, then cutadapt linked-primer trim, then `count_barcodes.py`, ≤1 mismatch).
**Count tables:** `../data/counts/<sample>_counts.csv`.
`count_barcodes.py` tallies unique reads before matching (~20× faster), which gives output identical to matching read by read.

## Headline numbers (`qc_summary.csv`)

| Samples | Trimmed | Exact | ≤1 mm (pipeline) | ≤2 mm |
|---|---|---|---|---|
| s5–s12 | ~99% | 86–90% | 93–96% | 97–98% |
| s1 / s3 | 99% | 79 / 74% | 85 / 80% | 88 / 84% |
| s2 / s4 | 99% | 66 / 67% | 71 / 73% | 75 / 77% |
| s13 baseline | 99% | 40% | 43% | 47% |

- **1 mismatch is safe:** the closest pair of reference barcodes is 4 mismatches apart, and no read was 1 mismatch from two barcodes. Allowing 1 mismatch recovers 5–7 percentage points of reads and 4–46 more barcodes detected per sample.
- **2 mismatches would add another ~3 points, but it is not provably unambiguous.** A
minimum pairwise distance of 4 rules out a read lying ≤1 from two barcodes (that would
require distance ≤2), but not one lying exactly 2 from each of a pair at distance 4,
and 9,572 such pairs exist. Raising the tolerance would need the ambiguous-read count
checked, not assumed.

## Low match rate in s1–s4 and s13: one off-library barcode

`TTTTGCGCTG` is **not in the 946-barcode reference** and appeared in no earlier run. Its exact-match read counts:

| s1 | s2 | s3 | s4 | s5–s12 | s13 |
|---|---|---|---|---|---|
| 175k (9%) | 376k (20%) | 259k (13%) | 364k (18%) | 21–39 reads (index-hopping level) | 817k (45%) |

Together with its sequencing-error variants it accounted for nearly all the unmatched reads in those samples. The amplicon backbone is identical to the library; only the barcode differs. So it's a barcoded plasmid in the same vector that got into those 5 samples (a spiked control, a single clone, or contamination), not a PCR artifact.

**Impact:** s1–s4 still have 1.3–1.9M library reads. Normalizing within matched (library) reads, which `Count` already is, removes the effect on GPCR proportions. The identity of that construct is not established; those samples are excluded from all
analysis regardless.
