# 20261015 GPCRome basal — barcode counting QC

> **Scope.** Two replicate sorts, four bins each: s5–s12.

**Pipeline:** `pipeline/process_one.sh` (fastp merge, then cutadapt linked-primer trim, then `count_barcodes.py`, ≤1 mismatch).
**Count tables:** `../data/counts/<sample>_counts.csv`.
`count_barcodes.py` tallies unique reads before matching (~20× faster), which gives output identical to matching read by read.

## Headline numbers (`qc_summary.csv`)

| Samples | Trimmed | Exact | ≤1 mm (pipeline) | ≤2 mm |
|---|---|---|---|---|
| s5–s12 | ~99% | 86–90% | 93–96% | 97–98% |

- **1 mismatch is safe:** the closest pair of reference barcodes is 4 mismatches apart, and no read was 1 mismatch from two barcodes. Allowing 1 mismatch recovers 5–7 percentage points of reads and 4–46 more barcodes detected per sample.
- **2 mismatches would add another ~3 points, but it is not provably unambiguous.** A
minimum pairwise distance of 4 rules out a read lying ≤1 from two barcodes (that would
require distance ≤2), but not one lying exactly 2 from each of a pair at distance 4,
and 9,572 such pairs exist. Raising the tolerance would need the ambiguous-read count
checked, not assumed.
