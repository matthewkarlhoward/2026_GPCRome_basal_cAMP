# Barcode counting pipeline

FASTQ → per-barcode counts. This is what produced `../data/counts/`.

```bash
bash process_one.sh <sample>_R1_001.fastq.gz          # one sample; R2 is derived
bash run_paired.sh  <dir-of-per-sample-folders>       # the whole run
```

Both run three steps:

1. **fastp** merges R1/R2 into a single consensus read. Reads are full-length ~695 bp
   amplicons, so the pair overlaps; adapter trimming, quality filtering and length
   filtering are all disabled, and unmerged pairs are kept.
2. **cutadapt** pulls the 10 bp barcode out using linked primer anchors, in both strand
   orientations, discarding reads where the anchors are not found (~1%):

   ```
   -g "GGGTTAGGCG...CTTCGCGATG"   -g "CATCGCGAAG...CGCCTTACCC"
   ```

3. **`count_barcodes.py`** matches the barcode to the 946-barcode reference allowing
   ≤1 mismatch, annotates it with its gene symbol, and writes one row per reference
   barcode sorted by count. Barcodes with zero counts are retained so dropouts are
   visible. It tallies unique reads before matching, which is ~20× faster than matching
   read by read and gives identical output.

`fastp` and `cutadapt` must be on `PATH`; set `PY` and `CUTADAPT` to override the
interpreter and binary. Count tables land in `output/`; the copies in `../data/counts/`
are those files.

Match rate on the retained samples is 93–96% of trimmed reads. One mismatch is
unambiguous: the closest pair of reference barcodes is 4 mismatches apart and no read in
this run was 1 mismatch from two barcodes. `barcode_qc.py` produced the tables in
`../data/qc/`; see `../docs/SEQUENCING_QC.md`.

## Reference files

`reference/barcode_mapping.csv` — two columns, **no header**: barcode sequence → UniProt
accession: 946 barcodes over 766 distinct UniProt accessions (one receptor can appear
as more than one construct).

```
CAGCAGAGAT,O14842
GCTCCAGTAT,O14843
```

`reference/barcode_annotation.csv` — the fuller per-barcode annotation that
`count_barcodes.py` joins onto the output. Note that the column named `BarcodeID`, here and
in every counts table, holds a **UniProt accession**, not a barcode.

## Notes

- Do not run `count_barcodes.py` on raw reads — match rates will be <1%, because the reads
  are full-length amplicons, not barcodes.
- Reads that match no barcode within the tolerance are discarded.
- A read equally close to two barcodes is resolved to whichever appears **first in the
  reference file** (the match loop keeps a strict `>` on the score). At tolerance 1 this
  never arises in this run: `ambiguous_1mm_reads` is 0 in all eight samples.
