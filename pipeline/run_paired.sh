#!/bin/bash
set -euo pipefail

# Pipeline for paired-end fastq.gz input:
#   1. fastp merges overlapping R1/R2 into a single consensus read
#   2. cutadapt extracts the 10 bp barcode via linked primer anchors
#   3. count_barcodes.py counts and annotates against the reference
#
# Usage:
#   bash run_paired.sh <input_root>
# where <input_root> contains subfolders per sample, each holding
# *_R1_001.fastq.gz and *_R2_001.fastq.gz.

INPUT_ROOT="${1:-../matt_howard}"
TOOL_DIR="$(cd "$(dirname "$0")" && pwd)"
PY="${PY:-python}"
CUTADAPT="${CUTADAPT:-cutadapt}"

mkdir -p "$TOOL_DIR/merged" "$TOOL_DIR/trimmed" "$TOOL_DIR/output" "$TOOL_DIR/logs"

cd "$TOOL_DIR"

shopt -s nullglob
for r1 in "$INPUT_ROOT"/*/*_R1_001.fastq.gz; do
  r2="${r1/_R1_/_R2_}"
  base="$(basename "$r1" _R1_001.fastq.gz)"
  echo "=== $base ==="

  merged="merged/${base}_merged.fastq.gz"
  trimmed="trimmed/${base}_trimmed.fastq"

  fastp \
    -i "$r1" -I "$r2" \
    --merge --merged_out "$merged" \
    --include_unmerged \
    --disable_adapter_trimming \
    --disable_quality_filtering \
    --disable_length_filtering \
    -j "logs/${base}_fastp.json" \
    -h "logs/${base}_fastp.html" \
    2> "logs/${base}_fastp.log"

  "$CUTADAPT" \
    -g "GGGTTAGGCG...CTTCGCGATG" \
    -g "CATCGCGAAG...CGCCTTACCC" \
    --discard-untrimmed \
    -o "$trimmed" \
    "$merged" \
    > "logs/${base}_cutadapt.log" 2>&1

  "$PY" -c "
from count_barcodes import count_barcodes
count_barcodes(
    fastq_file='${trimmed}',
    mapping_file='reference/barcode_mapping.csv',
    annotation_file='reference/barcode_annotation.csv',
    output_file='output/${base}_counts.csv',
    mismatch_tolerance=1,
)
"
done

echo ""
echo "Done. Results in output/"
