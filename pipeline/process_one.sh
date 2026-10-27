#!/bin/bash
# Process a single paired-end sample: fastp merge -> cutadapt -> count_barcodes.
# Usage: process_one.sh <R1.fastq.gz>
# R2 is derived by replacing _R1_ with _R2_.
set -euo pipefail

r1="$1"
r2="${r1/_R1_/_R2_}"
base="$(basename "$r1" _R1_001.fastq.gz)"

TOOL_DIR="$(cd "$(dirname "$0")" && pwd)"
PY="${PY:-python}"
CUTADAPT="${CUTADAPT:-cutadapt}"

cd "$TOOL_DIR"

merged="merged/${base}_merged.fastq.gz"
trimmed="trimmed/${base}_trimmed.fastq"
counts="output/${base}_counts.csv"

# Skip if already completed
if [[ -s "$counts" ]]; then
  echo "[$base] already done — skipping"
  exit 0
fi

echo "[$base] START $(date +%H:%M:%S)"

fastp \
  -i "$r1" -I "$r2" \
  --merge --merged_out "$merged" \
  --include_unmerged \
  --disable_adapter_trimming \
  --disable_quality_filtering \
  --disable_length_filtering \
  --thread 2 \
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
    output_file='${counts}',
    mismatch_tolerance=1,
)
" > "logs/${base}_count.log" 2>&1

match=$(grep "Matched reads:" "logs/${base}_count.log" || true)
echo "[$base] DONE  $(date +%H:%M:%S)  ${match}"
