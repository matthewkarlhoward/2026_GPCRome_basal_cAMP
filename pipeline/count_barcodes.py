import csv
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import Counter


def count_barcodes(fastq_file, mapping_file, annotation_file, output_file, mismatch_tolerance=1):
    # Read the barcode mapping (no header: barcode_sequence, uniprot_id)
    barcode_mapping = {}
    with open(mapping_file, 'r', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        for row in reader:
            barcode, barcode_id = row
            barcode_mapping[barcode] = [barcode_id, 0]

    # Function to find best matching barcode
    def find_best_match(sequence, barcodes, tolerance):
        best_match = None
        best_score = -float('inf')
        for barcode in barcodes:
            num_mismatches = sum(c1 != c2 for c1, c2 in zip(sequence, barcode))
            if num_mismatches <= tolerance:
                score = len(barcode) - num_mismatches
                if score > best_score:
                    best_score = score
                    best_match = barcode
        return best_match

    # Count barcodes in the fastq file. Tally identical read sequences first,
    # then match each unique sequence once (same result, far fewer comparisons).
    read_counts = Counter()
    total_reads = 0
    with open(fastq_file, 'r') as f:
        for _, seq, _ in FastqGeneralIterator(f):
            total_reads += 1
            read_counts[seq] += 1

    barcode_counts = Counter()
    for seq, n in read_counts.items():
        best_match = find_best_match(seq, barcode_mapping.keys(), mismatch_tolerance)
        if best_match:
            barcode_counts[best_match] += n

    # Update counts in mapping
    for barcode, count in barcode_counts.items():
        barcode_mapping[barcode][1] = count

    # Build counts dataframe
    counts_df = pd.DataFrame(
        [(barcode, barcode_id, count) for barcode, (barcode_id, count) in barcode_mapping.items()],
        columns=['Barcode', 'BarcodeID', 'Count']
    )

    # Join with per-barcode annotation (symbol, transcript, isoform flag)
    # Each barcode is unique so this join is 1-to-1 — no deduplication needed
    annot_df = pd.read_csv(annotation_file)  # columns: Barcode, symbol, BarcodeID, transcript_id, isoform
    result = counts_df.merge(annot_df[['Barcode', 'symbol', 'transcript_id', 'isoform']], on='Barcode', how='left')
    result = result[['symbol', 'BarcodeID', 'transcript_id', 'isoform', 'Barcode', 'Count']]
    result = result.sort_values('Count', ascending=False).reset_index(drop=True)

    result.to_csv(output_file, index=False)

    total_matches = sum(barcode_counts.values())
    print(f"Total reads:          {total_reads:,}")
    print(f"Matched reads:        {total_matches:,} ({100*total_matches/total_reads:.1f}%)")
    print(f"Unmatched reads:      {total_reads - total_matches:,}")
    print(f"Results written to:   {output_file}")


if __name__ == "__main__":
    count_barcodes(
        fastq_file='data/your_reads.fastq',
        mapping_file='reference/barcode_mapping.csv',
        annotation_file='reference/barcode_annotation.csv',
        output_file='output/barcode_counts.csv',
        mismatch_tolerance=1
    )
