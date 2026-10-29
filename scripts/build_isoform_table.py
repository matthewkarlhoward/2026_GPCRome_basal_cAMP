#!/usr/bin/env python3
"""Alternative isoforms against their canonical construct.

One row per alternative-isoform construct: its own basal score, the score of the
canonical construct for the same gene, and the difference. Supplemental Figure 2
is built from this table.

Provenance note. This table used to fall out of an exploratory script
(`plot_basal.py`) in the pre-repo working directory, as a side effect of drawing
an unrelated panel. It is rebuilt here from `receptor_construct_features.csv`
alone so that the supplement has a declared input rather than depending on a
figure script having been run first. Genes whose canonical construct did not
score are dropped, since there is nothing to compare against.

Usage: python build_isoform_table.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "processed" / "receptor_construct_features.csv"
OUT = ROOT / "data" / "processed" / "alt_vs_canonical_isoform.csv"

COLS = ["symbol", "transcript_id", "Barcode", "score", "canonical", "delta",
        "total_reads"]


def main():
    f = pd.read_csv(SRC)
    # A few genes carry two canonical constructs under different barcodes, so
    # "the canonical score" is ambiguous for them. Average the duplicates, which
    # is what the original exploratory build did; among the genes that also have
    # alternative isoforms this affects RXFP2 and GRM8.
    canon = (f[f.isoform == "canonical"]
             .dropna(subset=["score"])
             .groupby("symbol")["score"].mean())
    alt = f[(f.isoform != "canonical") & f.score.notna()].copy()
    alt["canonical"] = alt.symbol.map(canon)
    alt = alt.dropna(subset=["canonical"])
    alt["delta"] = alt.score - alt.canonical
    alt = alt.sort_values("delta")[COLS]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    alt.to_csv(OUT, index=False, float_format="%.4g")
    print(f"wrote {OUT}  ({len(alt)} alternative isoforms, "
          f"{alt.symbol.nunique()} genes)")
    print(f"  delta: median {alt.delta.median():+.3f}, "
          f"{(alt.delta < 0).sum()} below canonical, "
          f"{(alt.delta > 0).sum()} above")


if __name__ == "__main__":
    main()
