#!/usr/bin/env python
"""canonical_scores_ranked.csv -- the canonical barcodes ranked by basal score.

Attaches the hand-curated constitutive-activity evidence from
literature_constitutive_activity.csv (top 40 only -- a blank ca_evidence means
nobody looked, not that no evidence exists).

Reads:  data/processed/basal_scores_per_barcode.csv, coupling_annotation.csv,
        literature_constitutive_activity.csv
Writes: data/processed/canonical_scores_ranked.csv

Table-only. Extracted from plot_basal.py, whose plotting code is not part of this
release; the columns are that script's, less the ones no published panel reads.

Usage: python scripts/build_ranked_tables.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "data" / "processed"


def main():
    d = pd.read_csv(HERE / "basal_scores_per_barcode.csv")
    # corrected coupling: see build_coupling.py. Keyed on Barcode so isoforms inherit
    # their own construct's annotation rather than a gene-level guess.
    _ann = pd.read_csv(HERE / "coupling_annotation.csv")[
        ["Barcode", "coupling_display", "gs_primary", "gs_any"]]
    d = d.merge(_ann, on="Barcode", how="left")
    d["coupling"] = d["coupling_display"].fillna("unknown")
    d[["gs_primary", "gs_any"]] = d[["gs_primary", "gs_any"]].fillna(False)
    s = d[d["n_days"] == 2].sort_values("score", ascending=False).reset_index(drop=True)

    # ---- canonical barcodes, ranked ----------------------------------------
    can = s[s["isoform"] == "canonical"].reset_index(drop=True)
    # literature annotation: curated for the top 40 only. A blank ca_evidence
    # means "not reviewed here", NOT "no prior evidence exists".
    lit = pd.read_csv(HERE / "literature_constitutive_activity.csv")
    can = can.merge(lit, on="symbol", how="left")
    can["ca_evidence"] = can["ca_evidence"].fillna("not reviewed")
    can["rank"] = np.arange(1, len(can) + 1)
    cols = ["rank", "symbol", "coupling", "score", "score_sd", "total_reads",
            "ca_evidence", "ca_mechanism", "ca_caveat", "ca_source", "transcript_id", "Barcode",
            "gs_primary", "gs_any"]
    can[cols].to_csv(HERE / "canonical_scores_ranked.csv", index=False, float_format="%.4g")

    n_zoom = 40
    top = can.head(n_zoom)
    print(f"wrote canonical_scores_ranked.csv ({len(can)} canonical barcodes "
          f"over {can.symbol.nunique()} distinct receptors)")
    print(f"top {n_zoom} canonical: "
          f"{(top['ca_evidence'] == 'established').sum()} established, "
          f"{(top['ca_evidence'] == 'reported').sum()} reported, "
          f"{(top['ca_evidence'].isin(['weak', 'none found'])).sum()} weak/none, "
          f"{(top['ca_evidence'] == 'not assessed').sum()} not assessed; "
          f"{top['ca_caveat'].notna().sum()} carry a caveat")


if __name__ == "__main__":
    main()
