"""Supplementary table: basal cAMP score for every GPCRome construct.

One row per barcode -- the unit of analysis. Canonical and alternative isoforms of a
receptor are separate constructs and are never collapsed to the gene, so 3 genes (GRM8,
RXFP1, RXFP2) appear twice as canonical.

All 946 library constructs are included, so dropouts are visible rather than silently
absent. Filter `n_replicates == 2` to get the 901 constructs the figures use.

Writes: ../tables/supp_table1_GPCRome_basal.csv

Usage: python scripts/make_supp_table.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "data" / "processed"
OUT = ROOT / "tables" / "supp_table1_GPCRome_basal.csv"

D1, D2 = "20260924", "20260926"

# The two curated flags drawn on Figure 1c, carried here so the table and the figure
# cannot disagree. Both sets are duplicated verbatim in
# 2026_GPCRome_basal_cAMP/scripts/plot_figure1.py; the assertion in main() fails if the
# top-30 counts drift apart. `orphan_curated` is NOT the same as `ligand_family ==
# "orphan"`: GPR4 is annotated an orphan there but is a deorphanised proton sensor.
ORPHAN = {"GPR6", "GPR52", "GPR161", "GPR12", "GPR62", "GPR101", "GPR26", "GPR3",
          "GPR21", "GPR61", "GPR85", "GPR173"}
AUTOCRINE = {"PTGDR", "PTGER2", "PTGER4", "PTGIR", "ADORA2A", "ADORA2B",
             "GPBAR1", "P2RY1", "P2RY2", "S1PR1", "FFAR4"}
N_TOP = 30


def main():
    s = pd.read_csv(HERE / "basal_scores_per_barcode.csv")
    coup = pd.read_csv(HERE / "coupling_annotation.csv")[
        ["Barcode", "coupling_display", "gs_primary", "gs_any"]]
    feat = pd.read_csv(HERE / "receptor_construct_features.csv")[
        ["Barcode", "gpcr_class", "ligand_family"]]
    rank = pd.read_csv(HERE / "canonical_scores_ranked.csv")[
        ["Barcode", "rank", "ca_evidence", "ca_mechanism", "ca_caveat", "ca_source"]]
    ted = pd.read_csv(HERE / "basal_vs_tedman_pme.csv")[
        ["Barcode", "tedman_pme", "pme_compressed"]]

    d = (s.merge(coup, on="Barcode", how="left")
          .merge(feat, on="Barcode", how="left")
          .merge(rank, on="Barcode", how="left")
          .merge(ted, on="Barcode", how="left"))

    out = pd.DataFrame({
        "symbol":                     d["symbol"],
        "uniprot_id":                 d["BarcodeID"],
        "ensembl_transcript_id":      d["transcript_id"],
        "isoform":                    d["isoform"],
        "barcode":                    d["Barcode"],
        "basal_camp_score":           d["score"],
        "basal_camp_score_sd":        d["score_sd"],
        "n_replicates":               d["n_days"],
        "rep1_score":                 d[f"{D1}_wmb"],
        "rep2_score":                 d[f"{D2}_wmb"],
        "total_reads":                d["total_reads"],
        "canonical_rank":             d["rank"],
        "coupling":                   d["coupling_display"],
        "gs_primary":                 d["gs_primary"],
        "gs_any":                     d["gs_any"],
        "gpcr_class":                 d["gpcr_class"],
        "ligand_family":              d["ligand_family"],
        "orphan_curated":             d["symbol"].isin(ORPHAN),
        "agonist_in_situ":            d["symbol"].isin(AUTOCRINE),
        "surface_expression_tedman":  d["tedman_pme"],
        "surface_value_saturated":    d["pme_compressed"],
        "prior_constitutive_evidence": d["ca_evidence"],
        "prior_evidence_mechanism":   d["ca_mechanism"],
        "prior_evidence_caveat":      d["ca_caveat"],
        "prior_evidence_source":      d["ca_source"],
    })

    # booleans as TRUE/FALSE rather than a mix of 1.0/True; blank where not applicable
    for c in ("gs_primary", "gs_any", "surface_value_saturated",
              "orphan_curated", "agonist_in_situ"):
        out[c] = out[c].map({True: "TRUE", False: "FALSE", 1: "TRUE", 0: "FALSE",
                             1.0: "TRUE", 0.0: "FALSE"})
    for c in ("n_replicates", "total_reads", "canonical_rank"):
        out[c] = out[c].astype("Int64")

    # highest basal activity first; unscored constructs last
    out = out.sort_values(["basal_camp_score", "symbol"],
                          ascending=[False, True], na_position="last")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False, float_format="%.4g")

    scored = out["n_replicates"].eq(2).sum()
    print(f"wrote {OUT}")
    print(f"  {len(out)} constructs | {scored} scored in both replicates "
          f"({out.loc[out.n_replicates.eq(2), 'isoform'].value_counts().to_dict()})")
    print(f"  {out['canonical_rank'].notna().sum()} canonical ranked | "
          f"{out['prior_constitutive_evidence'].notna().sum()} with curated literature | "
          f"{out['surface_expression_tedman'].notna().sum()} with a Tedman surface value")
    print(f"  top 5: {', '.join(out['symbol'].head(5))}")

    # ---- Figure 1c reconciliation ----------------------------------------
    # These are the numbers drawn in the Fig. 1c keys. If this block stops matching the
    # figure, one of the two copies of ORPHAN/AUTOCRINE has drifted.
    top = out[out["canonical_rank"].notna()].nsmallest(N_TOP, "canonical_rank")
    grade = top["prior_constitutive_evidence"].value_counts()
    n_rep = int(top["prior_constitutive_evidence"].isin(["established", "reported"]).sum())
    n_orph = int(top["orphan_curated"].eq("TRUE").sum())
    n_lf = int(top["ligand_family"].eq("orphan").sum())
    print(f"  Fig. 1c top {N_TOP}: prior report {n_rep}, weak or not assessed "
          f"{N_TOP - n_rep} ({dict(grade)}); orphan_curated {n_orph} "
          f"(ligand_family=='orphan' would give {n_lf})")
    assert n_rep == 23 and n_orph == 8, (
        f"Fig. 1c keys drifted: prior report {n_rep} (expected 23), "
        f"orphan {n_orph} (expected 8)")


if __name__ == "__main__":
    main()
