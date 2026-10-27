#!/usr/bin/env python
"""basal_vs_tedman_pme.csv -- our basal scores joined to published surface expression.

Joins the GPCRome basal Gs scores to the plasma-membrane expression measured for the
SAME library by Tedman et al. 2026 Nat Commun (DOI 10.1038/s41467-026-76564-7),
Supplementary Data 1, on Ensembl transcript ID. Canonical barcodes take the "Canonical"
sheet and alternative isoforms the "Isoforms" sheet, each alt also carrying its own
canonical's basal score.

Scoring unit is the BARCODE. Canonical and alternative isoforms are separate constructs
and are never averaged together.

Table-only. Extracted verbatim from plot_basal_vs_tedman.py, whose plotting code is not
part of this release; the output is byte-identical to what that script wrote.

Optional, and NOT part of this release:
  data/inhouse/20260422_surface_screen_weighted_mean_bin.csv   in-house April 2026 surface
      screen. That screen shares experimental lineage with this one, so it is not an
      independent surface measurement and its per-construct values are not published
      here. Without it the april_* columns are left out; everything else is unchanged.

Usage: python scripts/build_surface_join.py
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "data" / "processed"      # derived tables
EXT = REPO / "data" / "external"        # third-party reference data
INHOUSE = REPO / "data" / "inhouse"     # other screens of the same library
TED = EXT / "tedman2026"
S3 = TED / "41467_2026_76564_MOESM3_ESM.xlsx"
APRIL = INHOUSE / "20260422_surface_screen_weighted_mean_bin.csv"
HAVE_APRIL = APRIL.exists()   # see the docstring; absent in the public release
APRIL_MIN_READS = 200
# Tedman PME compresses towards its 6515 maximum; above ~4000 the reported
# per-receptor SD collapses ~10-fold, so those values are a lower bound only.
PME_COMPRESSED = 4000
COLLAPSE = {}

ENST_RE = re.compile(r"(ENST\d+)")


def enst(value) -> str | None:
    m = ENST_RE.search(str(value))
    return m.group(1) if m else None


def load():
    basal = pd.read_csv(ROOT / "basal_scores_per_barcode.csv")
    basal = basal[basal.n_days == 2].copy()

    ranked = pd.read_csv(ROOT / "canonical_scores_ranked.csv")
    cmap = dict(zip(ranked.symbol, ranked.coupling))

    def coupling_of(sym):
        v = cmap.get(sym)
        if v is None or pd.isna(v):
            return "unassigned"
        return COLLAPSE.get(v, v)

    basal["coupling"] = basal.symbol.map(coupling_of)

    can = pd.read_excel(S3, sheet_name="Canonical")
    can["tid"] = can["Receptor Name"].map(enst)
    iso = pd.read_excel(S3, sheet_name="Isoforms")
    iso["tid"] = iso["Isoform Transcript"].map(enst)
    # Two alt transcripts (GRM8 ENST00000358373, RXFP2 ENST00000380314) appear twice in
    # the Isoforms sheet as a plain and an "_SP_HA" construct with very different PME,
    # and our library likewise carries two barcodes each. The barcode-to-construct
    # pairing is not recoverable, so prefer the plain row and flag these as ambiguous.
    iso["sp_ha"] = iso["Isoform Transcript"].astype(str).str.contains("_SP_HA")
    ambiguous = set(iso.loc[iso.tid.duplicated(keep=False), "tid"])
    iso_u = iso[~iso.sp_ha].drop_duplicates("tid")

    if HAVE_APRIL:
        april = pd.read_csv(APRIL)[["Barcode", "weighted_mean_bin", "total_count"]]
        april = april.rename(columns={"weighted_mean_bin": "april_surface_wmb",
                                      "total_count": "april_reads"})
    else:
        april = pd.DataFrame({"Barcode": pd.Series(dtype=str),
                              "april_surface_wmb": pd.Series(dtype=float),
                              "april_reads": pd.Series(dtype=float)})

    # ---- canonical barcodes -> Canonical sheet
    keep = ["tid", "Immunostaining Intensity", "Immunostaining Standard Deviation",
            "GPCR Class", "Protein Length", "GC Content",
            "Library RPKM RNAseq average", "Mock RPKM RNAseq average",
            "TE", "mono_TPM", "rna_total_coverage"]
    bc = basal[basal.isoform == "canonical"].merge(
        can[keep], left_on="transcript_id", right_on="tid", how="left")

    # ---- alt barcodes -> Isoforms sheet (paired with their own canonical)
    isokeep = ["tid", "Surface Expression", "Surface Expression Standard Deviation",
               "Canonical Surface Expression", "Canonical Surface Expression Standard Deviation",
               "% change in PME", "Difference in Isoform PME",
               "length (AA)", "canonical length", "length difference"]
    ba = basal[basal.isoform == "alt_isoform"].merge(
        iso_u[isokeep], left_on="transcript_id", right_on="tid", how="left")
    ba["ambiguous_construct"] = ba.transcript_id.isin(ambiguous)
    # each alt's own canonical basal score (mean over replicate canonical barcodes
    # of the same construct; 3 symbols have two canonical barcodes)
    canon_basal = basal[basal.isoform == "canonical"].groupby("symbol").score.mean()
    ba["canonical_basal"] = ba.symbol.map(canon_basal)
    ba["d_basal"] = ba.score - ba.canonical_basal
    ba = ba.merge(can[["tid", "GPCR Class"]].rename(columns={"tid": "ctid"}),
                  left_on="symbol", right_on="ctid", how="left").drop(columns=["ctid"])

    for d in (bc, ba):
        d.drop(columns=[c for c in ("tid",) if c in d.columns], inplace=True)
    bc = bc.merge(april, on="Barcode", how="left")
    ba = ba.merge(april, on="Barcode", how="left")

    bc["pme"] = bc["Immunostaining Intensity"]
    bc["pme_sd"] = bc["Immunostaining Standard Deviation"]
    ba["pme"] = ba["Surface Expression"]
    ba["pme_sd"] = ba["Surface Expression Standard Deviation"]
    for d in (bc, ba):
        d["log10_pme"] = np.log10(d["pme"])
        d["in_april_subset"] = d["april_reads"].fillna(0) >= APRIL_MIN_READS

    bc["olfactory"] = bc["GPCR Class"].astype(str).str.contains("Olfactory")
    return bc, ba


def write_csv(bc, ba):
    cols = ["symbol", "transcript_id", "isoform", "Barcode", "score", "score_sd",
            "total_reads", "pme", "log10_pme", "pme_sd", "coupling"]
    # the in-house April columns are written only when that screen is available locally
    apr = ["april_surface_wmb", "april_reads", "in_april_subset"] if HAVE_APRIL else []
    a = bc[cols + ["GPCR Class", "olfactory"] + apr].copy()
    for c in ("canonical_pme", "canonical_pme_sd", "pct_change_pme",
              "canonical_basal", "d_basal", "ambiguous_construct"):
        a[c] = np.nan
    b = ba.copy()
    b["GPCR Class"] = np.nan
    b["olfactory"] = np.nan
    b = b.rename(columns={"Canonical Surface Expression": "canonical_pme",
                          "Canonical Surface Expression Standard Deviation":
                              "canonical_pme_sd",
                          "% change in PME": "pct_change_pme"})
    b = b[cols + ["GPCR Class", "olfactory"] + apr
          + ["canonical_pme", "canonical_pme_sd", "pct_change_pme",
             "canonical_basal", "d_basal", "ambiguous_construct"]]
    out = pd.concat([a, b], ignore_index=True)
    out = out.rename(columns={"pme": "tedman_pme", "log10_pme": "tedman_log10_pme",
                              "pme_sd": "tedman_pme_sd", "GPCR Class": "gpcr_class"})
    # Tedman PME compresses toward its 6515 ceiling; flag values that are a lower bound
    out["pme_compressed"] = out.tedman_pme >= PME_COMPRESSED
    out["canonical_pme_compressed"] = out.canonical_pme >= PME_COMPRESSED
    out = out.sort_values(["symbol", "isoform", "Barcode"])
    p = ROOT / "basal_vs_tedman_pme.csv"
    out.to_csv(p, index=False, float_format="%.6g")
    print(f"\nwrote {p}  ({len(out)} rows, "
          f"{out.tedman_pme.notna().sum()} with Tedman PME)")

def main():
    bc, ba = load()
    write_csv(bc, ba)


if __name__ == "__main__":
    main()
