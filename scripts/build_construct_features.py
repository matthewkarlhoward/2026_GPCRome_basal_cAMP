#!/usr/bin/env python
"""receptor_construct_features.csv -- one row per barcoded construct.

Joins the basal scores to the library annotation (GPCR class, TMD count, signal
peptide, ORF length) and to the coupling annotation, and assigns each class-A
receptor a ligand family from its gene symbol.

Table-only. Extracted verbatim from plot_class_features.py, whose plotting code is not
part of this release; the output is byte-identical to what that script wrote.

Usage: python scripts/build_construct_features.py
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "data" / "processed"
LIB_CSV = ROOT / "data" / "library" / "gpcrome_library.csv"

# ---- class-A ligand family from the gene symbol -----------------------------
# Ordered rules: the first matching pattern wins. Prefix-based by design (the
# library has no ligand annotation), so these are families of *names*, and a few
# receptors with unconventional symbols fall into "other".
FAMILY_RULES = [
    ("olfactory",  r"^OR\d"),
    ("gustducin",      r"^TAS\d?R"),
    ("aminergic",  r"^(HTR|ADRA|ADRB|CHRM|DRD|HRH|TAAR)"),
    ("chemokine",  r"^(CCR|CXCR|CX3CR|XCR|ACKR|CCRL)"),
    ("opsin",      r"^(OPN|RHO$|RGR$|RRH$)"),
    ("nucleotide", r"^(P2RY|ADORA|GPR17$)"),
    ("lipid",      r"^(LPAR|S1PR|CNR|PTGER|PTGDR|PTGFR|PTGIR|TBXA2R|FFAR|HCAR|CYSLTR|OXER|"
                   r"LTB4R|PTAFR$|GPBAR1$|GPR18$|GPR55$|GPR119$|GPR120$)"),
    ("peptide",    r"^(AVPR|OXTR|NPY|NPFFR|NPSR|NPBWR|GHSR|MLNR|KISS1R|TACR|OPR|SSTR|MC\dR|"
                   r"GALR|CCKAR|CCKBR|AGTR|BDKRB|EDNR|GNRHR|TRHR|PRLHR|QRFPR|HCRTR|APLNR|"
                   r"C3AR|C5AR|CMKLR|FPR|RXFP|NMUR|MCHR|PROKR|UTS2R|GRPR|NMBR|BRS3|F2R|F2RL|"
                   r"NTSR|TSHR$|LHCGR$|FSHR$|LGR\d|MLN|RLN|NPFF|XPR)"),
    ("orphan",     r"^(GPR|MAS1|TPRA1|MRGPR)"),
]


def assign_family(sym):
    for name, pat in FAMILY_RULES:
        if re.match(pat, sym):
            return name
    return "other"


def load():
    s = pd.read_csv(HERE / "basal_scores_per_barcode.csv")
    s = s.loc[s["n_days"] == 2].copy()

    lib = pd.read_csv(LIB_CSV, encoding="utf-8-sig", low_memory=False)
    lib = lib.rename(columns={"unique_molecular_identifier": "Barcode",
                              "signal_peptide_topcons": "signal_peptide",
                              "gly_fix?": "gly_fix", "len_gene_bloc_seq": "orf_len"})
    keep = ["Barcode", "gpcr_class", "num_tmds", "signal_peptide", "gly_fix", "orf_len",
            "approved_name", "num_goldengate_sites_removed"]
    df = s.merge(lib[keep], on="Barcode", how="left", validate="one_to_one")

    coup = (pd.read_csv(HERE / "coupling_annotation.csv")[
        ["Barcode", "coupling_display", "gs_primary", "gs_any"]]
            .rename(columns={"coupling_display": "coupling"}))
    df = df.merge(coup, on="Barcode", how="left")
    df["coupling"] = df["coupling"].fillna("unknown")
    df[["gs_primary", "gs_any"]] = df[["gs_primary", "gs_any"]].fillna(False)

    df["ligand_family"] = np.where(df["gpcr_class"].eq("A"),
                                   df["symbol"].map(assign_family), "n/a (non-class-A)")
    df["is_orphan"] = (df["gpcr_class"].eq("A")
                       & (df["ligand_family"].eq("orphan")
                          | df["approved_name"].fillna("").str.contains("orphan", case=False)))
    df["log_reads"] = np.log10(df["total_reads"])
    df["read_q"] = pd.qcut(df["total_reads"], 4, labels=["Q1", "Q2", "Q3", "Q4"])
    return df



def main():
    df = load()
    out_cols = ["symbol", "transcript_id", "isoform", "Barcode", "score", "score_sd",
                "total_reads", "coupling", "gpcr_class", "ligand_family", "is_orphan",
                "num_tmds", "signal_peptide", "gly_fix", "orf_len", "gs_primary", "gs_any"]
    df[out_cols].to_csv(HERE / "receptor_construct_features.csv", index=False)
    print(f"wrote receptor_construct_features.csv ({len(df)} rows)")


if __name__ == "__main__":
    main()
