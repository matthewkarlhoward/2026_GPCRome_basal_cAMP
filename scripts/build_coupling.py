"""Rebuild the G-protein coupling annotation from a documented, citable source.

WHY NOT THE LEGACY ANNOTATION. An earlier hand-written coupling table (shipped as
data/inhouse/legacy_coupling.csv and still used here as a fallback)
covers 420 gene symbols with one primary transducer each. It agrees with IUPHAR curation
for 166/177 (94%) of the receptors both cover, but it is incomplete, and its gaps were
invisible because every uncovered receptor landed in a bucket called "unknown" rather than
being marked as unannotated. GPR4 is the clean example: a proton-sensing receptor that
IUPHAR documents as stimulating adenylate cyclase, and that scores 0.78 here, was "unknown"
only because it was absent from that table.

Note that Tedman et al. 2026 cannot supply this: their supplementary tables carry structural
class (A/B/C/F, olfactory, taste, vomeronasal) and no transducer assignment at all. Their
paper is about surface expression, not signalling. The only field we take from them is class.

SOURCE. GproteinDb's bulk coupling table (gproteindb/, see SOURCE.txt), sheet "GtoPdb" --
the IUPHAR/Guide to Pharmacology curation, 400 human receptors, keyed by gene symbol. The
other sheets in that workbook are individual screens (Bouvier GEMTA n=100, Inoue TGFa n=62,
...) keyed by UniProt short name; they are small subsets of GtoPdb's coverage, so they are
not used. If they are ever wanted, they carry a ranked "Primary (1') fam" column.

PROMISCUITY, AND WHY WE DO NOT COLLAPSE IT. 40 of the 242 curated receptors here (17%) list
more than one primary transducer family, and IUPHAR gives no rank among them, so there is no
non-arbitrary way to pick one. We therefore do not pick one: coupling is stored as four
independent booleans per tier (couples_gs/gi/gq/g12, and the same for secondary), which is
what the source actually asserts.

Keeping the tiers apart turns out to matter, because basal cAMP tracks how DEDICATED a
receptor is to Gs rather than whether Gs appears anywhere in its profile. Canonical medians:

    Gs sole primary      0.750  (n=54)
    Gs one of several    0.503  (n=19)
    Gs secondary only    0.475  (n=18)
    no documented Gs     0.34-0.46

and the separation is strongest under the strictest definition, weakening as it is loosened
(AUCs below are against non-sensory receptors only; see the comment at the AUC loop):

    Gs sole primary vs rest      AUC 0.853
    Gs among primaries vs rest   AUC 0.786
    Gs primary or secondary      AUC 0.750

These are recomputed by the script on every run -- trust its stdout over this docstring.

That ordering is the expected behaviour of a Gs/cAMP readout and is itself a result. It also
means a single "Gs vs not" split understates the signal. `coupling_display` therefore keeps
"Gs only", "Gs + other" and "Gs secondary" as separate categories; `gs_primary` and `gs_any`
are provided for whichever binary a given analysis needs.

LAYERS, most specific first, each recorded in `coupling_source`:
  1. GtoPdb primary transducer (IUPHAR curation)
  2. GtoPdb secondary transducer, where no primary is curated
  3. structural class from Tedman et al. 2026: olfactory -> Golf (Gs-family),
     taste -> gustducin (Gi-family), vomeronasal -> V1R. These are not in IUPHAR's table.
  4. receptors with documented Gs coupling that IUPHAR's table does not cover at all
     (GS_DOCUMENTED below). GPR161 is here -- it is the
     subject of this paper and is absent from the IUPHAR transducer table.
  5. the legacy curated dict, flagged as unverified, where 1-4 are all silent
  6. "unknown", which now means no documented coupling rather than "absent from a dict"

Writes coupling_annotation.csv. Usage: python build_coupling.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "data" / "processed"      # derived tables
EXT = ROOT / "data" / "external"        # third-party reference data
INHOUSE = ROOT / "data" / "inhouse"     # other screens of the same library
FIG = ROOT / "figures"
LEGACY = INHOUSE / "legacy_coupling.csv"
GTOPDB = EXT / "gproteindb" / "GPCR-G_protein_couplings.xlsx"
TEDMAN = EXT / "tedman2026" / "41467_2026_76564_MOESM3_ESM.xlsx"

# IUPHAR family tokens -> our short names. The two non-family tokens are kept out: a
# receptor whose only primary is "G protein (identity unknown)" couples something, but we
# cannot say what, so it is reported as coupling-unspecified rather than as a family.
FAM = {"Gs": "gs", "Gi/o": "gi", "Gq/11": "gq", "G12/13": "g12"}
VAGUE = {"G protein (identity unknown)", "G protein independent mechanism", "GPa1"}

# Gs coupling documented in the literature for receptors IUPHAR's transducer table omits
# entirely. Applied
# only where layers 1-3 are silent, so IUPHAR always wins where it speaks.
GS_DOCUMENTED = {"GPR161", "GPBAR1", "GHRHR", "LGR6", "LGR5", "ADGRD1", "ADGRG5"}

CLASS_RULE = {"A (Olfactory)": ("Golf", "Gs-family", "Golf"),
              "Taste": ("gustducin", "Gi-family", "gustducin"),
              "Vomeronasal": ("V1R", "unknown", "unassigned")}


def parse_transducers(s):
    """Split an IUPHAR transducer cell into our family short names, plus a vague flag."""
    if pd.isna(s):
        return set(), False
    toks = [t.strip() for t in str(s).split(",")]
    return {FAM[t] for t in toks if t in FAM}, any(t in VAGUE for t in toks)


def load_gtopdb():
    """Gene symbol -> (primary set, secondary set, vague flag, raw strings).

    Matches on the Gene column, then on each alias in GeneAlt for symbols not already taken,
    so that receptors our library names by an alias still resolve.
    """
    g = pd.read_excel(GTOPDB, "GtoPdb", header=0)
    g.columns = [str(c).replace("\n", " ").strip() for c in g.columns]
    out = {}
    for r in g.itertuples(index=False):
        rec = dict(zip(g.columns, r))
        prim, vague = parse_transducers(rec["Primary Transducer"])
        sec, _ = parse_transducers(rec["Secondary Transducer"])
        entry = (prim, sec, vague, rec["Primary Transducer"], rec["Secondary Transducer"])
        sym = str(rec["Gene"]).upper().strip()
        if sym and sym != "NAN":
            out[sym] = entry
    for r in g.itertuples(index=False):           # aliases only fill gaps
        rec = dict(zip(g.columns, r))
        if pd.isna(rec["GeneAlt"]):
            continue
        prim, vague = parse_transducers(rec["Primary Transducer"])
        sec, _ = parse_transducers(rec["Secondary Transducer"])
        entry = (prim, sec, vague, rec["Primary Transducer"], rec["Secondary Transducer"])
        for a in str(rec["GeneAlt"]).split(";"):
            a = a.upper().strip()
            if a and a not in out:
                out[a] = entry
    return out


def display_category(prim, sec, source, class_disp):
    """One exclusive label per receptor, for panels that can show only one colour.

    Ordered by how dedicated the receptor is to Gs, which is the axis the basal cAMP signal
    actually follows: sole primary > one of several primaries > secondary only > no Gs.
    """
    if class_disp:
        return class_disp
    if not prim:
        # "unassigned" covers every receptor we cannot put in a transducer family: no
        # documented coupling, a receptor-class fallback (Frizzled, adhesion, ACKR), or
        # IUPHAR recording only "G protein (identity unknown)". It is not a coupling class.
        return "unassigned"
    if "gs" in prim:
        return "Gs only" if len(prim) == 1 else "Gs + other"
    if "gs" in sec:
        return "Gs secondary"
    if len(prim) > 1:
        return "multi (no Gs)"
    # one receptor is G12/13-only; it sits with the other non-Gs multi-family receptors
    return {"gq": "Gq", "gi": "Gi", "g12": "multi (no Gs)"}[next(iter(prim))]


def main():
    d = pd.read_csv(HERE / "basal_scores_per_barcode.csv")
    d = d[d["n_days"] == 2][["symbol", "transcript_id", "isoform", "Barcode", "score"]]

    legacy = pd.read_csv(LEGACY)[["symbol", "coupling"]].rename(columns={"coupling": "legacy"})
    d = d.merge(legacy, on="symbol", how="left")

    tx = pd.ExcelFile(TEDMAN).parse("Canonical")
    tx["tid"] = tx["Receptor Name"].str.extract(r"(ENST\d+)")
    d["tedman_class"] = d["transcript_id"].map(tx.set_index("tid")["GPCR Class"].to_dict())

    gt = load_gtopdb()
    rows = []
    for r in d.itertuples():
        sym, lg, tc = r.symbol, r.legacy, str(r.tedman_class)
        prim = sec = set()
        vague = False
        raw_p = raw_s = None
        class_disp = None

        hit = gt.get(str(sym).upper())
        if hit:
            prim, sec, vague, raw_p, raw_s = hit

        if prim:
            source = "GtoPdb primary transducer (IUPHAR)"
        elif sec:
            source = "GtoPdb secondary transducer only (IUPHAR)"
            prim = sec                      # treat as the best available statement
        elif lg == "Olfactory" or tc in CLASS_RULE:
            key = "A (Olfactory)" if lg == "Olfactory" else tc
            _, _, class_disp = CLASS_RULE[key]
            source = "structural class (Tedman et al. 2026)"
        elif sym in GS_DOCUMENTED:
            prim = {"gs"}
            source = "documented in literature, absent from IUPHAR table"
        elif lg in ("Gs", "Gi", "Gq", "G12/13"):
            prim = {{"Gs": "gs", "Gi": "gi", "Gq": "gq", "G12/13": "g12"}[lg]}
            source = "legacy curated dict (unverified)"
        elif vague:
            source = "GtoPdb: couples a G protein, family not specified"
        elif lg in ("Frizzled", "Adhesion", "Atypical", "Orphan"):
            source = "receptor class, no transducer documented"
        else:
            source = "unresolved"

        rows.append({
            "couples_gs": "gs" in prim, "couples_gi": "gi" in prim,
            "couples_gq": "gq" in prim, "couples_g12": "g12" in prim,
            "gs_secondary": "gs" in sec,
            "n_primary_families": len(prim),
            "gtopdb_primary": raw_p, "gtopdb_secondary": raw_s,
            "coupling_source": source,
            "coupling_display": display_category(prim, sec, source, class_disp),
        })
    d = pd.concat([d.reset_index(drop=True), pd.DataFrame(rows)], axis=1)

    # gs_family keeps its old meaning: able to signal through a Gs-family protein at all,
    # which includes the olfactory receptors' Golf.
    d["gs_primary"] = d["couples_gs"]
    d["gs_any"] = d["couples_gs"] | d["gs_secondary"]
    d["gs_family"] = d["couples_gs"] | d["coupling_display"].eq("Golf")
    # back-compat columns for scripts that predate the multi-label scheme
    d["transducer"] = d["coupling_display"]
    d["coupling_family"] = d["coupling_display"].map(
        {"Gs only": "Gs-family", "Gs + other": "Gs-family", "Golf": "Gs-family",
         "Gq": "Gq-family", "Gi": "Gi-family", "gustducin": "Gi-family"}).fillna(
        d["coupling_display"].map({"unassigned": "unknown"})).fillna("other/mixed")

    can = d[d["isoform"] == "canonical"]
    print(f"{len(d)} barcodes, {len(can)} canonical\n")
    print("provenance:")
    print(d["coupling_source"].value_counts().to_string())
    print(f"\nstill unresolved: {(d['coupling_source'] == 'unresolved').sum()} barcodes "
          f"({(can['coupling_source'] == 'unresolved').sum()} canonical)")
    print("\ncoupling_display (what the figures show):")
    print(d["coupling_display"].value_counts().to_string())
    print("\npromiscuity among curated receptors:")
    cur = can[can["n_primary_families"] > 0]
    print(f"  {(cur['n_primary_families'] > 1).sum()}/{len(cur)} canonical "
          f"({100 * (cur['n_primary_families'] > 1).mean():.0f}%) couple >1 family")
    print("\nmedian score by display category (canonical):")
    print(can.groupby("coupling_display")["score"].agg(["size", "median"])
          .sort_values("median", ascending=False).round(3).to_string())
    from scipy.stats import mannwhitneyu
    print("\nhow the separation depends on how strictly Gs coupling is defined:")
    # Sensory receptors are excluded from the comparison set: they are 54% of the
    # canonical library and cannot signal through a Gs/cAMP readout, so leaving them in
    # inflates every AUC. These labels MUST match display_category()'s output.
    nonsens = can[~can["coupling_display"].isin(["Golf", "gustducin"])]
    for lab, mask in [("Gs sole primary", nonsens["coupling_display"].eq("Gs only")),
                      ("Gs among primaries", nonsens["gs_primary"]),
                      ("Gs primary or secondary", nonsens["gs_any"])]:
        a, b = nonsens.loc[mask, "score"], nonsens.loc[~mask, "score"]
        u, p = mannwhitneyu(a, b, alternative="greater")
        print(f"  {lab:24s} n={len(a):3d}  AUC {u / (len(a) * len(b)):.3f}  p={p:.2g}")

    d[["symbol", "transcript_id", "isoform", "Barcode", "score", "legacy", "tedman_class",
       "gtopdb_primary", "gtopdb_secondary", "couples_gs", "couples_gi", "couples_gq",
       "couples_g12", "gs_secondary", "n_primary_families", "gs_primary", "gs_any",
       "gs_family", "transducer", "coupling_family", "coupling_display",
       "coupling_source"]].to_csv(
        HERE / "coupling_annotation.csv", index=False, float_format="%.4g")


if __name__ == "__main__":
    main()
