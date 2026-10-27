"""Figure 1, assembled at 165 mm in three even columns.

    col 1   a  blank square, reserved for the question/assay cartoon
            d  basal score by G-protein coupling
    col 2   b  the whole library, ranked
            e  class A by ligand family
    col 3   c  the top 30, full height

Emits two versions, differing only in panel b:

    F1_main_canonical.pdf      727 canonical receptors
    F1_main_with_isoforms.pdf  all 901 scored barcodes (727 canonical + 174 alt isoforms)

Panel b carries the coupling colour key, which panel d also uses.

Usage: python plot_figure1.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy.stats import kruskal, mannwhitneyu

import figure_style as S

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "data" / "processed"      # derived tables
EXT = ROOT / "data" / "external"        # third-party reference data
INHOUSE = ROOT / "data" / "inhouse"     # other screens of the same library
FIG = ROOT / "figures"

COUP_C = {"Gs only": S.RED, "Gs + other": S.ROSE, "Gs secondary": S.PALEPURPLE,
          "Gq": S.GREEN, "Gi": S.BLUE, "multi (no Gs)": S.GREY1,
          "Golf": S.LAVENDER, "gustducin": S.PALEBLUE,
          "unassigned": S.GREY3}
COUP_ORDER = ["Gs only", "Gs + other", "Gs secondary", "Gq", "Gi", "multi (no Gs)",
              "Golf", "gustducin", "unassigned"]
AUTOCRINE = {"PTGDR", "PTGER2", "PTGER4", "PTGIR", "ADORA2A", "ADORA2B",
             "GPBAR1", "P2RY1", "P2RY2", "S1PR1", "FFAR4"}
CLASS_B = {"GIPR", "VIPR2", "PTH1R", "PTH2R", "SCTR", "GLP1R", "GHRHR", "CALCR",
           "ADCYAP1R1"}
# evidence strength is not a distinction this figure needs to draw
FILL = {"established": S.BLACK, "reported": S.BLACK}
# curated: the ligand_family column lists GPR4 as an orphan, but it is a proton-sensing
# receptor (deorphanised), so it is excluded here
ORPHAN = {"GPR6", "GPR52", "GPR161", "GPR12", "GPR62", "GPR101", "GPR26", "GPR3",
          "GPR21", "GPR61", "GPR85", "GPR173"}
N_TOP = 30


def panel(fig, ax, letter, dx=0.052, dy=0.012):
    bb = ax.get_position()
    fig.text(bb.x0 - dx, bb.y1 + dy, letter, va="bottom", ha="left")


def load():
    can = pd.read_csv(HERE / "canonical_scores_ranked.csv")
    bar = pd.read_csv(HERE / "basal_scores_per_barcode.csv")
    bar = bar[bar["n_days"] == 2].copy()
    ann = pd.read_csv(HERE / "coupling_annotation.csv")
    cmap = dict(zip(ann["Barcode"], ann["coupling_display"]))
    bar["coupling"] = bar["Barcode"].map(cmap).fillna("unknown")
    feat = pd.read_csv(HERE / "receptor_construct_features.csv")
    feat = feat[feat["isoform"] == "canonical"].drop_duplicates("symbol")
    can["ligand_family"] = can["symbol"].map(
        dict(zip(feat["symbol"], feat["ligand_family"])))
    return can, bar, feat


# ------------------------------------------------------------------ panel a --
def draw_cartoon_slot(fig, ax):
    """Reserve the slot and leave it blank; the cartoon is placed in Illustrator."""
    ax.set_axis_off()


# ------------------------------------------------------------------ panel b --
def draw_waterfall(fig, ax, d, med, with_isoforms):
    d = d.sort_values("score", ascending=False).reset_index(drop=True)
    y = np.arange(len(d))
    ax.barh(y, d["score"] - med, left=med, height=1.0, linewidth=0,
            color=d["coupling"].map(COUP_C))
    ax.axvline(med, c=S.BLACK, lw=0.5)
    ax.invert_yaxis()
    n = len(d)
    ylab = f"Receptor rank (n = {n})"
    ax.set(xlim=(0, 1.0), xticks=[0, 0.25, 0.5, 0.75, 1.0], ylim=(n, -n * 0.012),
           xlabel="Basal cAMP score", ylabel=ylab)
    # the median line is self-evident; its value belongs in the legend, not on the bars
    handles = [Line2D([], [], ls="", marker="s", ms=2.4, c=COUP_C[k], label=k)
               for k in COUP_ORDER]
    ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.12, 0.02), ncol=1,
              handlelength=0.5, handletextpad=0.3, labelspacing=0.3,
              borderpad=0.15, frameon=False)


# ------------------------------------------------------------------ panel d --
def draw_coupling(fig, ax, can, med):
    rng = np.random.default_rng(0)
    order = COUP_ORDER[::-1]
    for i, k in enumerate(order):
        v = can.loc[can["coupling"] == k, "score"].to_numpy()
        if not len(v):
            continue
        ax.scatter(v, i + rng.uniform(-0.26, 0.26, len(v)), s=1.5, c=COUP_C[k],
                   linewidths=0, rasterized=True, zorder=2)
        ax.boxplot([v], positions=[i], widths=0.62, showfliers=False, vert=False,
                   medianprops=dict(color=S.BLACK))
        ax.text(1.035, i, f"{len(v)}", va="center", ha="left", color=S.GREY1)
    ax.axvline(med, c=S.GREY2, ls="--", zorder=0)
    ax.set_yticks(range(len(order)), order)
    ax.set(xlim=(0, 1.0), xticks=[0, 0.25, 0.5, 0.75, 1], ylim=(-0.7, len(order) - 0.3),
           xlabel="Basal cAMP score")


# ------------------------------------------------------------------ panel e --
def draw_ligand(fig, ax, can, feat):
    A = can.merge(feat[["symbol", "gpcr_class"]], on="symbol", how="left")
    A = A[A["gpcr_class"] == "A"].dropna(subset=["ligand_family"])
    olf = A[A["ligand_family"] == "olfactory"]
    rest = A[A["ligand_family"] != "olfactory"]
    fams = (A.groupby("ligand_family")["score"].agg(["size", "median"])
            .query("size >= 5").sort_values("median"))
    order = list(fams.index)
    rng = np.random.default_rng(0)
    for i, k in enumerate(order):
        v = A.loc[A["ligand_family"] == k, "score"].to_numpy()
        ax.scatter(v, i + rng.uniform(-0.26, 0.26, len(v)), s=1.5, c=S.GREY2,
                   linewidths=0, rasterized=True, zorder=2)
        ax.boxplot([v], positions=[i], widths=0.62, showfliers=False, vert=False,
                   medianprops=dict(color=S.BLACK))
        ax.text(1.035, i, f"{len(v)}", va="center", ha="left", color=S.GREY1)
    ax.axvline(rest["score"].median(), c=S.GREY2, ls="--", zorder=0)
    u, pv = mannwhitneyu(rest["score"], olf["score"], alternative="greater")
    k_ = kruskal(*[v["score"].to_numpy() for _, v in rest.groupby("ligand_family")
                   if len(v) > 4])
    ax.set_yticks(range(len(order)), order)
    ax.set(xlim=(0, 1.0), xticks=[0, 0.25, 0.5, 0.75, 1], ylim=(-0.7, len(order) - 0.3),
           xlabel="Basal cAMP score")
    return (f"olfactory vs rest of class A: AUC {u / (len(olf) * len(rest)):.2f}, "
            f"p = {pv:.0e}; among the other nine families Kruskal p = {k_.pvalue:.2f}")


# ------------------------------------------------------------------ panel c --
def draw_top30(fig, ax, can, med):
    t = can.head(N_TOP).iloc[::-1].reset_index(drop=True)
    ax.axvline(med, c=S.GREY2, ls="--", zorder=0)
    for i, r in t.iterrows():
        ax.errorbar(r["score"], i, xerr=r["score_sd"], fmt="none", ecolor=S.GREY1,
                    elinewidth=0.5, zorder=2)
        ax.scatter([r["score"]], [i], s=9,
                   facecolor=FILL.get(r["ca_evidence"], S.WHITE), edgecolor=S.BLACK,
                   linewidths=0.5, zorder=3)
        if r["symbol"] in ORPHAN:
            ax.scatter([0.262], [i], s=6, marker="s", c=S.GREY1, linewidths=0,
                       clip_on=False, zorder=4)
        elif r["symbol"] in AUTOCRINE:
            ax.scatter([0.262], [i], s=7, marker="^", c=S.GREY1, linewidths=0,
                       clip_on=False, zorder=4)
    ax.set_yticks(range(len(t)), t["symbol"])
    # the axis stops at 1.0; the strip beyond the rule is annotation, not data
    ax.set(xlim=(0.3, 1.0), xticks=[0.4, 0.6, 0.8, 1.0],
           ylim=(-0.7, len(t) - 0.85), xlabel="Basal cAMP score")
    ax.tick_params(axis="y", length=0, pad=10)   # opens the lane for the flags
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

    bb = ax.get_position()

    n_rep = int(t["ca_evidence"].isin(["established", "reported"]).sum())
    evidence = [
        Line2D([], [], ls="", marker="o", ms=2.6, mfc=S.BLACK, mec=S.BLACK, mew=0.5,
               label=f"reported constitutive ({n_rep})"),
        Line2D([], [], ls="", marker="o", ms=2.6, mfc=S.WHITE, mec=S.BLACK, mew=0.5,
               label=f"not assessed ({N_TOP - n_rep})"),
    ]
    flags = [
        Line2D([], [], ls="", marker="s", ms=2.2, c=S.GREY1,
               label=f"orphan ({int(t['symbol'].isin(ORPHAN).sum())})"),
        Line2D([], [], ls="", marker="^", ms=2.4, c=S.GREY1,
               label="agonist in situ"),
    ]
    # both keys sit on the baseline of the d/e axis labels
    leg1 = ax.legend(handles=flags, loc="lower left",
                     bbox_to_anchor=(bb.x0 - 0.062, 0.050), bbox_transform=fig.transFigure,
                     handlelength=0.5, handletextpad=0.3, borderpad=0.1,
                     labelspacing=0.3, frameon=False)
    ax.add_artist(leg1)
    leg2 = ax.legend(handles=evidence, loc="lower right",
                     bbox_to_anchor=(bb.x1 + 0.045, 0.050), bbox_transform=fig.transFigure,
                     handlelength=0.5, handletextpad=0.3, borderpad=0.1, labelspacing=0.3,
                     frameon=False)
    return [(leg1, bb.x0 - 0.062), (leg2, bb.x1 + 0.045)]


def build(with_isoforms, can, bar, feat):
    S.apply()
    med = can["score"].median()
    fig = S.figure(165, 124)
    gs = GridSpec(2, 3, figure=fig, left=0.108, right=0.952, top=0.945, bottom=0.085,
                  wspace=0.95, hspace=0.42, width_ratios=[1, 1, 1],
                  height_ratios=[1, 1])

    ax_a = fig.add_subplot(gs[0, 0]); draw_cartoon_slot(fig, ax_a)
    ax_b = fig.add_subplot(gs[0, 1])
    ax_d = fig.add_subplot(gs[1, 0])
    ax_e = fig.add_subplot(gs[1, 1])
    ax_c = fig.add_subplot(gs[:, 2])
    bb = ax_c.get_position()                      # raise it to clear its key
    ax_c.set_position([bb.x0, bb.y0 + 0.080, bb.width, bb.height - 0.080])

    src = bar if with_isoforms else can
    wmed = src["score"].median()
    draw_waterfall(fig, ax_b, src, wmed, with_isoforms)
    draw_coupling(fig, ax_d, can, med)
    draw_ligand(fig, ax_e, can, feat)
    keys = draw_top30(fig, ax_c, can, med)


    # drop the keys so their baseline matches the bottom of the d/e axis labels
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    y_lab = min(inv.transform(a.xaxis.label.get_window_extent())[0][1]
                for a in (ax_d, ax_e))
    for leg, xa in keys:
        leg.set_bbox_to_anchor((xa, y_lab), transform=fig.transFigure)

    name = "F1_main_with_isoforms" if with_isoforms else "F1_main_canonical"
    S.save(fig, FIG / name)
    print(f"  {name}: waterfall n = {len(src)}, median {wmed:.3f}")


def main():
    can, bar, feat = load()
    print("Figure 1, 165 mm, three columns")
    build(False, can, bar, feat)
    build(True, can, bar, feat)
    print("done")


if __name__ == "__main__":
    main()
