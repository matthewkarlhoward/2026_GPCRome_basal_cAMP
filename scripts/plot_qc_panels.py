"""QC panels for the GPCRome basal screen, each emitted as its own PDF.

  Q1_reads_per_sample   library-barcode reads per sorted bin, log axis  (Supp Fig 1b)
  Q2_reads_per_barcode  coverage across the library, split by olfactory   (Supp Fig 1c)
  Q3_score_histogram    distribution of the basal cAMP score              (Supp Fig 1d)
  Q4_replicates         replicate 1 vs replicate 2                        (Supp Fig 1e)

Reads: data/processed/basal_scores_per_barcode.csv, coupling_annotation.csv,
       data/qc/read_fate.csv

The two replicate sorts, 09/24 and 09/26, are "replicate 1" and "replicate 2" below.

Q2 shows that olfactory receptors are not evenly spread through the coverage
distribution. They are 40% of the scored library, with median coverage 720 reads against
22,606 for everything else. That is the origin of the library's bimodality.

Usage: python plot_qc_panels.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import pearsonr, spearmanr

import figure_style as S

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "data" / "processed"      # derived tables
FIG = ROOT / "figures"
RUN = HERE.parent
DAY_A, DAY_B = "20260924", "20260926"
# s5-s8 are the 20260924 bins a-d, s9-s12 the 20260926 bins a-d
SORTED = [(f"s{i}_basal_S{i}", lab) for i, lab in
          zip(range(5, 13), ["09/24 a", "09/24 b", "09/24 c", "09/24 d",
                             "09/26 a", "09/26 b", "09/26 c", "09/26 d"])]


def panel(w=50, h=44, left=0.26, bottom=0.235, width=0.70, height=0.70):
    fig = S.figure(w, h)
    return fig, fig.add_axes([left, bottom, width, height])


def _save(fig, name):
    if fig is not None:
        S.save(fig, FIG / name)


def load():
    b = pd.read_csv(HERE / "basal_scores_per_barcode.csv")
    b = b[b["n_days"] == 2].copy()
    ann = pd.read_csv(HERE / "coupling_annotation.csv")
    b["olf"] = b["Barcode"].map(dict(zip(ann["Barcode"], ann["coupling_display"]))).eq("Golf")
    return b


# ----------------------------------------------------------------- Q1 --------
def q1_reads_per_sample(ax=None):
    fate = pd.read_csv(RUN / "qc" / "read_fate.csv").set_index("sample")
    v = [fate.loc[s, "library"] for s, _ in SORTED]
    fig = None
    if ax is None:
        fig, ax = panel(54, 44, left=0.215, bottom=0.30, width=0.74, height=0.62)
    # a gap between bins 4 and 5 keeps the two replicates legible without labelling them
    x = np.array([0, 1, 2, 3, 4.7, 5.7, 6.7, 7.7])
    ax.bar(x, v, color=S.BLACK, width=0.68, linewidth=0, zorder=3)
    ax.set_yscale("log")
    ax.set_ylim(1e5, 1e7)
    # bin letters, with the two sort days named underneath instead of rotated labels
    ax.set_xticks(x, [lab.split()[1] for _, lab in SORTED])
    for x0, lab in ((1.5, "rep 1"), (6.2, "rep 2")):
        ax.text(x0, -0.175, lab, transform=ax.get_xaxis_transform(), ha="center",
                va="top")
    ax.set(ylabel="library-barcode reads")
    ax.set_title("Sequencing depth")
    ax.set_axisbelow(True)
    ax.grid(axis="y", color=S.GREY5, lw=0.4)
    _save(fig, "Q1_reads_per_sample")
    print(f"  Q1: {min(v) / 1e6:.2f}-{max(v) / 1e6:.2f} M library reads per bin")


# ----------------------------------------------------------------- Q2 --------
def q2_reads_per_barcode(b, ax=None):
    fig = None
    if ax is None:
        fig, ax = panel(54, 46, left=0.215, bottom=0.245, width=0.74, height=0.645)
    r = b["total_reads"]
    edges = np.logspace(np.log10(max(r.min(), 1)), np.log10(r.max()), 34)
    series = [(r, S.GREY2, 1.1, f"all ({len(r)})"),
              (r[~b["olf"]], S.BLACK, 0.7, f"other ({(~b['olf']).sum()})"),
              (r[b["olf"]], S.PURPLE, 0.7, f"olfactory ({b['olf'].sum()})")]
    for v, col, lw, _ in series:
        ax.hist(v, bins=edges, histtype="step", color=col, linewidth=lw)
    # hist hands back Patch handles; these are line plots, so key them with lines
    keys = [Line2D([], [], color=col, lw=lw, label=lab) for _, col, lw, lab in series]
    ax.set_xscale("log")
    ax.set(xlabel="reads per barcode (4 bins)", ylabel="barcodes")
    ax.set_title("Library coverage")
    ax.legend(handles=keys, loc="upper left", handlelength=0.9, handletextpad=0.35,
              borderpad=0.15, labelspacing=0.32, frameon=False)
    _save(fig, "Q2_reads_per_barcode")
    print(f"  Q2: olfactory median {r[b['olf']].median():.0f} vs "
          f"{r[~b['olf']].median():.0f} reads")


# ----------------------------------------------------------------- Q3 --------
def q3_score_histogram(b, ax=None):
    fig = None
    if ax is None:
        fig, ax = panel(50, 44)
    ax.hist(b["score"], bins=np.linspace(0, 1, 41), color=S.GREY2, linewidth=0)
    ax.axvline(b["score"].median(), c=S.BLACK, ls="--")
    ax.text(b["score"].median() + 0.02, 0.97, f"median {b['score'].median():.2f}",
            transform=ax.get_xaxis_transform(), va="top", ha="left")
    ax.set(xlabel="Basal cAMP score", ylabel="barcodes", xlim=(0, 1),
           xticks=[0, 0.25, 0.5, 0.75, 1])
    ax.set_title("Score distribution")
    _save(fig, "Q3_score_histogram")
    print(f"  Q3: n={len(b)}, median {b['score'].median():.3f}")


# ----------------------------------------------------------------- Q4 --------
def q4_replicates(b, ax=None):
    fig = None
    if ax is None:
        fig, ax = panel(50, 50, left=0.25, bottom=0.21, width=0.71, height=0.735)
    x, y = b[f"{DAY_A}_wmb"], b[f"{DAY_B}_wmb"]
    ok = x.notna() & y.notna()
    x, y = x[ok], y[ok]
    ax.plot([0, 1], [0, 1], c=S.GREY2, lw=0.5, zorder=1)
    ax.scatter(x, y, s=1.8, c=S.GREY1, linewidths=0, rasterized=True, zorder=2)
    r = pearsonr(x, y)[0]
    rho = spearmanr(x, y).correlation
    ax.text(0.03, 0.97, f"r = {r:.2f}\n$\\rho$ = {rho:.2f}\nn = {len(x)}",
            transform=ax.transAxes, va="top")
    ax.set(xlim=(0, 1.02), ylim=(0, 1.02), xticks=[0, 0.5, 1], yticks=[0, 0.5, 1],
           xlabel="Basal cAMP score, rep 1",
           ylabel="Basal cAMP score, rep 2")
    ax.set_title("Replicate agreement")
    ax.set_aspect("equal")
    _save(fig, "Q4_replicates")
    print(f"  Q4: r={r:.3f}, rho={rho:.3f}, n={len(x)}")


def main():
    S.apply()
    b = load()
    print("QC panels")
    q1_reads_per_sample()
    q2_reads_per_barcode(b)
    q3_score_histogram(b)
    q4_replicates(b)
    print("done -> figures/Q*.pdf")


if __name__ == "__main__":
    main()
