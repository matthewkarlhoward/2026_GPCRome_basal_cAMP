"""Published cell-surface expression (Tedman et al. 2026) vs our basal cAMP score.

One standalone panel. Every scored construct with a Tedman value is plotted, canonical
and alternative isoform alike; the Gs-inclusive set (gs_any: Gs as a primary OR secondary
transducer in GtoPdb, plus the documented-Gs literature set) is drawn in colour on top of
the rest. Olfactory (Golf) receptors are NOT in gs_any.

The point of the figure is that the two measurements are close to independent: a receptor
sitting at the plasma membrane is not thereby signalling, and the Gs-inclusive receptors
spread across the full range of surface expression. Any correlation that does exist is
reported on the panel so it cannot be read off by eye alone.

x is the Tedman immunostaining intensity (plasma-membrane expression, PME) on a log axis.
Note the compression above ~4000 (docs/METHODS.md, "Surface expression"): the top of this
axis is not well resolved in the source data, so do not read fine structure there.

Scoring unit is the barcode. Canonical and alternative isoforms are separate constructs
and are plotted together here without distinction: of the 901 barcodes scored in both
replicates, 859 have a Tedman PME value (42 alternative isoforms do not appear in Tedman
Supplementary Data 1).

Usage: python plot_surface_vs_basal.py
"""
from pathlib import Path

import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

import figure_style as S

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "data" / "processed"      # derived tables
FIG = ROOT / "figures"
GS_COL = S.RED
OTHER_COL = S.GREY1   # same grey as the replicate-correlation scatter (Q4)


def load():
    t = pd.read_csv(HERE / "basal_vs_tedman_pme.csv")
    a = pd.read_csv(HERE / "coupling_annotation.csv")
    d = t.merge(a[["Barcode", "gs_any"]], on="Barcode", how="left")
    d = d[d["tedman_pme"].notna() & d["score"].notna()].copy()
    d["gs_any"] = d["gs_any"].fillna(False).astype(bool)
    return d


def panel(d, ax=None):
    """40 x 40 mm. At this size the key carries no counts; the n's belong in the
    caption (they are printed by main())."""
    fig = None
    if ax is None:
        fig = S.figure(40, 40)
        ax = fig.add_axes([0.235, 0.205, 0.735, 0.755])

    sets = [(d[~d["gs_any"]], OTHER_COL, "other"), (d[d["gs_any"]], GS_COL, "Gs")]
    for sub, col, _ in sets:
        ax.scatter(sub["tedman_pme"], sub["score"], s=2.0, c=col, linewidths=0,
                   rasterized=True, zorder=3 if col is OTHER_COL else 4)

    ax.set_xscale("log")
    ax.set(xlabel="Cell-surface expression (a.u.)", ylabel="Basal cAMP score",
           ylim=(0, 1.02), yticks=[0, 0.5, 1])

    keys = [Line2D([], [], marker="o", ls="none", ms=1.8, mfc=col, mec="none",
                   label=lab) for _, col, lab in sets]
    leg = ax.legend(handles=keys, loc="upper left", handlelength=0.6, handletextpad=0.2,
                    borderpad=0.25, labelspacing=0.22, frameon=True, borderaxespad=0.15)
    leg.get_frame().set(facecolor=S.WHITE, edgecolor=S.BLACK, linewidth=0.5, alpha=1.0)
    leg.set_zorder(6)

    rho = spearmanr(d["tedman_pme"], d["score"]).correlation
    pv = spearmanr(d["tedman_pme"], d["score"]).pvalue
    g = d[d["gs_any"]]
    rho_g = spearmanr(g["tedman_pme"], g["score"]).correlation
    ax.text(0.98, 0.03, f"all $\\rho$ = {rho:.2f}\nGs $\\rho$ = {rho_g:.2f}",
            transform=ax.transAxes, ha="right", va="bottom", color=S.BLACK,
            linespacing=1.4)

    if fig is not None:
        S.save(fig, FIG / "P_surface_vs_basal")
    print(f"  n={len(d)} constructs, Gs {len(g)}, other {(~d['gs_any']).sum()}")
    print(f"  surface vs basal: rho={rho:.3f} (p={pv:.1e})")
    print(f"  within Gs: rho={rho_g:.3f} (p={spearmanr(g['tedman_pme'], g['score']).pvalue:.1e})")
    print(f"  Gs median PME {g['tedman_pme'].median():.0f} vs "
          f"{d.loc[~d['gs_any'], 'tedman_pme'].median():.0f}")


def main():
    S.apply()
    d = load()
    print("surface expression vs basal score")
    panel(d)
    print("done -> figures/P_surface_vs_basal.pdf")


if __name__ == "__main__":
    main()
