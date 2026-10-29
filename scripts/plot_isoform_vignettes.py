"""Supplemental Figure 2 -- isoform vignettes: GPR161 and ADORA2A.

  a  every alternative isoform against its canonical, the two vignettes marked
  b  GPR161 isoform architecture      c  GPR161 basal vs surface
  d  ADORA2A isoform architecture     e  ADORA2A basal vs surface

Architecture is derived from the library gene blocks: the ORF is translated from
the ATG-HA start, and each isoform is aligned to its canonical by longest shared
N-terminal and C-terminal run. Lengths are reported natively, i.e. with the
10-residue Met + HA tag removed.

Inputs:  data/library/gpcrome_library.csv
         data/processed/alt_vs_canonical_isoform.csv (scripts/build_isoform_table.py)
         data/processed/basal_vs_tedman_pme.csv (scripts/build_surface_join.py)
Output:  figures/S2_isoform_vignettes.{pdf,png}

Usage: python plot_isoform_vignettes.py
"""
import sys
from pathlib import Path
from itertools import product
from textwrap import wrap

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_style as S  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
LIB = ROOT / "data" / "library" / "gpcrome_library.csv"
OUT = ROOT / "figures" / "S2_isoform_vignettes"

TAG = "ATGTACCCATACGACGTACCAGATTACGCT"   # Met + HA, 10 residues
TAGLEN = 10
VIG = {"GPR161": S.RED, "ADORA2A": S.BLUE}

_bases = "TCAG"
_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {"".join(c): _aas[i] for i, c in enumerate(product(_bases, repeat=3))}


def translate(seq):
    """Protein from the ATG-HA start, stopping at the first stop codon."""
    i = seq.find(TAG)
    if i < 0:
        return None
    p = []
    for c in wrap(seq[i:], 3):
        if len(c) < 3:
            break
        a = CODON.get(c, "X")
        if a == "*":
            break
        p.append(a)
    return "".join(p)


def architecture(sym, lib):
    """Native lengths and the shared run shared with the canonical."""
    s = lib[lib.hgnc_symbol == sym]
    can = s[s["alt_isoform?"] == "canonical"].iloc[0]
    cp = translate(can.gene_block_sequence)
    rows = [dict(tid=can.ensembl_transcript_id, canonical=True,
                 length=len(cp) - TAGLEN, n_shared=len(cp) - TAGLEN,
                 c_shared=len(cp) - TAGLEN)]
    for r in s[s["alt_isoform?"] != "canonical"].itertuples():
        ap = translate(r.gene_block_sequence)
        n = 0
        while n < min(len(ap), len(cp)) and ap[n] == cp[n]:
            n += 1
        c = 0
        while c < min(len(ap), len(cp)) - n and ap[-1 - c] == cp[-1 - c]:
            c += 1
        rows.append(dict(tid=r.ensembl_transcript_id, canonical=False,
                         length=len(ap) - TAGLEN, n_shared=max(n - TAGLEN, 0),
                         c_shared=c))
    return len(cp) - TAGLEN, pd.DataFrame(rows)


# ------------------------------------------------------------------ panel a --
def draw_scatter(ax, alt):
    rest = alt[~alt.symbol.isin(VIG)]
    ax.plot([0, 1], [0, 1], color=S.GREY3, lw=S.LW, zorder=0)
    ax.scatter(rest.canonical, rest.score, s=5, c=S.GREY2, linewidths=0,
               rasterized=True, zorder=1)
    for sym, col in VIG.items():
        d = alt[alt.symbol == sym]
        ax.scatter(d.canonical, d.score, s=12, facecolor=col, edgecolor=S.WHITE,
                   linewidths=0.35, zorder=3)
        # both labels point inward: nothing may cross the right spine
        off, ha = ((-4, 7), "right") if sym == "ADORA2A" else ((-4, -9), "right")
        ax.annotate(sym, (d.canonical.iloc[0], d.score.max()),
                    xytext=off, textcoords="offset points", ha=ha, color=col)
    ax.set(xlim=(0, 1), ylim=(0, 1), xticks=[0, 0.25, 0.5, 0.75, 1],
           yticks=[0, 0.25, 0.5, 0.75, 1],
           xlabel="canonical basal cAMP score",
           ylabel="alternative isoform score")
    ax.text(0.03, 0.97, f"n = {len(alt)} isoforms", transform=ax.transAxes,
            va="top", ha="left", color=S.GREY1)


# ------------------------------------------------------------ panels b and d --
def draw_architecture(ax, sym, can_len, arch, scores, col):
    """Isoforms as bars in canonical residue coordinates, aligned on whichever
    terminus they share. Transcript on the y axis, residue change at the right."""
    arch = arch.copy()
    arch["score"] = arch.tid.map(scores)
    arch["start"] = [0 if r.n_shared >= r.c_shared else can_len - r.length
                     for r in arch.itertuples()]
    alt = arch[~arch.canonical].sort_values("score", ascending=False)
    order = pd.concat([arch[arch.canonical], alt]).reset_index(drop=True)

    bar_max = int(max(can_len, (order.start + order.length).max()))
    pad = int(bar_max * 0.42)
    n = len(order)
    labs = []
    for i, r in order.iterrows():
        y = n - 1 - i
        if r.canonical:
            ax.barh(y, can_len, left=0, height=0.6, color=col, linewidth=0, zorder=3)
        else:
            ax.barh(y, r.length, left=r.start, height=0.6, color=S.GREY4,
                    linewidth=0, zorder=2)
            w, left = ((r.n_shared, 0) if r.n_shared >= r.c_shared
                       else (r.c_shared, can_len - r.c_shared))
            ax.barh(y, w, left=left, height=0.6, color=col, linewidth=0, zorder=3)
        labs.append("canonical" if r.canonical else "ENST\u2026" + r.tid[-5:])
        ax.text(bar_max + pad, y, f"{r.length:d} aa" if r.canonical
                else f"{r.length - can_len:+d} aa",
                va="center", ha="right", color=S.GREY1)

    ax.axvline(0, color=S.GREY3, lw=S.LW, zorder=1)
    ax.set_yticks(range(n), labs[::-1])
    ax.set(ylim=(-0.7, n - 0.3), xlim=(-25, bar_max + pad + 20),
           xticks=[0, 250, 500] if bar_max > 450 else [0, 200, 400],
           xlabel="residue")
    ax.tick_params(axis="y", length=0)
    ax.text(0.0, 1.06, sym, transform=ax.transAxes, ha="left", va="bottom",
            color=col)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)


# ------------------------------------------------------------ panels c and e --
def draw_surface(ax, sym, d, col):
    can = d[d.isoform == "canonical"].iloc[0]
    alt = d[d.isoform != "canonical"]
    ax.scatter(alt.tedman_pme, alt.score, s=11, facecolor=S.WHITE,
               edgecolor=col, linewidths=0.6, zorder=3)
    ax.scatter([can.tedman_pme], [can.score], s=16, marker="D", facecolor=col,
               edgecolor=S.BLACK, linewidths=0.5, zorder=4)
    ax.set_xscale("log")
    ax.set(xlim=(80, 1e4), ylim=(0, 1), yticks=[0, 0.5, 1],
           xlabel="surface expression (a.u.)", ylabel="basal cAMP score")


def ax_mm(fig, x, y, w, h, W, H):
    """Axes placed by millimetre on a W x H mm figure."""
    return fig.add_axes([x / W, y / H, w / W, h / H])


def main():
    S.apply()
    lib = pd.read_csv(LIB, encoding="utf-8-sig", low_memory=False)
    alt = pd.read_csv(PROC / "alt_vs_canonical_isoform.csv")
    ted = pd.read_csv(PROC / "basal_vs_tedman_pme.csv")

    W, H = 165.0, 88.0
    fig = S.figure(W, H)

    SC = 68.0                                # scatter, mm square
    GW_A, GW_S, GH = 28.0, 25.0, 26.0        # grid cells
    X_ARCH, X_SURF = 96.0, 137.0
    Y_TOP, Y_BOT = 55.0, 13.0

    axa = ax_mm(fig, 12, 13, SC, SC, W, H)
    draw_scatter(axa, alt)

    for i, sym in enumerate(VIG):
        col, y = VIG[sym], (Y_TOP if i == 0 else Y_BOT)
        can_len, arch = architecture(sym, lib)
        d = ted[ted.symbol == sym]

        axd = ax_mm(fig, X_ARCH, y, GW_A, GH, W, H)
        draw_architecture(axd, sym, can_len, arch,
                          dict(zip(d.transcript_id, d.score)), col)

        axs = ax_mm(fig, X_SURF, y, GW_S, GH, W, H)
        draw_surface(axs, sym, d, col)
        axs.text(0.0, 1.06, sym, transform=axs.transAxes, ha="left",
                 va="bottom", color=col)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    S.save(fig, OUT)
    print(f"wrote {OUT}.pdf / .png   ({W:.0f} x {H:.0f} mm)")
    print(f"  scatter {SC:.0f} x {SC:.0f} mm; grid cells "
          f"{GW_A:.0f}/{GW_S:.0f} wide x {GH:.0f} mm; "
          f"right edge {X_SURF + GW_S:.0f} mm")


if __name__ == "__main__":
    main()
