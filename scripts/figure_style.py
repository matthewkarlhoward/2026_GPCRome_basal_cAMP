"""House figure style for the manuscript figures.

- Helvetica, black, 6 pt everywhere
- Figure size: width <= 165 mm, height <= 180 mm
- 0.5 pt black lines (axes, ticks, plot lines, edges)
- Fixed palette below

HELVETICA IS REQUIRED and is looked up only at the macOS system path below. Where it is
not found, matplotlib falls back to DejaVu Sans, whose metrics differ: the panels still
build but their typography is not the published typography, and roughly 5% of pixels
change. Add your own Helvetica path to _register_helvetica() to reproduce the figures
off macOS.

Usage:
    import figure_style as S          # this module sits beside the scripts that use it
    S.apply()
    fig = S.figure(width_mm=165, height_mm=120)
    ...
    S.save(fig, ROOT / "figures" / "name")   # writes .pdf and .png, checks size
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager

MM = 1 / 25.4
MAX_W_MM, MAX_H_MM = 165.0, 180.0
FONT_PT = 6
LW = 0.5

# ---- palette (sampled) -------------------------------------------------------
BLACK = "#000000"
GREY1 = "#808180"  # dark grey
GREY2 = "#979798"
GREY3 = "#CDCECE"
GREY4 = "#E2E2E2"
GREY5 = "#EBEBEB"  # lightest grey
RED = "#8F2924"
BLUE = "#4877B5"
PURPLE = "#65159E"
LAVENDER = "#9A8BD0"
GREEN = "#316B3E"
ROSE = "#DBB8AB"        # pastel of RED
PALEBLUE = "#CCD5EF"    # pastel of BLUE
PALEPURPLE = "#F0EAF9"  # pastel of PURPLE
PALELAV = "#F5F3FB"     # pastel of LAVENDER
WHITE = "#FFFFFF"       # use with a black 0.5 pt edge

GREYS = [BLACK, GREY1, GREY2, GREY3, GREY4, GREY5]
CATEGORICAL = [BLACK, RED, BLUE, PURPLE, GREEN, LAVENDER, GREY2]


def _register_helvetica():
    for p in ("/System/Library/Fonts/Helvetica.ttc",):
        if Path(p).exists():
            try:
                font_manager.fontManager.addfont(p)
            except Exception:
                pass


def apply():
    _register_helvetica()
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": FONT_PT,
        "axes.titlesize": FONT_PT,
        "axes.labelsize": FONT_PT,
        "xtick.labelsize": FONT_PT,
        "ytick.labelsize": FONT_PT,
        "legend.fontsize": FONT_PT,
        "legend.title_fontsize": FONT_PT,
        "figure.titlesize": FONT_PT,
        "axes.titleweight": "normal",
        "text.color": BLACK,
        "axes.labelcolor": BLACK,
        "axes.edgecolor": BLACK,
        "xtick.color": BLACK,
        "ytick.color": BLACK,
        "axes.linewidth": LW,
        "lines.linewidth": LW,
        "patch.linewidth": LW,
        "hatch.linewidth": LW,
        "xtick.major.width": LW,
        "ytick.major.width": LW,
        "xtick.minor.width": LW,
        "ytick.minor.width": LW,
        "xtick.major.size": 2,
        "ytick.major.size": 2,
        "xtick.minor.size": 1,
        "ytick.minor.size": 1,
        "grid.linewidth": LW,
        "grid.color": GREY4,
        "lines.markeredgewidth": LW,
        "scatter.edgecolors": "none",
        "errorbar.capsize": 0,
        "boxplot.boxprops.linewidth": LW,
        "boxplot.whiskerprops.linewidth": LW,
        "boxplot.capprops.linewidth": LW,
        "boxplot.medianprops.linewidth": LW,
        "boxplot.medianprops.color": BLACK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.prop_cycle": mpl.cycler(color=CATEGORICAL),
        "legend.frameon": False,
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "pdf.fonttype": 42,   # editable text in Illustrator
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def figure(width_mm=MAX_W_MM, height_mm=120, **kw):
    if width_mm > MAX_W_MM + 1e-6 or height_mm > MAX_H_MM + 1e-6:
        raise ValueError(f"figure {width_mm}x{height_mm} mm exceeds {MAX_W_MM}x{MAX_H_MM} mm")
    return plt.figure(figsize=(width_mm * MM, height_mm * MM), **kw)


def save(fig, stem, png_dpi=600, svg=False):
    """Save <stem>.pdf and <stem>.png at the figure's physical size (no bbox_inches='tight',
    so the page size is exactly the declared size). Validates size limits."""
    w_in, h_in = fig.get_size_inches()
    if w_in / MM > MAX_W_MM + 0.5 or h_in / MM > MAX_H_MM + 0.5:
        raise ValueError(f"figure is {w_in/MM:.1f}x{h_in/MM:.1f} mm; limit {MAX_W_MM}x{MAX_H_MM} mm")
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=png_dpi)
    if svg:
        fig.savefig(stem.with_suffix(".svg"))
    return stem
