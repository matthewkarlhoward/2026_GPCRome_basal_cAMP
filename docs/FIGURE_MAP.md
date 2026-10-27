# Published panel → script

Cartoons (Figure 1a) were drawn by hand and FACS gating plots (Supplemental Figure 1a) were
exported from FlowJo. Neither is produced by this repository. Everything else below is.

Panels were emitted as vector PDFs and arranged in Illustrator; nothing was resized or
restyled after export, so the font sizes and line weights in `scripts/gpr161_style.py`
(Helvetica, 6 pt, 0.5 pt rules) are the published ones.

## Figure 1 — the basal GPCRome screen

| Panel | Script | Output |
|---|---|---|
| a | — | hand-drawn cartoon |
| b | `scripts/plot_figure1.py` | `figures/F1_main_with_isoforms.pdf`, panel b |
| c | `scripts/plot_figure1.py` | `figures/F1_main_with_isoforms.pdf`, panel c |
| d | `scripts/plot_figure1.py` | `figures/F1_main_with_isoforms.pdf`, panel d |
| e | `scripts/plot_figure1.py` | `figures/F1_main_with_isoforms.pdf`, panel e |

`plot_figure1.py` emits two versions that differ only in panel b:
`F1_main_canonical` (727 canonical receptors) and `F1_main_with_isoforms`
(all 901 scored barcodes). **The published figure is the with-isoforms version**; panels
c–e are canonical-only in both.

## Supplemental Figure 1 — sort and screen QC

| Panel | Script | Output |
|---|---|---|
| a | — | FlowJo gating export |
| b | `scripts/plot_qc_panels.py` | `figures/Q1_reads_per_sample.pdf` |
| c | `scripts/plot_qc_panels.py` | `figures/Q2_reads_per_barcode.pdf` |
| d | `scripts/plot_qc_panels.py` | `figures/Q3_score_histogram.pdf` |
| e | `scripts/plot_qc_panels.py` | `figures/Q4_replicates.pdf` |
| f | `scripts/plot_surface_vs_basal.py` | `figures/P_surface_vs_basal.pdf` |

`Q_qc_strip.pdf` is panels b–e pre-arranged as one 165 × 46 mm strip.

## The tables behind them

Five scripts build the tables the panels read. None of them plots anything.

| Script | Table it writes | Panel that needs it |
|---|---|---|
| `score_basal.py` | `basal_scores_per_barcode.csv` | Fig 1b, Supp 1c–e |
| `build_coupling.py` | `coupling_annotation.csv` | Fig 1b, 1d, Supp 1c, 1f |
| `build_ranked_tables.py` | `canonical_scores_ranked.csv` | Fig 1b–e |
| `build_construct_features.py` | `receptor_construct_features.csv` | Fig 1e |
| `build_surface_join.py` | `basal_vs_tedman_pme.csv` | Supp 1f |

The last three were extracted verbatim from the larger exploratory scripts that originally
wrote these tables (`plot_basal.py`, `plot_class_features.py`, `plot_basal_vs_tedman.py`).
Their output is byte-identical; the plotting code those scripts carried is not part of
this release.
