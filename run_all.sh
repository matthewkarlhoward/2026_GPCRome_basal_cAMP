#!/bin/bash
# Rebuild Figure 1 and Supplemental Figure 1 from data/counts/.
# All five derived tables in data/processed/ regenerate byte-identically.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PY:-python}"
export MPLBACKEND="${MPLBACKEND:-Agg}"

# --- counts -> scores -> the tables the two figures read ---
for s in score_basal build_coupling build_ranked_tables build_construct_features \
         build_surface_join; do
  echo "=== $s"
  "$PY" "scripts/$s.py"
done

# --- the published panels ---
for s in plot_figure1 plot_qc_panels plot_surface_vs_basal; do
  echo "=== $s"
  "$PY" "scripts/$s.py"
done

echo
echo "Done. Figure 1 panels b-e: figures/F1_main_with_isoforms.pdf"
echo "      Supplemental Figure 1 panels b-f: figures/Q1-Q4, figures/P_surface_vs_basal.pdf"
