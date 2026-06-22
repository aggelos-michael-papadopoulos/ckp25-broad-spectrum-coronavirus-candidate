#!/usr/bin/env bash
# Generate the main trajectory summaries and paper-style figures.
#
# Usage:
#   bash scripts/run_all_analysis.sh [results_dir] [figure_dir] [ns]

set -euo pipefail

RESULTS_DIR="${1:-runs}"
FIGURE_DIR="${2:-paper_figures}"
NS="${3:-100}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

SYSTEMS=(
  "MERS_B_Spike_md"
  "SARS-CoV-1_E_spike_md"
  "hCoV_NL63_E_Spike_md"
  "hCoV_229E_E_spike_md"
  "hCoV_HKU1_A_spike_md"
)

mkdir -p "$ROOT/$FIGURE_DIR"

for system in "${SYSTEMS[@]}"; do
  outdir="$ROOT/$RESULTS_DIR/$system"
  if [ ! -d "$outdir" ]; then
    echo "Skipping $system: $outdir does not exist."
    continue
  fi

  echo ""
  echo "Analysing $system"
  python "$ROOT/gromacs_files/analyze_gmx.py" --outdir "$outdir" --stride 10
  python "$ROOT/gromacs_files/make_viewer_gmx.py" --outdir "$outdir"
done

echo ""
echo "Writing summary table and aggregate figures..."
python "$ROOT/extract_metrics_all.py" --results-dir "$ROOT/$RESULTS_DIR" > "$ROOT/$FIGURE_DIR/md_metrics_summary.md"
python "$ROOT/plot_all_metrics.py" --results-dir "$ROOT/$RESULTS_DIR" --out-dir "$ROOT/$FIGURE_DIR" --ns "$NS"
python "$ROOT/plot_rmsf_rg_combined.py" --results-dir "$ROOT/$RESULTS_DIR" --out-dir "$ROOT/$FIGURE_DIR" --ns "$NS"
python "$ROOT/plot_polar_stability.py" --results-dir "$ROOT/$RESULTS_DIR" --out-dir "$ROOT/$FIGURE_DIR" --ns "$NS"
python "$ROOT/make_combined_heatmap.py" --results-dir "$ROOT/$RESULTS_DIR" --data-dir "$ROOT/Data" --out-dir "$ROOT/$FIGURE_DIR"

echo ""
echo "Analysis finished. Figures and tables are in: $FIGURE_DIR/"
