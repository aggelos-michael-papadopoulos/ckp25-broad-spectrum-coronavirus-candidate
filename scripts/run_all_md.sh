#!/usr/bin/env bash
# Run the five manuscript CKP-25 coronavirus spike MD systems.
#
# Usage:
#   bash scripts/run_all_md.sh [ns] [results_dir]
#
# Example:
#   bash scripts/run_all_md.sh 100 runs

set -euo pipefail

NS="${1:-100}"
RESULTS_DIR="${2:-runs}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUNNER="$ROOT/gromacs_files/run_MD_gromacs.sh"

SYSTEMS=(
  "MERS_B_Spike_md|Data/proteins/MERS_4KR0_B_spike.pdb|Data/docked_poses/MERS_CKP25_top_pose.sdf"
  "SARS-CoV-1_E_spike_md|Data/proteins/SARS1_2AJF_E_spike.pdb|Data/docked_poses/SARS1_CKP25_top_pose.sdf"
  "hCoV_NL63_E_Spike_md|Data/proteins/NL63_3KBH_E_spike.pdb|Data/docked_poses/NL63_CKP25_top_pose.sdf"
  "hCoV_229E_E_spike_md|Data/proteins/229E_6ATK_E_spike.pdb|Data/docked_poses/229E_CKP25_top_pose.sdf"
  "hCoV_HKU1_A_spike_md|Data/proteins/HKU1_8Y7Y_A_spike.pdb|Data/docked_poses/HKU1_CKP25_top_pose.sdf"
)

mkdir -p "$ROOT/$RESULTS_DIR"

for entry in "${SYSTEMS[@]}"; do
  IFS="|" read -r system protein pose <<< "$entry"
  echo ""
  echo "============================================================"
  echo "Running $system for ${NS} ns"
  echo "============================================================"

  bash "$RUNNER" \
    "$ROOT/$protein" \
    "$ROOT/$pose" \
    "$NS" \
    "$ROOT/$RESULTS_DIR/$system"
done

echo ""
echo "All MD simulations finished. Outputs are in: $RESULTS_DIR/"
