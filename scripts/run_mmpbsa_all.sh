#!/usr/bin/env bash
# Run AmberTools MMPBSA.py for the five manuscript systems.
#
# Usage:
#   bash scripts/run_mmpbsa_all.sh [results_dir]

set -euo pipefail

RESULTS_DIR="${1:-runs}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

SYSTEMS=(
  "MERS_B_Spike_md"
  "SARS-CoV-1_E_spike_md"
  "hCoV_229E_E_spike_md"
  "hCoV_HKU1_A_spike_md"
  "hCoV_NL63_E_Spike_md"
)

for system in "${SYSTEMS[@]}"; do
  sys_dir="$ROOT/$RESULTS_DIR/$system"
  echo ""
  echo "Running MM-GBSA for $system"

  if [ ! -f "$sys_dir/complex.nc" ] || [ ! -f "$sys_dir/complex.prmtop" ]; then
    echo "Missing complex.nc or complex.prmtop in $sys_dir; skipping."
    continue
  fi

  (
    cd "$sys_dir"
    MMPBSA.py -O \
      -i "$ROOT/mmpbsa.in" \
      -cp complex.prmtop \
      -rp receptor.prmtop \
      -lp ligand.prmtop \
      -y complex.nc \
      -eo mmpbsa_energies.csv \
      > mmpbsa.log 2>&1
  )

  echo "Completed $system"
done

echo ""
echo "All MM-GBSA calculations finished."
