#!/usr/bin/env bash
# Build Amber-compatible MM-GBSA inputs from finished GROMACS simulations.
#
# Usage:
#   bash scripts/prepare_mmpbsa_inputs.sh [results_dir]

set -euo pipefail

RESULTS_DIR="${1:-runs}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

python "$ROOT/convert_topologies.py" --results-dir "$ROOT/$RESULTS_DIR"
python "$ROOT/extract_mmpbsa_traj.py" --results-dir "$ROOT/$RESULTS_DIR"

echo ""
echo "MM-GBSA input preparation finished."
echo "Expected files per system: complex.prmtop, receptor.prmtop, ligand.prmtop, complex.nc"
