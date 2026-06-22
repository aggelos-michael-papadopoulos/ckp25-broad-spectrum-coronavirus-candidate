#!/usr/bin/env bash
# run_MD_gromacs.sh — One-liner GROMACS protein-ligand MD pipeline
#
# Usage:
#   bash run_MD_gromacs.sh <target.pdb> <compound.sdf> <ns> [outdir]
#
# Example:
#   bash run_MD_gromacs.sh 4aqp_receptor.pdb top_pose_1.sdf 100
#   bash run_MD_gromacs.sh 4aqp_receptor.pdb top_pose_1.sdf 10 results/my_run

set -euo pipefail

T_START=$(date +%s)

TARGET="${1:?Usage: $0 <target.pdb> <compound.sdf> <ns> [outdir]}"
COMPOUND="${2:?Usage: $0 <target.pdb> <compound.sdf> <ns> [outdir]}"
NS="${3:?Usage: $0 <target.pdb> <compound.sdf> <ns> [outdir]}"
OUTDIR="${4:-results/$(basename "${TARGET%.pdb}")_$(basename "${COMPOUND%.sdf}")_GMX}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "=============================================="
echo " GROMACS MD Pipeline"
echo " Target   : $TARGET"
echo " Compound : $COMPOUND"
echo " Sim time : ${NS} ns"
echo " Output   : $OUTDIR"
echo "=============================================="

# ── Step 1: Prepare complex ─────────────────────────────────────────────────
echo ""
echo "[1/4] Preparing complex..."
python "$SCRIPT_DIR/prepare_complex_gmx.py" \
    --target   "$TARGET"   \
    --compound "$COMPOUND" \
    --outdir   "$OUTDIR"

# ── Step 2: Run simulation ──────────────────────────────────────────────────
echo ""
echo "[2/4] Running simulation (${NS} ns)..."
bash "$SCRIPT_DIR/run_simulation_gmx.sh" "$OUTDIR" "$NS"

# ── Step 3: Analyse ─────────────────────────────────────────────────────────
echo ""
echo "[3/4] Analysing trajectory..."
python "$SCRIPT_DIR/analyze_gmx.py" --outdir "$OUTDIR"

# ── Step 4: HTML viewer ──────────────────────────────────────────────────────
echo ""
echo "[4/4] Generating HTML viewer..."
python "$SCRIPT_DIR/make_viewer_gmx.py" --outdir "$OUTDIR"

T_END=$(date +%s)
T_ELAPSED=$(( T_END - T_START ))
T_H=$(( T_ELAPSED / 3600 ))
T_M=$(( (T_ELAPSED % 3600) / 60 ))
T_S=$(( T_ELAPSED % 60 ))

echo ""
echo "=============================================="
echo " Done! Results in: $OUTDIR"
echo "   Plots  : plot_A_hbonds.png  plot_B_rmsd.png"
echo "            plot_C_rmsf.png    plot_D_rg.png"
echo "   Viewer : trajectory_viewer_gmx.html"
echo "   Traj   : md.xtc"
echo "----------------------------------------------"
printf "   Total time : %02dh %02dm %02ds\n" $T_H $T_M $T_S
echo "=============================================="
