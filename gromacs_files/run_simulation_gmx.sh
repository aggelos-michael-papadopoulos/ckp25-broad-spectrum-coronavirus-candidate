#!/usr/bin/env bash
# run_simulation_gmx.sh — EM → NVT → NPT → Production MD
#
# Usage:
#   bash run_simulation_gmx.sh <outdir> <ns>
#
# Expects outdir to already contain:
#   complex_ions.gro  topol.top  posre.itp
#
# GPU: single GPU (id 0), 1 MPI rank, 4 OpenMP threads.
# Adjust -ntomp to match your CPU core count if needed.

set -euo pipefail

OUTDIR="${1:?Usage: $0 <outdir> <ns>}"
NS="${2:?Usage: $0 <outdir> <ns>}"
NSTEPS=$(python3 -c "print(int(${NS} * 1e6 / 2))")   # 2 fs timestep

GMX="gmx"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="$(cd "$OUTDIR" && pwd)"
MDP_DIR="${OUTDIR}/mdp"
mkdir -p "$MDP_DIR"

# ── Write MDP files ────────────────────────────────────────────────────────

cat > "$MDP_DIR/em.mdp" << 'EOF'
; Energy minimisation
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000

nstlist         = 1
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
EOF

cat > "$MDP_DIR/nvt.mdp" << 'EOF'
; NVT equilibration — 100 ps
integrator              = md
dt                      = 0.002
nsteps                  = 50000
nstxout-compressed      = 500
nstlog                  = 500
nstenergy               = 500

cutoff-scheme           = Verlet
nstlist                 = 10
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres

tcoupl                  = V-rescale
tc-grps                 = Protein_LIG Water_and_ions
tau_t                   = 0.1  0.1
ref_t                   = 300  300

pcoupl                  = no

gen-vel                 = yes
gen-temp                = 300
gen-seed                = -1

constraints             = h-bonds
constraint-algorithm    = LINCS
lincs-iter              = 1
lincs-order             = 4

; Position restraints on protein heavy atoms
define                  = -DPOSRES
pbc                     = xyz
EOF

cat > "$MDP_DIR/npt.mdp" << 'EOF'
; NPT equilibration — 100 ps
integrator              = md
dt                      = 0.002
nsteps                  = 50000
nstxout-compressed      = 500
nstlog                  = 500
nstenergy               = 500

cutoff-scheme           = Verlet
nstlist                 = 10
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres

tcoupl                  = V-rescale
tc-grps                 = Protein_LIG Water_and_ions
tau_t                   = 0.1  0.1
ref_t                   = 300  300

pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5

gen-vel                 = no

constraints             = h-bonds
constraint-algorithm    = LINCS

define                  = -DPOSRES
pbc                     = xyz
EOF

# Production MDP — nsteps computed from NS argument
cat > "$MDP_DIR/md.mdp" << EOF
; Production MD
integrator              = md
dt                      = 0.002
nsteps                  = ${NSTEPS}
nstxout-compressed      = 5000
nstlog                  = 5000
nstenergy               = 5000

cutoff-scheme           = Verlet
nstlist                 = 10
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres

tcoupl                  = V-rescale
tc-grps                 = Protein_LIG Water_and_ions
tau_t                   = 0.1  0.1
ref_t                   = 300  300

pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5

gen-vel                 = no

constraints             = h-bonds
constraint-algorithm    = LINCS

pbc                     = xyz
EOF

cd "$OUTDIR"

# ── Make a Protein_LIG group for thermostats ──────────────────────────────
if ! grep -q '\[ Protein_LIG \]' index.ndx 2>/dev/null; then
    # Create default index if it doesn't exist yet
    if [ ! -f index.ndx ]; then
        echo "q" | $GMX make_ndx -f complex_ions.gro -o index.ndx 2>/dev/null
    fi
    # Combine Protein (1) and MOL/ligand (13) then rename
    echo -e "1 | 13\nname 21 Protein_LIG\nq" | \
        $GMX make_ndx -f complex_ions.gro -o index.ndx -n index.ndx 2>/dev/null
fi

PROGRESS="python3 $SCRIPT_DIR/_mdrun_progress.py"

# ── Energy minimisation  [1/4] ────────────────────────────────────────────
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  [1/4] Energy Minimisation  (steepest descent)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
$GMX grompp -f "$MDP_DIR/em.mdp" -c complex_ions.gro -p topol.top -o em.tpr -maxwarn 20 2>&1 | grep -E "NOTE|WARNING|Error|Fatal" || true
rm -f em.log
$GMX mdrun -deffnm em -ntmpi 1 -ntomp 4 -nb gpu -gpu_id 0 > /dev/null 2>&1 &
$PROGRESS em.log 50000 "  EM       " 0.01
wait

# ── NVT equilibration  [2/4] ──────────────────────────────────────────────
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  [2/4] NVT Equilibration  (100 ps, T=300 K)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
$GMX grompp -f "$MDP_DIR/nvt.mdp" -c em.gro -r em.gro -p topol.top \
            -n index.ndx -o nvt.tpr -maxwarn 20 2>&1 | grep -E "NOTE|WARNING|Error|Fatal" || true
rm -f nvt.log
$GMX mdrun -deffnm nvt -ntmpi 1 -ntomp 4 -nb gpu -pme gpu -gpu_id 0 > /dev/null 2>&1 &
$PROGRESS nvt.log 50000 "  NVT      " 0.002
wait

# ── NPT equilibration  [3/4] ──────────────────────────────────────────────
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  [3/4] NPT Equilibration  (100 ps, P=1 bar)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
$GMX grompp -f "$MDP_DIR/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt \
            -p topol.top -n index.ndx -o npt.tpr -maxwarn 20 2>&1 | grep -E "NOTE|WARNING|Error|Fatal" || true
rm -f npt.log
$GMX mdrun -deffnm npt -ntmpi 1 -ntomp 4 -nb gpu -pme gpu -gpu_id 0 > /dev/null 2>&1 &
$PROGRESS npt.log 50000 "  NPT      " 0.002
wait

# ── Production MD  [4/4] ──────────────────────────────────────────────────
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  [4/4] Production MD  (${NS} ns)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
$GMX grompp -f "$MDP_DIR/md.mdp" -c npt.gro -t npt.cpt \
            -p topol.top -n index.ndx -o md.tpr -maxwarn 20 2>&1 | grep -E "NOTE|WARNING|Error|Fatal" || true
rm -f md.log
$GMX mdrun -deffnm md -ntmpi 1 -ntomp 4 -nb gpu -pme gpu -bonded gpu -gpu_id 0 > /dev/null 2>&1 &
$PROGRESS md.log "$NSTEPS" "  MD ${NS} ns" 0.002
wait

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Simulation complete: $(pwd)"
echo "  Trajectory : md.xtc    Energy : md.edr"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
