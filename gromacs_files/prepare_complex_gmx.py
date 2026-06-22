#!/usr/bin/env python3
"""
prepare_complex_gmx.py
Prepare a protein-ligand complex for GROMACS MD simulation.

Steps:
  1. Fix protein PDB with PDBFixer
  2. Run pdb2gmx (AMBER99SB-ILDN + TIP3P) on the protein
  3. Parametrize ligand with ACPYPE
  4. Merge topologies and coordinates into complex.gro / topol.top
  5. Define simulation box and solvate
  6. Add counterions

Usage:
  python prepare_complex_gmx.py --target protein.pdb --compound ligand.sdf [--outdir results/run]
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd, cwd=None, check=True):
    print(f"  $ {cmd}")
    result = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    if check and result.returncode != 0:
        sys.exit(f"Command failed (exit {result.returncode}): {cmd}")
    return result


def fix_protein(pdb_in: Path, pdb_out: Path):
    """Use PDBFixer to add missing atoms/residues and remove HETATM."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    print("[1] Fixing protein PDB with PDBFixer...")
    fixer = PDBFixer(filename=str(pdb_in))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    with open(pdb_out, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    print(f"   → {pdb_out}")


def get_ligand_charge(sdf_path: Path) -> int:
    """Return formal charge of the first molecule in the SDF."""
    try:
        from rdkit import Chem
        mol = Chem.SDMolSupplier(str(sdf_path), removeHs=False)[0]
        if mol is None:
            return 0
        return sum(a.GetFormalCharge() for a in mol.GetAtoms())
    except Exception:
        return 0


def split_acpype_itp(itp_path: Path, atomtypes_itp: Path, molecule_itp: Path):
    """
    Split ACPYPE _GMX.itp into:
      - atomtypes_itp : [ atomtypes ] section only
      - molecule_itp  : everything else (no atomtypes)
    This avoids 'duplicate atomtype' errors when including into topol.top.
    """
    with open(itp_path) as f:
        content = f.read()

    # Extract [ atomtypes ] block
    at_match = re.search(r'(\[ atomtypes \].*?)(?=\[|\Z)', content, re.DOTALL)
    atomtypes_block = at_match.group(1) if at_match else ""

    # Remove atomtypes block from main itp
    mol_content = re.sub(r'\[ atomtypes \].*?(?=\[|\Z)', '', content, flags=re.DOTALL).strip()

    with open(atomtypes_itp, "w") as f:
        f.write(atomtypes_block)
    with open(molecule_itp, "w") as f:
        f.write(mol_content + "\n")


def merge_topology(topol_path: Path, atomtypes_itp: Path, molecule_itp: Path,
                   lig_name: str, outdir: Path):
    """Insert ligand includes and molecule entry into topol.top."""
    with open(topol_path) as f:
        lines = f.readlines()

    new_lines = []
    inserted_at = inserted_mol = False

    for line in lines:
        # Insert atomtypes include right after forcefield include
        if not inserted_at and '#include' in line and 'forcefield.itp' in line:
            new_lines.append(line)
            new_lines.append(f'#include "{atomtypes_itp.name}"\n')
            new_lines.append(f'#include "{molecule_itp.name}"\n')
            inserted_at = True
            continue

        # Append ligand molecule count at end of [ molecules ]
        if not inserted_mol and line.strip().startswith('[ molecules ]'):
            new_lines.append(line)
            # collect rest until EOF
            idx = lines.index(line)
            remaining = lines[idx + 1:]
            new_lines.extend(remaining)
            new_lines.append(f'{lig_name}             1\n')
            inserted_mol = True
            break

        new_lines.append(line)

    with open(topol_path, "w") as f:
        f.writelines(new_lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Prepare protein-ligand complex for GROMACS")
    parser.add_argument("--target",   required=True, help="Input protein PDB")
    parser.add_argument("--compound", required=True, help="Input ligand SDF")
    parser.add_argument("--outdir",   default=None,  help="Output directory (default: auto)")
    args = parser.parse_args()

    target_pdb  = Path(args.target).resolve()
    compound_sdf = Path(args.compound).resolve()

    if args.outdir:
        outdir = Path(args.outdir)
    else:
        stem_t = target_pdb.stem
        stem_c = compound_sdf.stem
        outdir = Path("results") / f"{stem_t}_{stem_c}_GMX"

    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {outdir}\n")

    # ------------------------------------------------------------------ #
    # Step 1 — Fix protein
    # ------------------------------------------------------------------ #
    fixed_pdb = outdir / "protein_fixed.pdb"
    fix_protein(target_pdb, fixed_pdb)

    # ------------------------------------------------------------------ #
    # Step 2 — pdb2gmx (non-interactive)
    # ------------------------------------------------------------------ #
    print("\n[2] Running pdb2gmx...")
    run(
        f"gmx pdb2gmx "
        f"-f protein_fixed.pdb "
        f"-o protein.gro "
        f"-p topol.top "
        f"-i posre.itp "
        f"-ff amber99sb-ildn "
        f"-water tip3p "
        f"-ignh",
        cwd=outdir,
    )

    # ------------------------------------------------------------------ #
    # Step 3 — Parametrize ligand with ACPYPE
    # ------------------------------------------------------------------ #
    print("\n[3] Running ACPYPE on ligand...")
    charge = get_ligand_charge(compound_sdf)
    print(f"   Detected ligand formal charge: {charge}")

    lig_sdf_local = outdir / compound_sdf.name
    shutil.copy(compound_sdf, lig_sdf_local)

    # Skip if ACPYPE already ran (saves ~5 min sqm re-computation)
    acpype_dir = next(outdir.glob("*.acpype"), None)
    if acpype_dir and next(acpype_dir.glob("*_GMX.itp"), None):
        print(f"   ACPYPE output already exists ({acpype_dir.name}) — skipping.")
    else:
        run(
            f"acpype -i {compound_sdf.name} -c bcc -n {charge} -a gaff2",
            cwd=outdir,
        )
        acpype_dir = next(outdir.glob("*.acpype"), None)
        if acpype_dir is None:
            sys.exit("ACPYPE output directory not found.")

    lig_stem = acpype_dir.stem.replace(".acpype", "")
    raw_itp  = next(acpype_dir.glob("*_GMX.itp"), None)
    raw_gro  = next(acpype_dir.glob("*_GMX.gro"), None)

    if raw_itp is None or raw_gro is None:
        sys.exit("ACPYPE did not produce expected _GMX.itp / _GMX.gro files.")

    # Detect molecule name inside itp
    # Find the actual molecule name — skip comment lines starting with ';'
    lig_name = lig_stem
    in_moltype = False
    for _line in raw_itp.read_text().splitlines():
        if '[ moleculetype ]' in _line:
            in_moltype = True
            continue
        if in_moltype:
            _s = _line.strip()
            if _s and not _s.startswith(';'):
                lig_name = _s.split()[0]
                break

    # Copy ligand gro
    lig_gro = outdir / "ligand.gro"
    shutil.copy(raw_gro, lig_gro)

    # Split itp
    atomtypes_itp = outdir / "lig_atomtypes.itp"
    molecule_itp  = outdir / "lig_molecule.itp"
    split_acpype_itp(raw_itp, atomtypes_itp, molecule_itp)
    print(f"   → {atomtypes_itp.name}, {molecule_itp.name}")

    # ------------------------------------------------------------------ #
    # Step 4 — Merge protein + ligand coordinates
    # ------------------------------------------------------------------ #
    print("\n[4] Merging protein and ligand coordinates...")

    def gro_lines(path):
        with open(path) as f:
            return f.readlines()

    prot_lines = gro_lines(outdir / "protein.gro")
    lig_lines  = gro_lines(lig_gro)

    # gro format: line 0 = title, line 1 = natoms, lines 2..-1 = atoms, last = box
    prot_atoms = prot_lines[2:-1]
    lig_atoms  = lig_lines[2:-1]
    box_line   = prot_lines[-1]

    n_total = len(prot_atoms) + len(lig_atoms)

    complex_gro = outdir / "complex.gro"
    with open(complex_gro, "w") as f:
        f.write("Protein-Ligand complex\n")
        f.write(f"{n_total}\n")
        f.writelines(prot_atoms)
        f.writelines(lig_atoms)
        f.write(box_line)
    print(f"   → {complex_gro}")

    # ------------------------------------------------------------------ #
    # Step 5 — Patch topology
    # ------------------------------------------------------------------ #
    print("\n[5] Patching topol.top with ligand...")
    merge_topology(outdir / "topol.top", atomtypes_itp, molecule_itp, lig_name, outdir)

    # ------------------------------------------------------------------ #
    # Step 6 — Define box and solvate
    # ------------------------------------------------------------------ #
    print("\n[6] Defining simulation box (1.2 nm margin)...")
    run(
        "gmx editconf -f complex.gro -o complex_box.gro -c -d 1.2 -bt dodecahedron",
        cwd=outdir,
    )

    print("\n[7] Solvating...")
    run(
        "gmx solvate -cp complex_box.gro -cs spc216.gro -o complex_solv.gro -p topol.top",
        cwd=outdir,
    )

    # ------------------------------------------------------------------ #
    # Step 7 — Add ions
    # ------------------------------------------------------------------ #
    print("\n[8] Adding ions (150 mM NaCl, neutralising)...")
    ions_mdp = outdir / "ions.mdp"
    ions_mdp.write_text(
        "; Minimal mdp for grompp before adding ions\n"
        "integrator = steep\n"
        "nsteps     = 0\n"
        "emtol      = 1000\n"
        "emstep     = 0.01\n"
        "nstlist    = 1\n"
        "cutoff-scheme = Verlet\n"
        "ns_ns_type    = grid\n"
        "coulombtype   = cutoff\n"
        "rcoulomb      = 1.0\n"
        "rvdw          = 1.0\n"
        "pbc           = xyz\n"
    )

    run(
        "gmx grompp -f ions.mdp -c complex_solv.gro -p topol.top "
        "-o ions.tpr -maxwarn 20",
        cwd=outdir,
    )
    # echo SOL | gmx genion  (non-interactive via stdin)
    run(
        'echo "SOL" | gmx genion -s ions.tpr -o complex_ions.gro '
        '-p topol.top -pname NA -nname CL -neutral -conc 0.15',
        cwd=outdir,
    )

    print("\n[DONE] Complex ready for simulation.")
    print(f"  Coordinates : {outdir}/complex_ions.gro")
    print(f"  Topology    : {outdir}/topol.top")


if __name__ == "__main__":
    main()
