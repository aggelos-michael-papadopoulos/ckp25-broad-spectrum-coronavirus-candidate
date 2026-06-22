#!/usr/bin/env python3
import argparse
from pathlib import Path

import mdtraj as md
import numpy as np

DEFAULT_SYSTEMS = [
    "MERS_B_Spike_md",
    "SARS-CoV-1_E_spike_md",
    "hCoV_229E_E_spike_md",
    "hCoV_HKU1_A_spike_md",
    "hCoV_NL63_E_Spike_md",
]


def detect_lig_resname(top) -> str:
    non_std = set()
    for res in top.residues:
        if not res.is_protein and res.name not in ("HOH", "SOL", "WAT", "NA", "CL", "K"):
            non_std.add(res.name)
    if len(non_std) == 1:
        return non_std.pop()
    if non_std:
        return sorted(non_std)[0]
    return "LIG"


def summarize_system(sys_path: Path, stride: int) -> str:
    traj_path = sys_path / "md.xtc"
    top_path = sys_path / "npt.gro"

    if not traj_path.exists() or not top_path.exists():
        return f"| {sys_path.name} | Files not found | - | - | - | - | - | - | - |"

    try:
        traj = md.load(str(traj_path), top=str(top_path), stride=stride)
        traj.image_molecules(inplace=True)
        traj.center_coordinates()

        top = traj.topology
        lig_resname = detect_lig_resname(top)

        ca_idx = top.select("name CA")
        ref = traj[0]
        traj.superpose(ref, atom_indices=ca_idx)
        rmsd_prot = md.rmsd(traj, ref, 0, atom_indices=ca_idx, precentered=True) * 10

        lig_idx = top.select(f"resname {lig_resname} and not element H")
        if len(lig_idx) > 0:
            rmsd_lig = md.rmsd(traj, ref, 0, atom_indices=lig_idx, precentered=True) * 10
        else:
            rmsd_lig = np.zeros(len(traj))

        half_idx = len(traj) // 2
        prot_plateau = np.mean(rmsd_prot[half_idx:])
        prot_mean = np.mean(rmsd_prot)
        lig_plateau = np.mean(rmsd_lig[half_idx:]) if len(lig_idx) > 0 else 0.0
        lig_mean = np.mean(rmsd_lig) if len(lig_idx) > 0 else 0.0

        rg = md.compute_rg(traj.atom_slice(top.select("protein"))) * 10
        rg_mean = np.mean(rg)

        lig_atoms = [a.index for a in top.atoms if a.residue.name == lig_resname]
        lig_set = set(lig_atoms)
        prot_atoms_set = set(a.index for a in top.atoms if a.residue.is_protein)

        polar_elements = {"N", "O", "S"}
        lig_polar = np.array([
            a.index for a in top.atoms
            if a.residue.name == lig_resname and a.element.symbol in polar_elements
        ])
        prot_polar = np.array([
            a.index for a in top.atoms
            if a.residue.is_protein and a.element.symbol in polar_elements
        ])

        hbonds = []
        contacts = []
        polar_pairs = None
        if len(lig_polar) > 0 and len(prot_polar) > 0:
            polar_pairs = np.array([[lp, pp] for lp in lig_polar for pp in prot_polar])

        for frame in traj:
            hb = md.baker_hubbard(frame, periodic=False)
            count = sum(
                1
                for d, h, a in hb
                if (d in prot_atoms_set and a in lig_set)
                or (d in lig_set and a in prot_atoms_set)
            )
            hbonds.append(count)

            if polar_pairs is not None:
                dists = md.compute_distances(frame, polar_pairs)[0]
                contacts.append(int((dists < 0.35).sum()))
            else:
                contacts.append(0)

        hb_mean = np.mean(hbonds)
        ct_mean = np.mean(contacts)

        return (
            f"| {sys_path.name} | {lig_resname} | {prot_plateau:.2f} | "
            f"{prot_mean:.2f} | {lig_plateau:.2f} | {lig_mean:.2f} | "
            f"{rg_mean:.2f} | {hb_mean:.2f} | {ct_mean:.2f} |"
        )

    except Exception as exc:
        return f"| {sys_path.name} | Error: {exc} | - | - | - | - | - | - | - |"


def main():
    parser = argparse.ArgumentParser(
        description="Print a Markdown summary table of MD stability/contact metrics."
    )
    parser.add_argument("--results-dir", default="runs")
    parser.add_argument("--systems", nargs="*", default=DEFAULT_SYSTEMS)
    parser.add_argument("--stride", type=int, default=50)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    print("| System | Ligand Resname | Prot Cα RMSD Plateau (Å) | Prot Cα RMSD Mean (Å) | Lig RMSD Plateau (Å) | Lig RMSD Mean (Å) | Rg Mean (Å) | H-bonds Mean | Polar Contacts Mean |")
    print("|---|---|---|---|---|---|---|---|---|")
    for sys_name in args.systems:
        print(summarize_system(results_dir / sys_name, args.stride))


if __name__ == "__main__":
    main()
