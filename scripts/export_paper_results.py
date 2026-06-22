#!/usr/bin/env python3
"""Export lightweight numerical source data for the CKP-25 manuscript.

The exporter reads local, ignored GROMACS/AmberTools output folders and writes
small CSV files suitable for GitHub. It deliberately does not copy raw docking
logs, trajectories, checkpoint files, binary run files, or full simulation
outputs.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import mdtraj as md
import numpy as np


SYSTEMS = [
    {
        "key": "MERS",
        "folder": "MERS_B_Spike_md",
        "virus": "MERS-CoV",
        "pdb_id": "4KR0",
        "chain": "B",
        "protein": "Data/proteins/MERS_4KR0_B_spike.pdb",
        "pose": "Data/docked_poses/MERS_CKP25_top_pose.sdf",
        "affinity": -6.61,
    },
    {
        "key": "SARS1",
        "folder": "SARS-CoV-1_E_spike_md",
        "virus": "SARS-CoV-1",
        "pdb_id": "2AJF",
        "chain": "E",
        "protein": "Data/proteins/SARS1_2AJF_E_spike.pdb",
        "pose": "Data/docked_poses/SARS1_CKP25_top_pose.sdf",
        "affinity": -6.18,
    },
    {
        "key": "NL63",
        "folder": "hCoV_NL63_E_Spike_md",
        "virus": "HCoV-NL63",
        "pdb_id": "3KBH",
        "chain": "E",
        "protein": "Data/proteins/NL63_3KBH_E_spike.pdb",
        "pose": "Data/docked_poses/NL63_CKP25_top_pose.sdf",
        "affinity": -6.32,
    },
    {
        "key": "229E",
        "folder": "hCoV_229E_E_spike_md",
        "virus": "HCoV-229E",
        "pdb_id": "6ATK",
        "chain": "E",
        "protein": "Data/proteins/229E_6ATK_E_spike.pdb",
        "pose": "Data/docked_poses/229E_CKP25_top_pose.sdf",
        "affinity": -5.93,
    },
    {
        "key": "HKU1",
        "folder": "hCoV_HKU1_A_spike_md",
        "virus": "HCoV-HKU1",
        "pdb_id": "8Y7Y",
        "chain": "A",
        "protein": "Data/proteins/HKU1_8Y7Y_A_spike.pdb",
        "pose": "Data/docked_poses/HKU1_CKP25_top_pose.sdf",
        "affinity": -6.76,
    },
]

DOCKING_BOXES = {
    "MERS": {"center": (-43.0, 25.0, 29.0), "size": (59.0, 35.0, 23.0)},
    "SARS1": {"center": (9.0, -17.0, 70.0), "size": (17.0, 39.0, 29.0)},
    "NL63": {"center": (7.47, 0.6, 43.0), "size": (17.0, 36.0, 18.0)},
    "229E": {"center": (115.8, 99.0, 55.0), "size": (23.0, 27.0, 34.0)},
    "HKU1": {"center": (140.0, 131.0, 148.0), "size": (20.0, 25.0, 81.0)},
}

STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "HID", "HIE", "HIP", "CYX",
}
SOLVENT_AND_IONS = {"SOL", "WAT", "HOH", "NA", "CL", "K", "Na", "Cl", "ION"}
POLAR_ELEMENTS = {"N", "O", "S"}


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def detect_lig_resname(topology: md.Topology) -> str:
    names = set()
    for residue in topology.residues:
        if residue.name not in STANDARD_AMINO_ACIDS and residue.name not in SOLVENT_AND_IONS:
            names.add(residue.name)
    return sorted(names)[0] if names else "MOL"


def residue_name_matches(left: str, right: str) -> bool:
    if left == right:
        return True
    if left == "HIS" and right in {"HID", "HIE", "HIP", "HIS"}:
        return True
    if left == "CYS" and right in {"CYX", "CYS"}:
        return True
    return False


def original_residue_map(sim_residues, original_pdb: Path) -> dict[int, tuple[str, int]]:
    original = md.load(str(original_pdb)).topology
    original_residues = [res for res in original.residues if res.is_protein]
    sim_seq = [res.name for res in sim_residues]
    orig_seq = [res.name for res in original_residues]

    n, m = len(sim_seq), len(orig_seq)
    dp = np.zeros((n + 1, m + 1), dtype=int)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if residue_name_matches(sim_seq[i - 1], orig_seq[j - 1]):
                dp[i, j] = dp[i - 1, j - 1] + 1
            else:
                dp[i, j] = max(dp[i - 1, j], dp[i, j - 1])

    mapping: dict[int, tuple[str, int]] = {}
    i, j = n, m
    while i > 0 and j > 0:
        match = residue_name_matches(sim_seq[i - 1], orig_seq[j - 1])
        if match and dp[i, j] == dp[i - 1, j - 1] + 1:
            orig = original_residues[j - 1]
            mapping[i - 1] = (orig.name, int(orig.resSeq))
            i -= 1
            j -= 1
        elif dp[i, j] == dp[i - 1, j]:
            i -= 1
        else:
            j -= 1
    return mapping


def load_solute_trajectory(system: dict, results_dir: Path, stride: int):
    run_dir = results_dir / system["folder"]
    top_path = run_dir / "npt.gro"
    traj_path = run_dir / "md.xtc"
    if not top_path.exists() or not traj_path.exists():
        raise FileNotFoundError(f"Missing trajectory/topology for {system['folder']}")

    full_top = md.load_topology(str(top_path))
    lig_resname = detect_lig_resname(full_top)
    solute_atoms = full_top.select(f"protein or resname {lig_resname}")
    traj = md.load(str(traj_path), top=full_top, stride=stride, atom_indices=solute_atoms)
    traj.image_molecules(inplace=True)
    traj.center_coordinates()
    return traj, lig_resname


def export_docking(out_dir: Path) -> None:
    score_rows = []
    box_rows = []
    for system in SYSTEMS:
        score_rows.append({
            "system": system["key"],
            "virus": system["virus"],
            "pdb_id": system["pdb_id"],
            "chain": system["chain"],
            "selected_pose_file": system["pose"],
            "reported_docking_affinity_kcal_mol": f"{system['affinity']:.2f}",
        })
        box = DOCKING_BOXES[system["key"]]
        center_x, center_y, center_z = box["center"]
        size_x, size_y, size_z = box["size"]
        box_rows.append({
            "system": system["key"],
            "virus": system["virus"],
            "pdb_id": system["pdb_id"],
            "chain": system["chain"],
            "center_x_angstrom": f"{center_x:.2f}",
            "center_y_angstrom": f"{center_y:.2f}",
            "center_z_angstrom": f"{center_z:.2f}",
            "size_x_angstrom": f"{size_x:.2f}",
            "size_y_angstrom": f"{size_y:.2f}",
            "size_z_angstrom": f"{size_z:.2f}",
            "min_x_angstrom": f"{center_x - size_x / 2:.2f}",
            "max_x_angstrom": f"{center_x + size_x / 2:.2f}",
            "min_y_angstrom": f"{center_y - size_y / 2:.2f}",
            "max_y_angstrom": f"{center_y + size_y / 2:.2f}",
            "min_z_angstrom": f"{center_z - size_z / 2:.2f}",
            "max_z_angstrom": f"{center_z + size_z / 2:.2f}",
        })

    write_csv(
        out_dir / "docking" / "docking_scores.csv",
        ["system", "virus", "pdb_id", "chain", "selected_pose_file", "reported_docking_affinity_kcal_mol"],
        score_rows,
    )
    write_csv(
        out_dir / "docking" / "docking_boxes.csv",
        [
            "system",
            "virus",
            "pdb_id",
            "chain",
            "center_x_angstrom",
            "center_y_angstrom",
            "center_z_angstrom",
            "size_x_angstrom",
            "size_y_angstrom",
            "size_z_angstrom",
            "min_x_angstrom",
            "max_x_angstrom",
            "min_y_angstrom",
            "max_y_angstrom",
            "min_z_angstrom",
            "max_z_angstrom",
        ],
        box_rows,
    )


def export_md(results_dir: Path, data_dir: Path, out_dir: Path, stride: int, bins: int) -> None:
    protein_rmsd_rows = []
    ligand_rmsd_rows = []
    rmsf_rows = []
    rg_rows = []
    close_contact_rows = []
    polar_contact_rows = []
    residue_freq_rows = []
    residue_heatmap_rows = []
    summary_rows = []

    for system in SYSTEMS:
        print(f"Exporting MD processed data for {system['key']}...")
        traj, lig_resname = load_solute_trajectory(system, results_dir, stride)
        time_ns = traj.time / 1000.0
        top = traj.topology

        protein_ca = top.select("protein and name CA")
        ligand_heavy = top.select(f"resname {lig_resname} and not element H")
        protein_heavy = top.select("protein and not element H")

        ref = traj[0]
        traj.superpose(ref, atom_indices=protein_ca)
        protein_rmsd = md.rmsd(traj, ref, 0, atom_indices=protein_ca, precentered=True) * 10
        ligand_rmsd = (
            md.rmsd(traj, ref, 0, atom_indices=ligand_heavy, precentered=True) * 10
            if len(ligand_heavy)
            else np.zeros(len(traj))
        )
        rg = md.compute_rg(traj.atom_slice(top.select("protein"))) * 10
        rmsf = md.rmsf(traj, traj, 0, atom_indices=protein_ca) * 10

        for frame_idx, t_ns in enumerate(time_ns):
            common = {
                "system": system["key"],
                "frame_index": frame_idx,
                "time_ns": f"{t_ns:.4f}",
            }
            protein_rmsd_rows.append({**common, "protein_ca_rmsd_angstrom": f"{protein_rmsd[frame_idx]:.6f}"})
            ligand_rmsd_rows.append({**common, "ligand_heavy_atom_rmsd_angstrom": f"{ligand_rmsd[frame_idx]:.6f}"})
            rg_rows.append({**common, "radius_of_gyration_angstrom": f"{rg[frame_idx]:.6f}"})

        protein_residues = [res for res in top.residues if res.is_protein]
        residue_map = original_residue_map(protein_residues, data_dir / system["protein"].replace("Data/", ""))

        for rel_idx, atom_idx in enumerate(protein_ca, start=1):
            residue = top.atom(atom_idx).residue
            sim_res_idx = protein_residues.index(residue)
            mapped_name, mapped_number = residue_map.get(
                sim_res_idx, (residue.name, int(residue.resSeq))
            )
            rmsf_rows.append({
                "system": system["key"],
                "relative_residue_index": rel_idx,
                "residue_name": mapped_name,
                "residue_number": mapped_number,
                "rmsf_ca_angstrom": f"{rmsf[rel_idx - 1]:.6f}",
            })

        if len(ligand_heavy) and len(protein_heavy):
            pairs = np.array([[lig, prot] for lig in ligand_heavy for prot in protein_heavy])
            distances = md.compute_distances(traj, pairs)
            close_contacts = np.sum(distances < 0.35, axis=1)

            atom_to_residue_idx = {}
            for res_idx, residue in enumerate(protein_residues):
                for atom in residue.atoms:
                    if atom.element.symbol != "H":
                        atom_to_residue_idx[atom.index] = res_idx

            residue_contacts = np.zeros((len(traj), len(protein_residues)), dtype=bool)
            for pair_idx, (_, protein_atom) in enumerate(pairs):
                res_idx = atom_to_residue_idx.get(int(protein_atom))
                if res_idx is not None:
                    residue_contacts[:, res_idx] |= distances[:, pair_idx] <= 0.40
        else:
            close_contacts = np.zeros(len(traj))
            residue_contacts = np.zeros((len(traj), len(protein_residues)), dtype=bool)

        ligand_polar = np.array([
            atom.index for atom in top.atoms
            if atom.residue.name == lig_resname and atom.element.symbol in POLAR_ELEMENTS
        ])
        protein_polar = np.array([
            atom.index for atom in top.atoms
            if atom.residue.is_protein and atom.element.symbol in POLAR_ELEMENTS
        ])
        if len(ligand_polar) and len(protein_polar):
            polar_pairs = np.array([[lig, prot] for lig in ligand_polar for prot in protein_polar])
            polar_distances = md.compute_distances(traj, polar_pairs)
            polar_contacts = np.sum(polar_distances < 0.35, axis=1)
        else:
            polar_contacts = np.zeros(len(traj))

        for frame_idx, t_ns in enumerate(time_ns):
            common = {
                "system": system["key"],
                "frame_index": frame_idx,
                "time_ns": f"{t_ns:.4f}",
            }
            close_contact_rows.append({
                **common,
                "close_heavy_atom_contacts_3p5A": int(close_contacts[frame_idx]),
            })
            polar_contact_rows.append({
                **common,
                "polar_contacts_3p5A": int(polar_contacts[frame_idx]),
            })

        contact_freq = residue_contacts.mean(axis=0) * 100.0
        for sim_res_idx, freq in enumerate(contact_freq):
            residue = protein_residues[sim_res_idx]
            mapped_name, mapped_number = residue_map.get(
                sim_res_idx, (residue.name, int(residue.resSeq))
            )
            residue_freq_rows.append({
                "system": system["key"],
                "residue_name": mapped_name,
                "residue_number": mapped_number,
                "relative_residue_index": sim_res_idx + 1,
                "contact_frequency_percent_4A": f"{freq:.6f}",
            })

        top_residue_indices = np.argsort(contact_freq)[::-1][:10]
        frame_bins = np.array_split(np.arange(len(traj)), bins)
        for bin_idx, frame_indices in enumerate(frame_bins):
            if len(frame_indices) == 0:
                continue
            bin_start = float(time_ns[frame_indices[0]])
            bin_end = float(time_ns[frame_indices[-1]])
            for sim_res_idx in top_residue_indices:
                residue = protein_residues[int(sim_res_idx)]
                mapped_name, mapped_number = residue_map.get(
                    int(sim_res_idx), (residue.name, int(residue.resSeq))
                )
                label = f"{mapped_name}{mapped_number}"
                residue_heatmap_rows.append({
                    "system": system["key"],
                    "bin_index": bin_idx,
                    "bin_start_ns": f"{bin_start:.4f}",
                    "bin_end_ns": f"{bin_end:.4f}",
                    "residue_label": label,
                    "contact_frequency_percent_4A": f"{residue_contacts[frame_indices, sim_res_idx].mean() * 100.0:.6f}",
                })

        half = len(traj) // 2
        summary_rows.append({
            "system": system["key"],
            "ligand_resname": lig_resname,
            "protein_ca_rmsd_mean_angstrom": f"{protein_rmsd.mean():.6f}",
            "protein_ca_rmsd_plateau_second_half_angstrom": f"{protein_rmsd[half:].mean():.6f}",
            "ligand_rmsd_mean_angstrom": f"{ligand_rmsd.mean():.6f}",
            "ligand_rmsd_plateau_second_half_angstrom": f"{ligand_rmsd[half:].mean():.6f}",
            "radius_of_gyration_mean_angstrom": f"{rg.mean():.6f}",
            "close_contacts_mean_3p5A": f"{close_contacts.mean():.6f}",
            "polar_contacts_mean_3p5A": f"{polar_contacts.mean():.6f}",
        })

    md_out = out_dir / "md_processed"
    write_csv(md_out / "protein_rmsd.csv", ["system", "frame_index", "time_ns", "protein_ca_rmsd_angstrom"], protein_rmsd_rows)
    write_csv(md_out / "ligand_rmsd.csv", ["system", "frame_index", "time_ns", "ligand_heavy_atom_rmsd_angstrom"], ligand_rmsd_rows)
    write_csv(md_out / "rmsf.csv", ["system", "relative_residue_index", "residue_name", "residue_number", "rmsf_ca_angstrom"], rmsf_rows)
    write_csv(md_out / "radius_of_gyration.csv", ["system", "frame_index", "time_ns", "radius_of_gyration_angstrom"], rg_rows)
    write_csv(md_out / "close_contacts.csv", ["system", "frame_index", "time_ns", "close_heavy_atom_contacts_3p5A"], close_contact_rows)
    write_csv(md_out / "polar_contacts.csv", ["system", "frame_index", "time_ns", "polar_contacts_3p5A"], polar_contact_rows)
    write_csv(md_out / "residue_contact_frequencies.csv", ["system", "residue_name", "residue_number", "relative_residue_index", "contact_frequency_percent_4A"], residue_freq_rows)
    write_csv(md_out / "residue_contact_heatmap_top10.csv", ["system", "bin_index", "bin_start_ns", "bin_end_ns", "residue_label", "contact_frequency_percent_4A"], residue_heatmap_rows)
    write_csv(md_out / "summary.csv", [
        "system", "ligand_resname", "protein_ca_rmsd_mean_angstrom",
        "protein_ca_rmsd_plateau_second_half_angstrom", "ligand_rmsd_mean_angstrom",
        "ligand_rmsd_plateau_second_half_angstrom", "radius_of_gyration_mean_angstrom",
        "close_contacts_mean_3p5A", "polar_contacts_mean_3p5A",
    ], summary_rows)


def parse_mmgbsa_delta_frames(path: Path) -> list[dict]:
    lines = path.read_text().splitlines()
    start = lines.index("DELTA Energy Terms")
    header = next(csv.reader([lines[start + 1]]))
    rows = []
    for line in lines[start + 2:]:
        if not line.strip():
            break
        values = next(csv.reader([line]))
        rows.append(dict(zip(header, values)))
    return rows


def parse_mmgbsa_summary(path: Path) -> dict[str, tuple[str, str, str]]:
    components = [
        "VDWAALS", "EEL", "EGB", "ESURF",
        "DELTA G gas", "DELTA G solv", "DELTA TOTAL",
    ]
    parsed = {}
    in_delta = False
    for line in path.read_text().splitlines():
        if "Differences (Complex - Receptor - Ligand)" in line:
            in_delta = True
            continue
        if not in_delta:
            continue
        for component in components:
            if line.startswith(component):
                values = line[len(component):].split()
                if len(values) >= 3:
                    parsed[component] = (values[0], values[1], values[2])
    return parsed


def export_mmgbsa(results_dir: Path, out_dir: Path) -> None:
    per_frame_rows = []
    summary_rows = []
    for system in SYSTEMS:
        run_dir = results_dir / system["folder"]
        frame_path = run_dir / "mmpbsa_energies.csv"
        summary_path = run_dir / "FINAL_RESULTS_MMPBSA.dat"
        if frame_path.exists():
            for row in parse_mmgbsa_delta_frames(frame_path):
                per_frame_rows.append({
                    "system": system["key"],
                    "frame_index": row["Frame #"],
                    "vdwaals_kcal_mol": row["VDWAALS"],
                    "eel_kcal_mol": row["EEL"],
                    "egb_kcal_mol": row["EGB"],
                    "esurf_kcal_mol": row["ESURF"],
                    "delta_g_gas_kcal_mol": row["DELTA G gas"],
                    "delta_g_solv_kcal_mol": row["DELTA G solv"],
                    "delta_total_kcal_mol": row["DELTA TOTAL"],
                })
        if summary_path.exists():
            parsed = parse_mmgbsa_summary(summary_path)
            for component, values in parsed.items():
                avg, std, sem = values
                summary_rows.append({
                    "system": system["key"],
                    "energy_component": component,
                    "average_kcal_mol": avg,
                    "std_dev_kcal_mol": std,
                    "std_err_mean_kcal_mol": sem,
                })

    mmgbsa_out = out_dir / "mmgbsa"
    write_csv(mmgbsa_out / "per_frame_energies.csv", [
        "system", "frame_index", "vdwaals_kcal_mol", "eel_kcal_mol",
        "egb_kcal_mol", "esurf_kcal_mol", "delta_g_gas_kcal_mol",
        "delta_g_solv_kcal_mol", "delta_total_kcal_mol",
    ], per_frame_rows)
    write_csv(mmgbsa_out / "summary.csv", [
        "system", "energy_component", "average_kcal_mol",
        "std_dev_kcal_mol", "std_err_mean_kcal_mol",
    ], summary_rows)


def write_readmes(out_dir: Path, stride: int, bins: int) -> None:
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)
    (out_dir / "README.md").write_text(
        "# Processed Paper Results\n\n"
        "This directory contains lightweight numerical source data exported from "
        "the local MD/MM-GBSA analyses. Raw trajectories, binary GROMACS files, "
        "checkpoints, raw docking logs, and full intermediate outputs are not "
        "included.\n\n"
        f"MD trajectory-derived tables were exported with stride={stride}. "
        f"Residue contact heatmap tables use {bins} time bins and a 4 angstrom "
        "protein-ligand heavy-atom cutoff. Close-contact and polar-contact time "
        "series use a 3.5 angstrom cutoff.\n\n"
        "Subdirectories:\n\n"
        "- `docking/`: processed selected-pose docking scores and docking "
        "search-box coordinates; raw docking logs are not included.\n"
        "- `md_processed/`: numerical source data for RMSD, RMSF, radius of "
        "gyration, close contacts, polar contacts, and residue contact maps.\n"
        "- `mmgbsa/`: final MM-GBSA summary terms and per-frame delta energies.\n"
        "- `figures/`: placeholder for regenerated figures; publication figures "
        "can be recreated from the workflow scripts.\n"
    )
    (out_dir / "figures" / "README.md").write_text(
        "# Figures\n\n"
        "Final figures are generated by the analysis scripts from the processed "
        "tables and/or local trajectories. They are not duplicated here to keep "
        "the repository compact.\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--data-dir", default="Data")
    parser.add_argument("--out-dir", default="paper_results")
    parser.add_argument("--stride", type=int, default=10)
    parser.add_argument("--contact-bins", type=int, default=50)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    export_docking(out_dir)
    export_md(results_dir, data_dir, out_dir, args.stride, args.contact_bins)
    export_mmgbsa(results_dir, out_dir)
    write_readmes(out_dir, args.stride, args.contact_bins)
    print(f"Processed paper results written to {out_dir}")


if __name__ == "__main__":
    main()
