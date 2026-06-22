#!/usr/bin/env python3
import argparse
from pathlib import Path

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np

DEFAULT_SYSTEMS = [
    ("MERS_B_Spike_md", "MERS-CoV B-spike (PDB: 4KR0)", "mers_metrics.png"),
    ("SARS-CoV-1_E_spike_md", "SARS-CoV-1 E-spike (PDB: 2AJF)", "sars_metrics.png"),
    ("hCoV_NL63_E_Spike_md", "HCoV-NL63 E-spike (PDB: 3KBH)", "nl63_metrics.png"),
    ("hCoV_229E_E_spike_md", "HCoV-229E E-spike (PDB: 6ATK)", "229e_metrics.png"),
    ("hCoV_HKU1_A_spike_md", "HCoV-HKU1 A-spike (PDB: 8Y7Y)", "hku1_metrics.png"),
]

STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "HID", "HIE", "HIP", "CYX",
}

SOLVENT_AND_IONS = {"SOL", "WAT", "HOH", "NA", "CL", "K", "Na", "Cl", "ION"}


plt.style.use("default")
plt.rcParams.update({
    "font.size": 14,
    "axes.linewidth": 1.5,
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "xtick.major.width": 1.5,
    "ytick.major.width": 1.5,
    "lines.linewidth": 1.5,
})


def detect_lig_resname(top) -> str:
    lig_resnames = set()
    for residue in top.residues:
        if residue.name not in STANDARD_AMINO_ACIDS and residue.name not in SOLVENT_AND_IONS:
            lig_resnames.add(residue.name)
    return sorted(lig_resnames)[0] if lig_resnames else "MOL"


def plot_system_metrics(sys_dir: Path, display_title: str, out_img_path: Path,
                        total_ns: float, stride: int) -> bool:
    traj_file = sys_dir / "md.xtc"
    top_file = sys_dir / "npt.gro"

    if not traj_file.exists() or not top_file.exists():
        print(f"Skipping {sys_dir.name}: trajectory or topology file not found.")
        return False

    print(f"Loading topology for {sys_dir.name}...")
    top = md.load_topology(str(top_file))
    lig_resname = detect_lig_resname(top)
    print(f"  Detected ligand name: {lig_resname}")

    print(f"Loading trajectory for {sys_dir.name}...")
    sel = top.select(f"protein or resname {lig_resname}")
    traj = md.load(str(traj_file), top=top, stride=stride, atom_indices=sel)
    time_ns = np.linspace(0, total_ns, len(traj))

    fig, axs = plt.subplots(2, 2, figsize=(14, 10))

    print("  Computing RMSD...")
    lig_sel = traj.topology.select(f"resname {lig_resname} and not element H")
    prot_sel = traj.topology.select("protein and name CA")
    traj.superpose(traj, 0, atom_indices=prot_sel)
    prot_rmsd = md.rmsd(traj, traj, 0, atom_indices=prot_sel) * 10
    lig_rmsd = md.rmsd(traj, traj, 0, atom_indices=lig_sel) * 10 if len(lig_sel) else np.zeros(len(traj))

    axs[0, 0].plot(time_ns, prot_rmsd, color="#1f77b4", linewidth=1.2, alpha=0.85, label="Protein Cα")
    axs[0, 0].plot(time_ns, lig_rmsd, color="#ff7f0e", linewidth=1.2, alpha=0.85, label="Ligand")
    axs[0, 0].set_title("RMSD from Starting Structure", fontweight="bold")
    axs[0, 0].set_ylabel("RMSD (Å)")
    axs[0, 0].set_xlabel("Time (ns)")
    y_limit = max(float(np.max(prot_rmsd)) * 1.1, 2.4)
    max_lig_rmsd = float(np.max(lig_rmsd)) if len(lig_rmsd) else 0.0
    axs[0, 0].set_ylim(0, y_limit)
    if max_lig_rmsd > y_limit:
        axs[0, 0].text(
            0.98, 0.95, f"max ligand spike: {max_lig_rmsd:.1f} Å (clipped)",
            transform=axs[0, 0].transAxes, ha="right", va="top",
            fontsize=10, color="gray", style="italic",
        )
    axs[0, 0].legend(loc="lower right", frameon=True, fontsize=12)

    print("  Computing RMSF...")
    rmsf = md.rmsf(traj, traj, 0, atom_indices=prot_sel) * 10
    res_indices = [traj.topology.atom(a).residue.resSeq for a in prot_sel]
    axs[0, 1].plot(res_indices, rmsf, color="#1f77b4", linewidth=1.5)
    axs[0, 1].set_title("Protein Cα RMSF", fontweight="bold")
    axs[0, 1].set_ylabel("RMSF (Å)")
    axs[0, 1].set_xlabel("Residue Index")
    axs[0, 1].set_ylim(0, max(float(np.max(rmsf)) * 1.1, 2.0))

    print("  Computing Radius of Gyration...")
    traj_prot = traj.atom_slice(traj.topology.select("protein"))
    rg = md.compute_rg(traj_prot) * 10
    axs[1, 0].plot(time_ns, rg, color="#2ca02c", linewidth=1.5)
    axs[1, 0].set_title("Protein Radius of Gyration", fontweight="bold")
    axs[1, 0].set_ylabel("Rg (Å)")
    axs[1, 0].set_xlabel("Time (ns)")
    axs[1, 0].set_ylim(float(np.min(rg)) - 0.5, float(np.max(rg)) + 0.5)

    print("  Computing Contact Frequencies...")
    lig_atoms = traj.topology.select(f"resname {lig_resname} and not element H")
    prot_atoms = traj.topology.select("protein and not element H")
    if len(lig_atoms) and len(prot_atoms):
        pairs = np.array([[la, pa] for la in lig_atoms for pa in prot_atoms])
        dists = md.compute_distances(traj, pairs)
        contacts = np.sum(dists < 0.35, axis=1)
    else:
        contacts = np.zeros(len(traj))

    axs[1, 1].plot(time_ns, contacts, color="#9467bd", linewidth=1.0, alpha=0.5, label="Raw")
    window = min(30, max(1, len(contacts)))
    rolling_contacts = np.convolve(contacts, np.ones(window) / window, mode="valid")
    axs[1, 1].plot(time_ns[window - 1:], rolling_contacts, color="#4b0082", linewidth=2.0, label="Moving Average")
    axs[1, 1].set_title("Protein-Ligand Contacts (<3.5 Å)", fontweight="bold")
    axs[1, 1].set_ylabel("Number of Contacts")
    axs[1, 1].set_xlabel("Time (ns)")
    axs[1, 1].legend(loc="upper right", frameon=True, fontsize=10)
    axs[1, 1].set_ylim(0, max(float(np.max(contacts)) * 1.1, 5))

    plt.suptitle(f"Molecular Dynamics Metrics: {display_title}", fontsize=18, fontweight="bold", y=1.02)
    plt.tight_layout()
    out_img_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_img_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved figure to {out_img_path}\n")
    return True


def main():
    parser = argparse.ArgumentParser(description="Generate per-system MD metric panels.")
    parser.add_argument("--results-dir", default="runs")
    parser.add_argument("--out-dir", default="paper_figures")
    parser.add_argument("--ns", type=float, default=100.0)
    parser.add_argument("--stride", type=int, default=10)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    out_dir = Path(args.out_dir)
    for sys_name, display_title, out_file in DEFAULT_SYSTEMS:
        plot_system_metrics(results_dir / sys_name, display_title, out_dir / out_file, args.ns, args.stride)

    print("Batch figure generation completed.")


if __name__ == "__main__":
    main()
