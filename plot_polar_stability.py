#!/usr/bin/env python3
import argparse
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def get_polar_contacts_data(results_dir, sys_folder, standard_amino_acids, polar_elements, stride):
    sys_dir = results_dir / sys_folder
    traj_file = sys_dir / "md.xtc"
    top_file = sys_dir / "npt.gro"
    
    if not traj_file.exists() or not top_file.exists():
        return None
        
    traj = md.load(str(traj_file), top=str(top_file), stride=stride)
    traj.image_molecules(inplace=True)
    top = traj.topology
    
    # Detect ligand
    lig_resnames = set()
    for residue in top.residues:
        if residue.name not in standard_amino_acids:
            if residue.name not in {'SOL', 'WAT', 'HOH', 'NA', 'CL', 'K', 'Na', 'Cl', 'ION'}:
                lig_resnames.add(residue.name)
    lig_resname = list(lig_resnames)[0] if lig_resnames else "MOL"
    
    # Identify polar atoms
    lig_polar = np.array([a.index for a in top.atoms 
                           if a.residue.name == lig_resname and a.element.symbol in polar_elements])
    prot_polar = np.array([a.index for a in top.atoms 
                            if a.residue.is_protein and a.element.symbol in polar_elements])
    
    if len(lig_polar) == 0 or len(prot_polar) == 0:
        return np.zeros(len(traj))
        
    pairs = np.array([[lp, pp] for lp in lig_polar for pp in prot_polar])
    distances = md.compute_distances(traj, pairs)  # shape: [n_frames, n_pairs]
    
    # Polar contacts count per frame (distance <= 3.5 A / 0.35 nm)
    contacts_per_frame = np.sum(distances < 0.35, axis=1)
    return contacts_per_frame

def main():
    parser = argparse.ArgumentParser(description="Plot polar contact distributions and time traces.")
    parser.add_argument("--results-dir", default="runs")
    parser.add_argument("--out-dir", default="paper_figures")
    parser.add_argument("--ns", type=float, default=100.0)
    parser.add_argument("--stride", type=int, default=10)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    standard_amino_acids = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'CYX'}
    polar_elements = {"N", "O", "S"}
    
    systems = [
        {"name": "MERS-CoV (4KR0)", "folder": "MERS_B_Spike_md", "color": "#800080", "key": "MERS"},
        {"name": "HCoV-NL63 (3KBH)", "folder": "hCoV_NL63_E_Spike_md", "color": "#FF8C00", "key": "NL63"},
        {"name": "SARS-CoV-1 (2AJF)", "folder": "SARS-CoV-1_E_spike_md", "color": "#1F77B4", "key": "SARS1"},
        {"name": "HCoV-HKU1 (8Y7Y)", "folder": "hCoV_HKU1_A_spike_md", "color": "#008080", "key": "HKU1"},
        {"name": "HCoV-229E (6ATK)", "folder": "hCoV_229E_E_spike_md", "color": "#2CA02C", "key": "229E"}
    ]
    
    # Extract data
    data = {}
    for sys in systems:
        print(f"Extracting polar contacts for {sys['name']}...")
        contacts = get_polar_contacts_data(results_dir, sys["folder"], standard_amino_acids, polar_elements, args.stride)
        if contacts is not None:
            data[sys["key"]] = contacts

    missing = [sys["key"] for sys in systems if sys["key"] not in data]
    if missing:
        raise SystemExit(f"Missing trajectory data for: {', '.join(missing)}")
            
    # Plotting setup
    plt.style.use('default')
    plt.rcParams.update({
        'font.size': 12,
        'axes.linewidth': 1.2,
        'axes.labelsize': 13,
        'axes.titlesize': 14,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'xtick.major.width': 1.2,
        'ytick.major.width': 1.2,
        'lines.linewidth': 1.5
    })
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
    
    # -------------------------------------------------------------
    # Panel A: Violin Plot of Polar Contacts Distribution
    # -------------------------------------------------------------
    violin_data = [data[sys["key"]] for sys in systems]
    labels = [sys["name"].split(" ")[0] for sys in systems]  # Short labels: MERS-CoV, HCoV-NL63, etc.
    colors = [sys["color"] for sys in systems]
    
    # Create violin plot with means and extrema
    parts = ax1.violinplot(violin_data, showmeans=True, showmedians=False, showextrema=True)
    
    # Style the violin bodies
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_edgecolor('black')
        pc.set_alpha(0.6)
        
    # Style the mean and extrema lines
    parts['cmeans'].set_edgecolor('black')
    parts['cmeans'].set_linewidth(2.0)
    parts['cmins'].set_edgecolor('black')
    parts['cmins'].set_linewidth(1.0)
    parts['cmaxes'].set_edgecolor('black')
    parts['cmaxes'].set_linewidth(1.0)
    parts['cbars'].set_edgecolor('black')
    parts['cbars'].set_linewidth(1.0)
                
    ax1.set_xticks(np.arange(1, len(labels) + 1))
    ax1.set_xticklabels(labels)
    ax1.set_title("A — Polar Contacts Distribution Comparison", fontweight="bold", pad=12)
    ax1.set_ylabel(r"Number of Polar Contacts (d $\leq$ 3.5 Å)")
    ax1.set_xlabel("Target System")
    ax1.grid(axis='y', linestyle='--', alpha=0.5)
    
    # -------------------------------------------------------------
    # Panel B: Time-series of Smoothed Polar Contacts
    # -------------------------------------------------------------
    time_ns = np.linspace(0, args.ns, len(violin_data[0]))
    window = 25  # 5 ns smoothing window (since stride=10, 500 frames total, 1 frame = 0.2 ns)
    
    for sys in systems:
        raw_trace = data[sys["key"]]
        # Compute moving average
        smoothed = np.convolve(raw_trace, np.ones(window)/window, mode='valid')
        smoothed_time = time_ns[window-1:]
        
        ax2.plot(smoothed_time, smoothed, color=sys["color"], label=f"{sys['name'].split(' ')[0]} (mean {np.mean(raw_trace):.2f})", linewidth=2.0)
        
    ax2.set_title("B — Dynamic Polar Contact Stability (5-ns Smooth)", fontweight="bold", pad=12)
    ax2.set_ylabel(r"Number of Polar Contacts (d $\leq$ 3.5 Å)")
    ax2.set_xlabel("Time (ns)")
    ax2.set_xlim(0, args.ns)
    ax2.legend(loc="upper right", fontsize=10, frameon=True)
    ax2.grid(linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "md_polar_contacts_stability.jpg"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    
    print(f"Successfully generated and saved high-impact polar contacts stability plot to {out_path}")

if __name__ == "__main__":
    main()
