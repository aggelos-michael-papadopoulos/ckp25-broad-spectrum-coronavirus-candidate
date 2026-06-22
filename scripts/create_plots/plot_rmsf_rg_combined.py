#!/usr/bin/env python3
import argparse
import os
import sys
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Plot combined RMSF and Rg stability panels.")
    parser.add_argument("--results-dir", default="runs")
    parser.add_argument("--out-dir", default="paper_figures")
    parser.add_argument("--ns", type=float, default=100.0)
    parser.add_argument("--stride", type=int, default=10)
    args = parser.parse_args()

    systems = [
        {"name": "MERS-CoV (4KR0)", "folder": "MERS_B_Spike_md", "color": "#9B59B6"},
        {"name": "HCoV-NL63 (3KBH)", "folder": "hCoV_NL63_E_Spike_md", "color": "#E67E22"},
        {"name": "SARS-CoV-1 (2AJF)", "folder": "SARS-CoV-1_E_spike_md", "color": "#3498DB"},
        {"name": "HCoV-HKU1 (8Y7Y)", "folder": "hCoV_HKU1_A_spike_md", "color": "#1ABC9C"},
        {"name": "HCoV-229E (6ATK)", "folder": "hCoV_229E_E_spike_md", "color": "#2ECC71"}
    ]
    
    results_dir = Path(args.results_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Dictionary to store calculated data
    data = {}
    
    for sys in systems:
        sys_folder = sys["folder"]
        sys_path = results_dir / sys_folder
        traj_file = sys_path / "md.xtc"
        top_file = sys_path / "npt.gro"
        
        if not traj_file.exists() or not top_file.exists():
            print(f"Error: files not found for {sys_folder}")
            continue
            
        print(f"Loading trajectory for {sys['name']}...")
        traj = md.load(str(traj_file), top=str(top_file), stride=args.stride)
        traj.image_molecules(inplace=True)
        traj.center_coordinates()
        
        top = traj.topology
        
        # 1. RMSF calculation (Protein CA)
        print(f"  Computing RMSF for {sys['name']}...")
        ca_idx = top.select("protein and name CA")
        traj.superpose(traj, 0, atom_indices=ca_idx)
        rmsf = md.rmsf(traj, traj, 0, atom_indices=ca_idx) * 10  # convert to Å
        
        # 2. Radius of gyration (Rg) calculation (Protein)
        print(f"  Computing Rg for {sys['name']}...")
        traj_prot = traj.atom_slice(top.select("protein"))
        rg = md.compute_rg(traj_prot) * 10  # convert to Å
        
        data[sys["folder"]] = {
            "rmsf": rmsf,
            "rg": rg,
            "time_ns": np.linspace(0, args.ns, len(traj)),
            "n_residues": len(ca_idx)
        }

    # Setup publication quality plotting parameters
    plt.style.use('default')
    plt.rcParams.update({
        'font.size': 12,
        'axes.linewidth': 1.5,
        'axes.labelsize': 13,
        'axes.titlesize': 14,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'xtick.major.width': 1.5,
        'ytick.major.width': 1.5,
        'lines.linewidth': 1.5
    })
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.8))
    
    # -------------------------------------------------------------
    # Panel A: RMSF (aligned by relative residue position)
    # -------------------------------------------------------------
    for sys in systems:
        folder = sys["folder"]
        if folder not in data:
            continue
        rmsf = data[folder]["rmsf"]
        x_vals = np.arange(1, len(rmsf) + 1)
        ax1.plot(x_vals, rmsf, color=sys["color"], label=sys["name"], alpha=0.9)
        
    ax1.set_title("A — Per-Residue Cα RMSF Profile", fontweight="bold", pad=12)
    ax1.set_xlabel("Relative Residue Index (N- to C-terminus)")
    ax1.set_ylabel("RMSF (Å)")
    ax1.set_xlim(1, 310)
    ax1.set_ylim(0, 6)
    ax1.set_xticks([1, 50, 100, 150, 200, 250, 300])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.grid(linestyle='--', alpha=0.4)
    ax1.legend(loc="upper right", frameon=True, fontsize=10)
    
    # -------------------------------------------------------------
    # Panel B: Radius of Gyration (Rg) over 100 ns
    # -------------------------------------------------------------
    for sys in systems:
        folder = sys["folder"]
        if folder not in data:
            continue
        rg = data[folder]["rg"]
        time_ns = data[folder]["time_ns"]
        ax2.plot(time_ns, rg, color=sys["color"], label=sys["name"], alpha=0.9)
        
    ax2.set_title("B — Protein Radius of Gyration (Rg) Stability", fontweight="bold", pad=12)
    ax2.set_xlabel("Time (ns)")
    ax2.set_ylabel("Rg (Å)")
    ax2.set_xlim(0, args.ns)
    ax2.set_ylim(14, 21)
    ax2.set_xticks([0, 20, 40, 60, 80, 100])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.grid(linestyle='--', alpha=0.4)
    ax2.legend(loc="upper right", frameon=True, fontsize=10)
    
    plt.tight_layout()
    
    out_path = out_dir / "md_rmsf_rg_stability.jpg"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    
    print(f"\nSuccessfully generated and saved combined RMSF/Rg plot to {out_path}!")

if __name__ == "__main__":
    main()
