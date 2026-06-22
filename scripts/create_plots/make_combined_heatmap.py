#!/usr/bin/env python3
import argparse
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def get_contact_data(results_dir, sys_folder, orig_pdb_path, standard_amino_acids, stride):
    sys_dir = results_dir / sys_folder
    traj_file = sys_dir / "md.xtc"
    top_file = sys_dir / "npt.gro"
    
    if not traj_file.exists() or not top_file.exists():
        return None
        
    traj = md.load(str(traj_file), top=str(top_file), stride=stride)
    traj.image_molecules(inplace=True)
    top = traj.topology
    
    # Load original PDB for sequence alignment
    orig = md.load(str(orig_pdb_path)).topology
    orig_res = [r for r in orig.residues if r.is_protein]
    
    # Detect ligand
    lig_resnames = set()
    for residue in top.residues:
        if residue.name not in standard_amino_acids:
            if residue.name not in {'SOL', 'WAT', 'HOH', 'NA', 'CL', 'K', 'Na', 'Cl', 'ION'}:
                lig_resnames.add(residue.name)
    lig_resname = list(lig_resnames)[0] if lig_resnames else "MOL"
    
    lig_atoms = top.select(f"resname {lig_resname} and not element H")
    protein_residues = [r for r in top.residues if r.is_protein]
    
    # Align simulation sequence to original PDB sequence
    sim_seq = [r.name for r in protein_residues]
    orig_seq = [r.name for r in orig_res]
    
    n, m = len(sim_seq), len(orig_seq)
    dp = np.zeros((n + 1, m + 1), dtype=int)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = (sim_seq[i-1] == orig_seq[j-1]) or \
                    (sim_seq[i-1] == 'HIS' and orig_seq[j-1] in ('HID', 'HIE', 'HIP', 'HIS')) or \
                    (sim_seq[i-1] == 'CYS' and orig_seq[j-1] in ('CYX', 'CYS'))
            if match:
                dp[i, j] = dp[i-1, j-1] + 1
            else:
                dp[i, j] = max(dp[i-1, j], dp[i, j-1])
                
    i, j = n, m
    alignment = []
    while i > 0 and j > 0:
        match = (sim_seq[i-1] == orig_seq[j-1]) or \
                (sim_seq[i-1] == 'HIS' and orig_seq[j-1] in ('HID', 'HIE', 'HIP', 'HIS')) or \
                (sim_seq[i-1] == 'CYS' and orig_seq[j-1] in ('CYX', 'CYS'))
        if match and dp[i, j] == dp[i-1, j-1] + 1:
            alignment.append((i-1, j-1))
            i -= 1
            j -= 1
        elif dp[i, j] == dp[i-1, j]:
            i -= 1
        else:
            j -= 1
            
    alignment.reverse()
    sim_to_orig_map = {}
    for sim_idx, orig_idx in alignment:
        sim_to_orig_map[sim_idx] = orig_res[orig_idx]
        
    prot_heavy_atoms = []
    atom_to_res_idx = {}
    for r_idx, r in enumerate(protein_residues):
        for a in r.atoms:
            if a.element.symbol != 'H':
                prot_heavy_atoms.append(a.index)
                atom_to_res_idx[a.index] = r_idx
                
    prot_heavy_atoms = np.array(prot_heavy_atoms)
    pairs = np.array([[la, pa] for la in lig_atoms for pa in prot_heavy_atoms])
    distances = md.compute_distances(traj, pairs)
    
    n_frames = len(traj)
    n_residues = len(protein_residues)
    
    contact_matrix = np.zeros((n_frames, n_residues), dtype=bool)
    for i, (la, pa) in enumerate(pairs):
        res_idx = atom_to_res_idx[pa]
        contact_matrix[:, res_idx] |= (distances[:, i] <= 0.4)
        
    freqs = np.mean(contact_matrix, axis=0) * 100
    sorted_indices = np.argsort(freqs)[::-1]
    top_n = 10
    top_res_indices = sorted_indices[:top_n]
    
    # Bin into 50 windows (2 ns each)
    n_bins = 50
    frames_per_bin = n_frames // n_bins
    binned_matrix = np.zeros((top_n, n_bins))
    for r_idx_in_top, res_idx in enumerate(top_res_indices):
        for b in range(n_bins):
            start_frame = b * frames_per_bin
            end_frame = (b + 1) * frames_per_bin
            binned_matrix[r_idx_in_top, b] = np.mean(contact_matrix[start_frame:end_frame, res_idx]) * 100
            
    y_labels = []
    for idx in top_res_indices:
        orig_r = sim_to_orig_map.get(idx)
        if orig_r:
            y_labels.append(f"{orig_r.name}{orig_r.resSeq}")
        else:
            r = protein_residues[idx]
            y_labels.append(f"{r.name}{r.resSeq}")
            
    return binned_matrix, y_labels

def main():
    parser = argparse.ArgumentParser(
        description="Generate residue contact heatmaps with original PDB residue numbering."
    )
    parser.add_argument("--results-dir", default="runs")
    parser.add_argument("--data-dir", default="Data")
    parser.add_argument("--out-dir", default="paper_figures")
    parser.add_argument("--stride", type=int, default=10)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    data_dir = Path(args.data_dir)
    standard_amino_acids = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'CYX'}
    
    systems = [
        {
            "name": "MERS-CoV (4KR0)",
            "folder": "MERS_B_Spike_md",
            "orig_pdb": data_dir / "proteins" / "MERS_4KR0_B_spike.pdb",
            "key": "MERS"
        },
        {
            "name": "HCoV-NL63 (3KBH)",
            "folder": "hCoV_NL63_E_Spike_md",
            "orig_pdb": data_dir / "proteins" / "NL63_3KBH_E_spike.pdb",
            "key": "NL63"
        },
        {
            "name": "SARS-CoV-1 (2AJF)",
            "folder": "SARS-CoV-1_E_spike_md",
            "orig_pdb": data_dir / "proteins" / "SARS1_2AJF_E_spike.pdb",
            "key": "SARS1"
        },
        {
            "name": "HCoV-HKU1 (8Y7Y)",
            "folder": "hCoV_HKU1_A_spike_md",
            "orig_pdb": data_dir / "proteins" / "HKU1_8Y7Y_A_spike.pdb",
            "key": "HKU1"
        },
        {
            "name": "HCoV-229E (6ATK)",
            "folder": "hCoV_229E_E_spike_md",
            "orig_pdb": data_dir / "proteins" / "229E_6ATK_E_spike.pdb",
            "key": "229E"
        }
    ]
    
    data = {}
    for sys in systems:
        print(f"Extracting data for {sys['name']}...")
        res = get_contact_data(results_dir, sys["folder"], sys["orig_pdb"], standard_amino_acids, args.stride)
        if res:
            data[sys["key"]] = res

    missing = [sys["key"] for sys in systems if sys["key"] not in data]
    if missing:
        raise SystemExit(f"Missing trajectory/contact data for: {', '.join(missing)}")
            
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Generate 4-panel combined figure (MERS, NL63, SARS1, 229E)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True)
    targets_4 = ["MERS", "NL63", "SARS1", "229E"]
    names_4 = {
        "MERS": "MERS-CoV (4KR0)",
        "NL63": "HCoV-NL63 (3KBH)",
        "SARS1": "SARS-CoV-1 (2AJF)",
        "229E": "HCoV-229E (6ATK)"
    }
    
    for idx, key in enumerate(targets_4):
        ax = axes[idx // 2, idx % 2]
        binned, y_labels = data[key]
        
        im = ax.imshow(binned, aspect='auto', cmap='RdYlBu_r', origin='lower', extent=[0, 100, -0.5, 9.5])
        ax.set_title(names_4[key], fontsize=13, fontweight="bold", pad=8)
        ax.set_yticks(np.arange(10))
        ax.set_yticklabels(y_labels, fontsize=10)
        
        if idx >= 2:
            ax.set_xlabel("Time (ns)", fontsize=11)
            
        cbar = fig.colorbar(im, ax=ax, pad=0.02)
        cbar.ax.tick_params(labelsize=9)
        if idx % 2 == 1:
            cbar.set_label("Contact Frequency (%)", fontsize=10)
            
    plt.tight_layout()
    out_path_4 = out_dir / "md_contact_heatmap_4panel.jpg"
    plt.savefig(out_path_4, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved 4-panel contact heatmap to {out_path_4}")
    
    # 2. Generate 5-panel combined figure (All 5 targets)
    fig, axes = plt.subplots(3, 2, figsize=(14, 14), sharex=True)
    targets_5 = ["MERS", "NL63", "SARS1", "HKU1", "229E"]
    names_5 = {
        "MERS": "MERS-CoV (4KR0)",
        "NL63": "HCoV-NL63 (3KBH)",
        "SARS1": "SARS-CoV-1 (2AJF)",
        "HKU1": "HCoV-HKU1 (8Y7Y)",
        "229E": "HCoV-229E (6ATK)"
    }
    
    for idx, key in enumerate(targets_5):
        ax = axes[idx // 2, idx % 2]
        binned, y_labels = data[key]
        
        im = ax.imshow(binned, aspect='auto', cmap='RdYlBu_r', origin='lower', extent=[0, 100, -0.5, 9.5])
        ax.set_title(names_5[key], fontsize=13, fontweight="bold", pad=8)
        ax.set_yticks(np.arange(10))
        ax.set_yticklabels(y_labels, fontsize=10)
        
        if idx >= 3:
            ax.set_xlabel("Time (ns)", fontsize=11)
            
        cbar = fig.colorbar(im, ax=ax, pad=0.02)
        cbar.ax.tick_params(labelsize=9)
        if idx % 2 == 1:
            cbar.set_label("Contact Frequency (%)", fontsize=10)
            
    # Remove the empty 6th subplot
    fig.delaxes(axes[2, 1])
    # Adjust last left subplot to look nice
    axes[2, 0].set_xlabel("Time (ns)", fontsize=11)
    
    plt.tight_layout()
    out_path_5 = out_dir / "md_contact_heatmap_5panel.jpg"
    plt.savefig(out_path_5, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved 5-panel contact heatmap to {out_path_5}")

if __name__ == "__main__":
    main()
