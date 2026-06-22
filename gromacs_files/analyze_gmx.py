#!/usr/bin/env python3
"""
analyze_gmx.py — Post-simulation analysis for GROMACS MD.

Produces four plots matching the OpenMM pipeline:
  plot_A_hbonds.png  — Protein–ligand hydrogen bonds over time
  plot_B_rmsd.png    — Protein Cα + ligand RMSD vs time
  plot_C_rmsf.png    — Per-residue Cα RMSF
  plot_D_rg.png      — Protein radius of gyration vs time

Usage:
  python analyze_gmx.py [--outdir results/run]
  python analyze_gmx.py --traj md.xtc --top md.tpr --outdir .

Requires: mdtraj, numpy, matplotlib
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import mdtraj as md
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── Helpers ─────────────────────────────────────────────────────────────────

def load_trajectory(traj_path: Path, top_path: Path, stride: int = 1) -> md.Trajectory:
    print(f"Loading trajectory: {traj_path}  (stride={stride})")
    t = md.load(str(traj_path), top=str(top_path), stride=stride)
    # Remove PBC artefacts and centre on protein
    t.image_molecules(inplace=True)
    t.center_coordinates()
    return t


def time_ns(traj: md.Trajectory) -> np.ndarray:
    """Return time axis in nanoseconds."""
    return traj.time / 1000.0


# ── Plot A — Hydrogen bonds ──────────────────────────────────────────────────

def plot_hbonds(traj: md.Trajectory, lig_resname: str, outdir: Path):
    print("Computing protein–ligand H-bonds + polar contacts...")

    top = traj.topology
    lig_atoms  = [a.index for a in top.atoms if a.residue.name == lig_resname]
    prot_atoms = set(a.index for a in top.atoms if a.residue.is_protein)

    if not lig_atoms:
        print(f"  WARNING: no residue named '{lig_resname}' found — skipping H-bond plot.")
        return

    lig_set = set(lig_atoms)

    # Polar heavy atoms (H-bond capable) for contact metric
    polar_elements = {"N", "O", "S"}
    lig_polar  = np.array([a.index for a in top.atoms
                           if a.residue.name == lig_resname and a.element.symbol in polar_elements])
    prot_polar = np.array([a.index for a in top.atoms
                           if a.residue.is_protein and a.element.symbol in polar_elements])

    hbonds_per_frame  = []
    contacts_per_frame = []

    for frame in traj:
        # Strict H-bonds (baker_hubbard geometry)
        hb = md.baker_hubbard(frame, periodic=False)
        count = sum(
            1 for d, h, a in hb
            if (d in prot_atoms and a in lig_set) or
               (d in lig_set    and a in prot_atoms)
        )
        hbonds_per_frame.append(count)

        # Relaxed polar contacts: N/O/S pairs within 3.5 Å
        if len(lig_polar) > 0 and len(prot_polar) > 0:
            pairs = np.array([[lp, pp] for lp in lig_polar for pp in prot_polar])
            dists = md.compute_distances(frame, pairs)[0]
            contacts_per_frame.append(int((dists < 0.35).sum()))  # 0.35 nm = 3.5 Å
        else:
            contacts_per_frame.append(0)

    hb_arr  = np.array(hbonds_per_frame)
    ct_arr  = np.array(contacts_per_frame)
    t       = time_ns(traj)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(t, ct_arr, lw=0.8, color="steelblue", alpha=0.5,
            label=f"Polar contacts ≤3.5 Å  (mean {ct_arr.mean():.2f})")
    ax.plot(t, hb_arr, lw=0.8, color="darkorange", alpha=0.8,
            label=f"H-bonds (strict)  (mean {hb_arr.mean():.2f})")
    ax.axhline(ct_arr.mean(), color="steelblue", lw=1.2, ls="--", alpha=0.7)
    ax.axhline(hb_arr.mean(), color="crimson",   lw=1.5, ls="--",
               label=f"Mean H-bonds = {hb_arr.mean():.2f}")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("# interactions")
    ax.set_title("A — Protein–Ligand Hydrogen Bonds & Polar Contacts")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = outdir / "plot_A_hbonds.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  → {out}  (H-bonds mean {hb_arr.mean():.2f}, polar contacts mean {ct_arr.mean():.2f})")


# ── Plot B — RMSD ────────────────────────────────────────────────────────────

def plot_rmsd(traj: md.Trajectory, lig_resname: str, outdir: Path):
    print("Computing RMSD...")
    top = traj.topology

    # Protein Cα
    ca_idx = top.select("name CA")
    ref = traj[0]
    traj.superpose(ref, atom_indices=ca_idx)
    rmsd_prot = md.rmsd(traj, ref, 0, atom_indices=ca_idx, precentered=True) * 10  # nm → Å

    # Ligand heavy atoms
    lig_idx = top.select(f"resname {lig_resname} and not element H")
    if len(lig_idx) > 0:
        rmsd_lig = md.rmsd(traj, ref, 0, atom_indices=lig_idx, precentered=True) * 10
    else:
        rmsd_lig = None
        print(f"  WARNING: ligand '{lig_resname}' not found for RMSD.")

    t = time_ns(traj)

    # Cap y-axis at 99th percentile to avoid spikes dominating the scale
    all_vals = rmsd_prot if rmsd_lig is None else np.concatenate([rmsd_prot, rmsd_lig])
    ymax = np.percentile(all_vals, 99) * 1.3
    spike_val = all_vals.max()

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(t, rmsd_prot, lw=0.8, label="Protein Cα", color="steelblue")
    if rmsd_lig is not None:
        ax.plot(t, rmsd_lig, lw=0.8, label=f"Ligand ({lig_resname})", color="darkorange")
    ax.set_ylim(0, ymax)
    # Annotate if there are spikes above the cap
    if spike_val > ymax:
        ax.annotate(f"max spike: {spike_val:.1f} Å (clipped)",
                    xy=(0.99, 0.97), xycoords="axes fraction",
                    ha="right", va="top", fontsize=8, color="grey",
                    style="italic")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title("B — RMSD from Starting Structure")
    ax.legend()
    fig.tight_layout()
    out = outdir / "plot_B_rmsd.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  → {out}")


# ── Plot C — RMSF ────────────────────────────────────────────────────────────

def plot_rmsf(traj: md.Trajectory, outdir: Path):
    print("Computing per-residue RMSF...")
    top = traj.topology
    ca_idx = top.select("name CA")
    rmsf = md.rmsf(traj, traj, 0, atom_indices=ca_idx) * 10  # Å

    res_ids = [top.atom(i).residue.resSeq for i in ca_idx]

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(res_ids, rmsf, lw=0.8, color="mediumseagreen")
    ax.fill_between(res_ids, rmsf, alpha=0.3, color="mediumseagreen")
    ax.set_xlabel("Residue number")
    ax.set_ylabel("RMSF (Å)")
    ax.set_title("C — Per-Residue Cα RMSF")
    fig.tight_layout()
    out = outdir / "plot_C_rmsf.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  → {out}")


# ── Plot D — Radius of gyration ──────────────────────────────────────────────

def plot_rg(traj: md.Trajectory, outdir: Path):
    print("Computing radius of gyration...")
    top    = traj.topology
    prot   = top.select("protein")
    rg     = md.compute_rg(traj.atom_slice(prot)) * 10  # Å
    t      = time_ns(traj)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(t, rg, lw=0.8, color="mediumpurple")
    ax.axhline(rg.mean(), color="crimson", lw=1.5, ls="--",
               label=f"Mean = {rg.mean():.2f} Å")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Rg (Å)")
    ax.set_title("D — Protein Radius of Gyration")
    ax.legend()
    fig.tight_layout()
    out = outdir / "plot_D_rg.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  → {out}")


# ── Auto-detect ligand residue name ─────────────────────────────────────────

def detect_lig_resname(traj: md.Trajectory) -> str:
    top = traj.topology
    non_std = set()
    for res in top.residues:
        if not res.is_protein and res.name not in ("HOH", "SOL", "WAT", "NA", "CL", "K"):
            non_std.add(res.name)
    if len(non_std) == 1:
        name = non_std.pop()
        print(f"  Auto-detected ligand residue name: {name}")
        return name
    if non_std:
        name = sorted(non_std)[0]
        print(f"  Multiple non-standard residues {non_std}; using first: {name}")
        return name
    print("  WARNING: no ligand residue detected.")
    return "LIG"


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Analyse GROMACS MD trajectory")
    parser.add_argument("--outdir", default=".", help="Directory with md.xtc / md.tpr (and output)")
    parser.add_argument("--traj",   default=None, help="Trajectory file (default: <outdir>/md.xtc)")
    parser.add_argument("--top",    default=None, help="Topology file   (default: <outdir>/md.tpr)")
    parser.add_argument("--lig",    default=None, help="Ligand residue name (auto-detected if omitted)")
    parser.add_argument("--stride", type=int, default=1,
                        help="Load every Nth frame (default 1 = all). Use 10+ for large trajectories.")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    traj_p = Path(args.traj)  if args.traj else outdir / "md.xtc"
    top_p  = Path(args.top)   if args.top  else outdir / "npt.gro"

    for p in (traj_p, top_p):
        if not p.exists():
            sys.exit(f"File not found: {p}")

    traj = load_trajectory(traj_p, top_p, args.stride)

    lig_resname = args.lig if args.lig else detect_lig_resname(traj)

    plot_hbonds(traj, lig_resname, outdir)
    plot_rmsd(traj, lig_resname, outdir)
    plot_rmsf(traj, outdir)
    plot_rg(traj, outdir)

    print("\nAll plots saved.")


if __name__ == "__main__":
    main()
