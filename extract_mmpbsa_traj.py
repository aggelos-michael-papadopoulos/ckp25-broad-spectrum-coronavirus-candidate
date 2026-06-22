#!/usr/bin/env python3
import argparse
import os
import mdtraj as md
from pathlib import Path

DEFAULT_SYSTEMS = [
    "MERS_B_Spike_md",
    "SARS-CoV-1_E_spike_md",
    "hCoV_229E_E_spike_md",
    "hCoV_HKU1_A_spike_md",
    "hCoV_NL63_E_Spike_md"
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

def process_system(sys_dir: Path, stride: int, frames: int):
    traj_file = sys_dir / "md.xtc"
    top_file = sys_dir / "npt.gro"

    if not traj_file.exists() or not top_file.exists():
        print(f"  Missing md.xtc or npt.gro in {sys_dir}, skipping.")
        return

    try:
        traj = md.load(str(traj_file), top=str(top_file), stride=stride)

        traj.image_molecules(inplace=True)
        traj.center_coordinates()

        top = traj.topology
        lig_name = detect_lig_resname(top)
        sel = top.select(f"protein or resname {lig_name}")
        traj_stripped = traj.atom_slice(sel)

        if frames > 0 and len(traj_stripped) > frames:
            traj_stripped = traj_stripped[-frames:]

        out_nc = sys_dir / "complex.nc"
        traj_stripped.save_netcdf(str(out_nc))
        print(f"  Saved {len(traj_stripped)} frames to {out_nc.name}")

    except Exception as e:
        print(f"  Error: {str(e)}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract protein-ligand NetCDF trajectories for AmberTools MMPBSA.py."
    )
    parser.add_argument("--results-dir", default="runs")
    parser.add_argument("--systems", nargs="*", default=DEFAULT_SYSTEMS)
    parser.add_argument("--stride", type=int, default=50)
    parser.add_argument(
        "--frames",
        type=int,
        default=100,
        help="Keep the last N frames after striding; use 0 to keep all.",
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    for sys_name in args.systems:
        print(f"Processing {sys_name}...")
        process_system(results_dir / sys_name, args.stride, args.frames)


if __name__ == "__main__":
    main()
