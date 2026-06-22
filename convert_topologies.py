#!/usr/bin/env python3
import argparse
import os
import sys
import parmed as pmd
from pathlib import Path

DEFAULT_SYSTEMS = [
    "MERS_B_Spike_md",
    "SARS-CoV-1_E_spike_md",
    "hCoV_229E_E_spike_md",
    "hCoV_HKU1_A_spike_md",
    "hCoV_NL63_E_Spike_md",
]

def fix_dihedrals(top):
    for dihedral in top.dihedrals:
        if dihedral.type is not None:
            if isinstance(dihedral.type, pmd.topologyobjects.DihedralTypeList):
                for dt in dihedral.type:
                    if dt.per < 1 or dt.per > 6:
                        dt.per = 1
            else:
                if dihedral.type.per < 1 or dihedral.type.per > 6:
                    dihedral.type.per = 1
    return top

def convert_topology(sys_dir):
    print(f"Processing {sys_dir}...")
    
    top_file = sys_dir / "topol.top"
    gro_file = sys_dir / "npt.gro"
    
    if not top_file.exists() or not gro_file.exists():
        print(f"  Missing topol.top or npt.gro in {sys_dir}")
        return False
        
    try:
        # 1. Complex
        complex_top = pmd.load_file(str(top_file), xyz=str(gro_file))
        complex_top.strip(':SOL,WAT,HOH,NA,CL,K,Na,Cl')
        complex_top = fix_dihedrals(complex_top)
        complex_top.box = None  # Remove PBC for MM-PBSA
        pmd.tools.actions.changeRadii(complex_top, 'mbondi2').execute()
        complex_out = sys_dir / "complex.prmtop"
        complex_top.save(str(complex_out), overwrite=True)
        print(f"  Saved {complex_out.name}")
        
        # Find ligand residue name
        lig_resnames = set()
        for res in complex_top.residues:
            if not res.name in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'CYX']:
                lig_resnames.add(res.name)
                
        lig_name = list(lig_resnames)[0] if lig_resnames else "LIG"
        print(f"  Detected ligand: {lig_name}")
            
        # 2. Receptor
        receptor_top = pmd.load_file(str(top_file), xyz=str(gro_file))
        receptor_top.strip(':SOL,WAT,HOH,NA,CL,K,Na,Cl')
        receptor_top.strip(f':{lig_name}')
        receptor_top = fix_dihedrals(receptor_top)
        receptor_top.box = None  # Remove PBC for MM-PBSA
        pmd.tools.actions.changeRadii(receptor_top, 'mbondi2').execute()
        receptor_out = sys_dir / "receptor.prmtop"
        receptor_top.save(str(receptor_out), overwrite=True)
        print(f"  Saved {receptor_out.name}")
        
        # 3. Ligand
        ligand_top = pmd.load_file(str(top_file), xyz=str(gro_file))
        ligand_top.strip(':SOL,WAT,HOH,NA,CL,K,Na,Cl')
        ligand_top.strip(f'!(:{lig_name})')
        ligand_top = fix_dihedrals(ligand_top)
        ligand_top.box = None  # Remove PBC for MM-PBSA
        pmd.tools.actions.changeRadii(ligand_top, 'mbondi2').execute()
        ligand_out = sys_dir / "ligand.prmtop"
        ligand_top.save(str(ligand_out), overwrite=True)
        print(f"  Saved {ligand_out.name}")
        
        return True
        
    except Exception as e:
        print(f"  Error processing {sys_dir}: {str(e)}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Convert GROMACS topologies to Amber prmtop files for MM-GBSA."
    )
    parser.add_argument(
        "--results-dir",
        default="runs",
        help="Directory containing one subfolder per simulated system.",
    )
    parser.add_argument(
        "--systems",
        nargs="*",
        default=DEFAULT_SYSTEMS,
        help="System subfolders to process.",
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    for sys_name in args.systems:
        convert_topology(results_dir / sys_name)

    print("\nTopology conversion completed!")


if __name__ == "__main__":
    main()
