# Curated Input Data

This folder contains the minimal starting structures needed to reproduce the molecular dynamics workflow from the manuscript.

## Contents

| System | Virus | PDB | Spike chain | Protein input | CKP-25 pose |
|---|---|---:|:---:|---|---|
| MERS | MERS-CoV | 4KR0 | B | `proteins/MERS_4KR0_B_spike.pdb` | `docked_poses/MERS_CKP25_top_pose.sdf` |
| SARS1 | SARS-CoV-1 | 2AJF | E | `proteins/SARS1_2AJF_E_spike.pdb` | `docked_poses/SARS1_CKP25_top_pose.sdf` |
| NL63 | HCoV-NL63 | 3KBH | E | `proteins/NL63_3KBH_E_spike.pdb` | `docked_poses/NL63_CKP25_top_pose.sdf` |
| 229E | HCoV-229E | 6ATK | E | `proteins/229E_6ATK_E_spike.pdb` | `docked_poses/229E_CKP25_top_pose.sdf` |
| HKU1 | HCoV-HKU1 | 8Y7Y | A | `proteins/HKU1_8Y7Y_A_spike.pdb` | `docked_poses/HKU1_CKP25_top_pose.sdf` |

`metadata.tsv` provides the same mapping in machine-readable form. `ckp25.smiles` contains the ligand SMILES used for CKP-25.
