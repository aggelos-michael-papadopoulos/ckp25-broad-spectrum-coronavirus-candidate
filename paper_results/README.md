# Processed Paper Results

This directory contains lightweight numerical source data exported from the local MD/MM-GBSA analyses. Raw trajectories, binary GROMACS files, checkpoints, raw docking logs, and full intermediate outputs are not included.

MD trajectory-derived tables were exported with stride=10. Residue contact heatmap tables use 50 time bins and a 4 angstrom protein-ligand heavy-atom cutoff. Close-contact and polar-contact time series use a 3.5 angstrom cutoff.

Subdirectories:

- `docking/`: processed selected-pose docking score table only; raw docking logs and docking-box logs are not included.
- `md_processed/`: numerical source data for RMSD, RMSF, radius of gyration, close contacts, polar contacts, and residue contact maps.
- `mmgbsa/`: final MM-GBSA summary terms and per-frame delta energies.
- `figures/`: placeholder for regenerated figures; publication figures can be recreated from the workflow scripts.
