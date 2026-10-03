# CADD Trypsin Benchmark: Automated Redocking & MD Pipeline

![Python](https://img.shields.io/badge/Python-3.10-blue.svg)
![GROMACS](https://img.shields.io/badge/GROMACS-2024+-green.svg)
![Snakemake](https://img.shields.io/badge/Snakemake-Pipeline-009485.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

A rigorous, reproducible, and fully automated end-to-end workflow for protein-ligand complex validation. This project addresses parameter variability and manual friction in Computer-Aided Drug Design (CADD) by implementing the Benzamidine-Trypsin system (PDB: 3PTB) as a methodological benchmark.

## 🎯 Motivation

Historically, CADD and Molecular Dynamics pipelines suffer from critical reproducibility issues:
1. **Structural Fragmentation:** Heavy reliance on manual GUI interactions and disjointed ad-hoc scripts.
2. **Force Field Collisions:** Inconsistencies in mixing rules when combining heterogeneous topologies.
3. **Naive Preprocessing:** Improper treatment of key structural ions and physiological protonation states during apoprotein preparation.

This repository encapsulates an *in silico* software architecture orchestrated via **Snakemake**, eliminating manual preprocessing and ensuring deterministic traceability from the raw crystallographic file to the final thermodynamic evaluation of the docking pose.

## 🧬 Scientific Architecture & Ecosystem

To ensure strict reproducibility and isolate the pipeline against dependency conflicts, the workflow employs a unified approach built on the GROMACS/AMBER standard:

* **Receptor Structure:** `PDBFixer` repairs missing atoms, retains stabilizing structural ions (Ca2+), and assigns biological protonation states at pH 7.4.
* **Ligand Topology:** `RDKit` and `Meeko` extract native coordinates, applying 3D-aware bond order assignment and defining formal cationic (+1) charges.
* **Docking Engine:** `GNINA` utilizes 3D Convolutional Neural Network scoring (`CNNscore`) coupled with a crystallographic autobox for rigorous redocking.
* **Quantum Mechanical Parameterization:** `ACPYPE` and `AmberTools` dynamically assign AM1-BCC charges and generate GROMACS-compatible (GAFF2) topologies.
* **Molecular Dynamics:** `GROMACS` simulates the complex in NVT and NPT ensembles utilizing rigorous thermal coupling (`V-rescale`) and pressure coupling (`c-rescale` / `Parrinello-Rahman`).
* **In-Memory Analysis:** Direct trajectory processing on un-pbc `.xtc` binaries is performed via `MDAnalysis` to compute structural RMSD and dynamic anchoring distances.

## 📁 Project Structure

```text
cadd-trypsin-benchmark/
├── bin/                       # Standalone binaries (e.g., GNINA)
├── config/
│   └── config.yaml            # Centralized system and biophysics parameters
├── data/
│   └── raw/                   # Raw crystallographic inputs (e.g., 3PTB.pdb)
├── scripts/                   # Atomic Python modules
│   ├── fix_protein.py
│   ├── prepare_ligand.py
│   ├── run_docking.py
│   ├── build_topology.py
│   ├── build_complex.py
│   ├── write_mdps.py
│   ├── run_md.py
│   └── analyze_trajectory.py
├── Snakefile                  # Snakemake DAG orchestrator
└── environment.yaml           # Conda/Mamba dependency tree
```

## ⚙️ Environment Setup (Linux / WSL2)

A native POSIX environment is required. The ecosystem is managed via `mamba` to efficiently resolve complex dependency trees (such as GROMACS, RDKit, and AmberTools).

### 1. Clone the Repository
```bash
git clone https://github.com/YOUR_USERNAME/cadd-trypsin-benchmark.git
cd cadd-trypsin-benchmark
```

### 2. Install the Ecosystem
```bash
mamba env create -f environment.yaml
mamba activate cadd-env
```
*(Note: If `mamba activate` fails, initialize your shell with `mamba shell init --shell bash` and source your `.bashrc`)*.

### 3. Fetch Standalone Binaries
Due to HPC version control standards, large compiled binaries like GNINA are excluded via `.gitignore` and must be fetched during initialization.
```bash
mkdir -p bin
wget https://github.com/gnina/gnina/releases/download/v1.1/gnina -O bin/gnina
chmod +x bin/gnina
```

## 🚀 Usage & Pipeline Execution

The entire Directed Acyclic Graph (DAG) is managed by Snakemake. Global parameters (system ID, target pH, simulation time in ns, and docking exhaustiveness) are centralized in `config/config.yaml`.

**1. Dry-Run (Validate DAG connectivity mathematically):**
```bash
snakemake --dry-run
```

**2. Execute Full Pipeline:**
Allocate logical cores according to your hardware limits.
```bash
snakemake --cores 6
```

## 📊 Outputs & Validation

Upon successful execution, intermediate data is securely routed to `results/`. The final rule generates a comprehensive analytical panel (`results/analysis/md_stability_panel.png`) detailing:
* **Global Conformational Stability:** Protein Backbone RMSD over time.
* **Pose Stability:** Ligand RMSD aligned strictly to the S1 binding pocket.
* **Thermodynamic Distance:** Center of Geometry integration between the Asp189 carboxyl anchor and the ligand centroid.
