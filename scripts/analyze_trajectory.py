import argparse
import os
import numpy as np
import warnings
import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis as mda
from MDAnalysis.analysis import rms

def run_analysis(tpr_file, xtc_file, out_dir):
    print(f"[+] Cargando Universo con PBC pre-corregido por GROMACS: {tpr_file}")
    warnings.filterwarnings('ignore', category=UserWarning)

    u = mda.Universe(tpr_file, xtc_file)
    ref = mda.Universe(tpr_file)

    os.makedirs(out_dir, exist_ok=True)
    sns.set_context("paper", font_scale=1.2)
    sns.set_style("whitegrid")

    # --- 1. RESOLUCIÓN TOPOLÓGICA IN SILICO ---
    print("[*] Mapeando topología dinámica...")
    ligand_res = u.select_atoms("not protein and not resname SOL WATER WAT NA CL ION").residues
    if len(ligand_res) == 0:
         raise ValueError("Error topológico: Ligando no detectado.")
    ligand_name = ligand_res[0].resname

    asp_candidates = u.select_atoms(f"resname ASP and around 5.0 resname {ligand_name}").residues
    if len(asp_candidates) == 0:
         raise ValueError("Error topológico: Aspartato no detectado.")
    s1_asp = asp_candidates[0]

    print(f"    -> Ligando resuelto como: {ligand_name}")
    print(f"    -> Anclaje S1 (Asp) resuelto en el resid: {s1_asp.resid}")

    # --- 2. EXTRACCIÓN DE MÉTRICAS (RMSD) ---
    print("[*] Computando desviación cuadrática media (RMSD)...")
    rmsd_prot = rms.RMSD(u, ref, select='backbone', groupselections=['backbone'])
    rmsd_prot.run()

    pocket_sel = f"backbone and (byres (around 6.5 resname {ligand_name}))"
    # Fallback riguroso: Alinear al bolsillo y medir TODOS los átomos del ligando si no se detectan hidrógenos
    lig_sel = u.select_atoms(f"resname {ligand_name} and not type H")
    if len(lig_sel) == 0:
        ligand_rmsd_sel = f"resname {ligand_name}"
    else:
        ligand_rmsd_sel = f"resname {ligand_name} and not type H"

    rmsd_lig = rms.RMSD(u, ref, select=pocket_sel, groupselections=[ligand_rmsd_sel])
    rmsd_lig.run()

    # --- 3. INTEGRACIÓN DE DISTANCIA ROBUSTA ---
    print("[*] Integrando vectores de distancia (Center of Geometry)...")
    asp_carboxyl = u.select_atoms(f"resname ASP and resid {s1_asp.resid} and (name OD1 or name OD2)")

    if len(asp_carboxyl) == 0: # Último recurso: usar el C-alpha del Asp
        asp_carboxyl = u.select_atoms(f"resname ASP and resid {s1_asp.resid} and name CA")

    # Inferencia garantizada: Usar todos los átomos del residuo ligando
    lig_n = u.select_atoms(f"resname {ligand_name}")

    n_frames = len(u.trajectory)
    times = np.zeros(n_frames)
    distances = np.zeros(n_frames)

    for i, ts in enumerate(u.trajectory):
        times[i] = ts.time / 1000.0
        # Reemplazo táctico: center_of_geometry elude el vector nulo de masas (evita el warning numérico)
        distances[i] = np.linalg.norm(asp_carboxyl.center_of_geometry() - lig_n.center_of_geometry())

    # --- 4. RENDERIZADO ANALÍTICO ---
    print("[+] Exportando panel de correlación...")
    fig, axes = plt.subplots(3, 1, figsize=(8, 12), sharex=True)

    t_ns = rmsd_prot.results.rmsd[:, 1] / 1000.0

    axes[0].plot(t_ns, rmsd_prot.results.rmsd[:, 2] / 10.0, color='#1f77b4', lw=1.5)
    axes[0].set_ylabel(r"RMSD Backbone (nm)")
    axes[0].set_title("Estabilidad Conformacional Global (Tripsina)")

    axes[1].plot(t_ns, rmsd_lig.results.rmsd[:, 3] / 10.0, color='#ff7f0e', lw=1.5)
    axes[1].set_ylabel(r"RMSD Ligando (nm)")
    axes[1].set_title(f"Estabilidad de Pose ({ligand_name} alineado al Bolsillo S1)")

    axes[2].plot(times, distances / 10.0, color='#2ca02c', lw=1.5)
    axes[2].set_xlabel("Tiempo de Simulación (ns)")
    axes[2].set_ylabel(r"Distancia (nm)")
    axes[2].set_title(f"Distancia: Asp{s1_asp.resid} \(\leftrightarrow\) Centroide de {ligand_name}")

    plt.tight_layout()
    plot_path = os.path.join(out_dir, "md_stability_panel.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')

    print(f"[✓] Análisis completado sin fricción numérica. Gráfico en: {plot_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tpr", default="results/md_run/md_ref.pdb")
    parser.add_argument("--xtc", default="results/md_run/md_nopbc.xtc")
    parser.add_argument("--out-dir", default="results/analysis")
    args = parser.parse_args()

    run_analysis(args.tpr, args.xtc, args.out_dir)
