import argparse
import subprocess
import os
from rdkit import Chem
from rdkit.Chem import rdMolAlign

def run_gnina_docking(gnina_bin, receptor_pdb, ligand_sdf, autobox_sdf, 
                      output_sdf, out_log, exhaustiveness=16, num_modes=9, seed=123, buffer_box=4.0):
    print("[+] Ejecutando redocking con GNINA...")
    
    cmd = [
        gnina_bin,
        "--receptor", receptor_pdb,
        "--ligand", ligand_sdf,          # Usamos SDF nativo para preservar la topología exacta
        "--autobox_ligand", autobox_sdf, # Usamos el mismo SDF para centrar la caja
        "--autobox_add", str(buffer_box),
        "--out", output_sdf,
        "--log", out_log,
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--seed", str(seed),
        "--cnn_scoring", "rescore",
        "--cpu", "4"
    ]
    
    res = subprocess.run(cmd, capture_output=True, text=True)
    
    # Defensa algorítmica: Validar generación real de poses
    if res.returncode != 0 or not os.path.exists(output_sdf) or os.path.getsize(output_sdf) == 0:
        print("[!] Log de Error Crítico de GNINA:")
        print(res.stderr if res.stderr else res.stdout)
        raise RuntimeError("GNINA falló o produjo una salida vacía.")
        
    print(f"[✓] Docking finalizado. Poses guardadas en: {output_sdf}")

def evaluate_best_pose(output_sdf, ref_sdf, best_pose_sdf):
    print("[+] Analizando poses y calculando RMSD cristalográfico...")
    suppl = Chem.SDMolSupplier(output_sdf, removeHs=False)
    poses = [mol for mol in suppl if mol is not None]
    
    if not poses:
        raise ValueError("El archivo SDF de GNINA no contiene poses válidas.")
        
    ref_suppl = Chem.SDMolSupplier(ref_sdf, removeHs=False)
    ref_mol = next(ref_suppl)
    best_pose = poses[0]
    
    # Extracción de métricas de afinidad y CNN
    affinity = best_pose.GetProp("minimizedAffinity") if best_pose.HasProp("minimizedAffinity") else "N/A"
    cnn_score = best_pose.GetProp("CNNscore") if best_pose.HasProp("CNNscore") else "N/A"
    cnn_affinity = best_pose.GetProp("CNNaffinity") if best_pose.HasProp("CNNaffinity") else "N/A"
    
    try:
        # RMSD riguroso exclusivo de átomos pesados
        pose_heavy = Chem.RemoveHs(best_pose)
        ref_heavy = Chem.RemoveHs(ref_mol)
        rmsd = rdMolAlign.GetBestRMS(pose_heavy, ref_heavy)
    except Exception as e:
        print(f"[!] Error analítico topológico al calcular RMSD: {e}")
        rmsd = float('inf')
    
    print("\n" + "=" * 50)
    print(" RESULTADOS DEL REDOCKING (Top Pose)")
    print("=" * 50)
    print(f" - Afinidad (Vina)   : {affinity} kcal/mol")
    print(f" - CNN Pose Score    : {cnn_score}")
    print(f" - CNN Affinity      : {cnn_affinity}")
    if rmsd != float('inf'):
        print(f" - RMSD vs Cristal   : {rmsd:.3f} Å")
    print("=" * 50 + "\n")
    
    if rmsd <= 2.0:
        print("[✓] Criterio de validación CADD SUPERADO (RMSD <= 2.0 Å)")
    else:
        print("[!] Advertencia: Desviación estructural significativa (> 2.0 Å)")

    with Chem.SDWriter(best_pose_sdf) as writer:
        writer.write(best_pose)
    print(f"[✓] Mejor pose exportada para Dinámica Molecular en: {best_pose_sdf}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Docking GNINA unificado por SDF")
    parser.add_argument("--gnina", default="bin/gnina", help="Ruta al binario de GNINA")
    parser.add_argument("--receptor", required=True, help="Receptor PDB")
    parser.add_argument("--ligand", required=True, help="Ligando de entrada (SDF nativo)")
    parser.add_argument("--autobox", required=True, help="SDF cristalográfico para definir caja")
    parser.add_argument("--out-sdf", default="results/docking_poses.sdf", help="Salida multimodelo SDF")
    parser.add_argument("--out-log", default="results/docking.log", help="Log de GNINA")
    parser.add_argument("--best-pose", default="results/best_pose.sdf", help="Salida de la mejor pose")
    parser.add_argument("--exhaustiveness", type=int, default=16)
    parser.add_argument("--seed", type=int, default=123)
    args = parser.parse_args()

    run_gnina_docking(
        args.gnina, args.receptor, args.ligand, args.autobox, 
        args.out_sdf, args.out_log, args.exhaustiveness, seed=args.seed
    )
    evaluate_best_pose(args.out_sdf, args.autobox, args.best_pose)
