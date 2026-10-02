import argparse
import subprocess
import os
import shutil
from rdkit import Chem
from rdkit.Chem import AllChem

def reconstruct_pose_topology(ref_sdf, pose_sdf, out_sdf):
    print(f"[+] Reconstruyendo topología all-atom rigurosa...")
    ref_mol = next(Chem.SDMolSupplier(ref_sdf, removeHs=False))
    pose_mol = next(Chem.SDMolSupplier(pose_sdf, removeHs=False))
    
    if ref_mol is None or pose_mol is None:
        raise ValueError("Error al leer los SDF de referencia o docking.")

    # 1. Aislar esqueleto de átomos pesados de la pose de docking
    pose_heavy = Chem.RemoveHs(pose_mol)
    
    # 2. Encontrar el isomorfismo de subgrafos (mapeo topológico)
    match_indices = ref_mol.GetSubstructMatch(pose_heavy)
    if not match_indices:
        raise ValueError("Discordancia topológica: No se pudo mapear la pose sobre la referencia.")

    # 3. Transferencia de coordenadas estricta
    ref_conf = ref_mol.GetConformer()
    pose_conf = pose_heavy.GetConformer()
    
    for pose_idx, ref_idx in enumerate(match_indices):
        pos = pose_conf.GetAtomPosition(pose_idx)
        ref_conf.SetAtomPosition(ref_idx, pos)

    # 4. Relajación de hidrógenos explícitos (Force Field MMFF94)
    # Se genera el campo de fuerza y se anclan los átomos pesados (grados de libertad = 0)
    mp = AllChem.MMFFGetMoleculeProperties(ref_mol)
    ff = AllChem.MMFFGetMoleculeForceField(ref_mol, mp)
    
    for ref_idx in match_indices:
        ff.AddFixedPoint(ref_idx) # Fricción atómica infinita para átomos pesados
        
    print("[*] Minimizando geometría de hidrógenos...")
    ff.Minimize(maxIts=1000)

    # 5. Exportación del ligando híbrido saneado
    with Chem.SDWriter(out_sdf) as writer:
        writer.write(ref_mol)
    print(f"[✓] Pose reconstruida y saneada exportada en: {out_sdf}")
    return out_sdf

def run_acpype(input_sdf, basename, net_charge, out_dir):
    print(f"\n[+] Generando topología GAFF2/AM1-BCC para: {input_sdf}")
    os.makedirs(out_dir, exist_ok=True)
    
    cmd = [
        "acpype",
        "-i", input_sdf,
        "-b", basename,
        "-c", "bcc",
        "-n", str(net_charge),
        "-a", "gaff2",
        "-o", "gmx"
    ]
    
    print(f"[*] Ejecutando: {' '.join(cmd)}")
    res = subprocess.run(cmd, capture_output=True, text=True)
    
    if res.returncode != 0:
        print("[!] Error en ACPYPE/Antechamber:")
        print(res.stderr if res.stderr else res.stdout)
        raise RuntimeError("Fallo al generar la topología del ligando.")
        
    acpype_folder = f"{basename}.acpype"
    if not os.path.exists(acpype_folder):
        raise FileNotFoundError(f"Carpeta de salida no encontrada: {acpype_folder}")
        
    itp_file = os.path.join(acpype_folder, f"{basename}_GMX.itp")
    gro_file = os.path.join(acpype_folder, f"{basename}_GMX.gro")
    posre_file = os.path.join(acpype_folder, f"posre_{basename}.itp")
    
    shutil.copy(itp_file, os.path.join(out_dir, "ligand.itp"))
    shutil.copy(gro_file, os.path.join(out_dir, "ligand.gro"))
    if os.path.exists(posre_file):
        shutil.copy(posre_file, os.path.join(out_dir, "posre_ligand.itp"))
        
    shutil.rmtree(acpype_folder)
    print(f"[✓] Topología GROMACS generada exitosamente en: {out_dir}/")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reconstrucción topológica rigurosa y parametrización con ACPYPE")
    parser.add_argument("--ref-sdf", required=True, help="SDF de referencia (topología all-atom perfecta)")
    parser.add_argument("--pose-sdf", required=True, help="SDF salido del docking (coordenadas pesadas de alta confianza)")
    parser.add_argument("--basename", default="BEN", help="Nombre del residuo para GROMACS")
    parser.add_argument("--net-charge", type=int, default=1, help="Carga formal")
    parser.add_argument("--out-dir", default="results/topology", help="Directorio principal de salida")
    args = parser.parse_args()

    # Archivo temporal/definitivo para el ligando saneado
    sanitized_sdf = os.path.join(args.out_dir, "ligand_sanitized.sdf")
    os.makedirs(args.out_dir, exist_ok=True)

    # 1. Pipeline de reconstrucción topológica
    reconstruct_pose_topology(args.ref_sdf, args.pose_sdf, sanitized_sdf)
    
    # 2. Pipeline de parametrización mecano-cuántica
    run_acpype(sanitized_sdf, args.basename, args.net_charge, args.out_dir)
