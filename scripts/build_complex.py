import argparse
import subprocess
import os
import shutil

def run_cmd(cmd, work_dir=None, input_str=None):
    print(f"[*] Ejecutando en {work_dir or 'cwd'}: {' '.join(cmd)}")
    res = subprocess.run(cmd, input=input_str, capture_output=True, text=True, cwd=work_dir)
    if res.returncode != 0:
        print(f"[!] ERROR en: {cmd[0]}")
        print(res.stderr if res.stderr else res.stdout)
        raise RuntimeError(f"Fallo en la ejecución de {cmd[0]}")
    return res.stdout

def merge_gro_files(protein_gro, ligand_gro, output_gro):
    print("[+] Fusionando coordenadas de proteína y ligando...")
    with open(protein_gro, 'r') as f: prot_lines = f.readlines()
    with open(ligand_gro, 'r') as f: lig_lines = f.readlines()
    prot_count = int(prot_lines[1].strip())
    lig_count = int(lig_lines[1].strip())
    with open(output_gro, 'w') as f:
        f.write("Complex: Protein + Ligand\n")
        f.write(f"{prot_count + lig_count:>5}\n")
        f.writelines(prot_lines[2:-1])
        f.writelines(lig_lines[2:-1])
        f.write(prot_lines[-1])

def inject_ligand_topology(top_file, ligand_itp_name, ligand_resname):
    print(f"[+] Inyectando topología del ligando en {top_file}...")
    with open(top_file, 'r') as f: lines = f.readlines()
    with open(top_file, 'w') as f:
        for line in lines:
            f.write(line)
            if "forcefield.itp" in line:
                f.write(f'\n#include "{ligand_itp_name}"\n')
        f.write(f"{ligand_resname:<15} 1\n")

def create_ions_mdp(filepath):
    mdp_content = "integrator = steep\nemtol = 1000.0\nnsteps = 50000\nnstlist = 1\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\n"
    with open(filepath, 'w') as f: f.write(mdp_content)

def build_system(protein_pdb, ligand_dir, out_dir, forcefield="amber99sb-ildn", water="tip3p"):
    # Resolución estricta de rutas absolutas para el aislamiento
    abs_out_dir = os.path.abspath(out_dir)
    abs_prot_pdb = os.path.abspath(protein_pdb)
    abs_lig_dir = os.path.abspath(ligand_dir)
    os.makedirs(abs_out_dir, exist_ok=True)
    
    # Inyectar archivos del ligando en el directorio de ensamblaje
    shutil.copy(os.path.join(abs_lig_dir, "ligand.itp"), os.path.join(abs_out_dir, "ligand.itp"))
    lig_gro_src = os.path.join(abs_lig_dir, "ligand.gro")
    
    print("\n[1/5] Construyendo topología (pdb2gmx)...")
    # Al usar cwd=abs_out_dir, GROMACS genera rutas relativas impecables en el topol.top
    run_cmd([
        "gmx", "pdb2gmx",
        "-f", abs_prot_pdb,
        "-o", "protein.gro",
        "-p", "topol.top",
        "-i", "posre.itp",
        "-ff", forcefield,
        "-water", water,
        "-ignh"
    ], work_dir=abs_out_dir)
    
    print("\n[2/5] Ensamblando Complejo P-L...")
    prot_gro = os.path.join(abs_out_dir, "protein.gro")
    complex_gro = os.path.join(abs_out_dir, "complex.gro")
    top_file = os.path.join(abs_out_dir, "topol.top")
    
    merge_gro_files(prot_gro, lig_gro_src, complex_gro)
    inject_ligand_topology(top_file, "ligand.itp", "BEN")
    
    print("\n[3/5] Definiendo caja (editconf)...")
    run_cmd(["gmx", "editconf", "-f", "complex.gro", "-o", "complex_box.gro", "-d", "1.0", "-bt", "cubic"], work_dir=abs_out_dir)
    
    print("\n[4/5] Solvatando (solvate)...")
    run_cmd(["gmx", "solvate", "-cp", "complex_box.gro", "-cs", "spc216.gro", "-o", "complex_solv.gro", "-p", "topol.top"], work_dir=abs_out_dir)
    
    print("\n[5/5] Neutralizando (genion)...")
    ions_mdp = os.path.join(abs_out_dir, "ions.mdp")
    create_ions_mdp(ions_mdp)
    
    run_cmd(["gmx", "grompp", "-f", "ions.mdp", "-c", "complex_solv.gro", "-p", "topol.top", "-o", "ions.tpr", "-maxwarn", "1"], work_dir=abs_out_dir)
    run_cmd(["gmx", "genion", "-s", "ions.tpr", "-o", "system_ready.gro", "-p", "topol.top", "-pname", "NA", "-nname", "CL", "-neutral"], work_dir=abs_out_dir, input_str="SOL\n")
    
    print(f"\n[✓] SISTEMA ENSAMBLADO EN: {abs_out_dir}/")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein", required=True)
    parser.add_argument("--ligand-dir", required=True)
    parser.add_argument("--out-dir", default="results/complex")
    args = parser.parse_args()
    build_system(args.protein, args.ligand_dir, args.out_dir)
