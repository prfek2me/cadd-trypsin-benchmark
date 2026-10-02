import argparse
import subprocess
import sys
import os

def run_step(cmd, work_dir=None):
    """Ejecuta el comando confinado en el directorio de trabajo especificado."""
    print(f"\n[{'-'*50}]\n[+] Ejecutando en: {work_dir or 'cwd'}\n    {' '.join(cmd)}\n[{'-'*50}]")
    res = subprocess.run(cmd, cwd=work_dir)
    if res.returncode != 0:
        sys.exit(f"\n[!] Error crítico abortando el pipeline en: {cmd[0]} {cmd[1] if len(cmd)>1 else ''}")

def execute_md_pipeline(complex_dir, md_dir, mdp_dir):
    # 0. Resolución topológica rigurosa (Rutas absolutas)
    abs_topol = os.path.abspath(os.path.join(complex_dir, "topol.top"))
    abs_sys_gro = os.path.abspath(os.path.join(complex_dir, "system_ready.gro"))
    abs_mdp = os.path.abspath(mdp_dir)
    abs_md = os.path.abspath(md_dir)
    
    os.makedirs(abs_md, exist_ok=True)

    # Fase 0: Minimización de Energía (EM)
    print("\n[Fase 0/3] Minimización de Energía (EM)")
    run_step([
        "gmx", "grompp", 
        "-f", os.path.join(abs_mdp, "em.mdp"), 
        "-c", abs_sys_gro, 
        "-p", abs_topol, 
        "-o", "em.tpr",
        "-maxwarn", "1"
    ], work_dir=abs_md)
    run_step(["gmx", "mdrun", "-deffnm", "em", "-v"], work_dir=abs_md)

    em_gro = os.path.join(abs_md, "em.gro")
    if not os.path.exists(em_gro):
        raise FileNotFoundError(f"Fallo crítico: No se generó la estructura minimizada {em_gro}")

    # Fase 1: Equilibración Canónica (NVT - Volumen Constante)
    print("\n[Fase 1/3] Calentamiento (NVT)")
    run_step([
        "gmx", "grompp", 
        "-f", os.path.join(abs_mdp, "nvt.mdp"), 
        "-c", "em.gro", 
        "-r", "em.gro", 
        "-p", abs_topol, 
        "-o", "nvt.tpr",
        "-maxwarn", "1"
    ], work_dir=abs_md)
    run_step(["gmx", "mdrun", "-deffnm", "nvt", "-v"], work_dir=abs_md)

    # Fase 2: Equilibración Isobárica-Isotérmica (NPT - Presión Constante)
    print("\n[Fase 2/3] Relajación de Densidad (NPT)")
    run_step([
        "gmx", "grompp", 
        "-f", os.path.join(abs_mdp, "npt.mdp"), 
        "-c", "nvt.gro", 
        "-r", "nvt.gro", 
        "-t", "nvt.cpt", 
        "-p", abs_topol, 
        "-o", "npt.tpr",
        "-maxwarn", "1"
    ], work_dir=abs_md)
    run_step(["gmx", "mdrun", "-deffnm", "npt", "-v"], work_dir=abs_md)

    # Fase 3: Producción (MD)
    print("\n[Fase 3/3] Dinámica Molecular de Producción")
    run_step([
        "gmx", "grompp", 
        "-f", os.path.join(abs_mdp, "md.mdp"), 
        "-c", "npt.gro", 
        "-r", "npt.gro", 
        "-t", "npt.cpt", 
        "-p", abs_topol, 
        "-o", "md.tpr",
        "-maxwarn", "1"
    ], work_dir=abs_md)
    run_step(["gmx", "mdrun", "-deffnm", "md", "-v"], work_dir=abs_md)
    
    print("\n[✓] PIPELINE DE DINÁMICA MOLECULAR COMPLETADO")
    print(f"    Trayectoria final: {os.path.join(abs_md, 'md.xtc')}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ejecución Secuencial de DM (EM -> NVT -> NPT -> MD)")
    parser.add_argument("--complex-dir", default="results/complex")
    parser.add_argument("--md-dir", default="results/md_run")
    parser.add_argument("--mdp-dir", default="config/mdp")
    args = parser.parse_args()
    
    execute_md_pipeline(args.complex_dir, args.md_dir, args.mdp_dir)
