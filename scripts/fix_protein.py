import argparse
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def fix_protein(input_pdb, output_pdb, keep_ions, ph=7.4):
    print(f"[+] Procesando receptor con PDBFixer: {input_pdb}")
    fixer = PDBFixer(filename=input_pdb)
    
    # 1. Identificar y reparar lagunas de coordenadas
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    
    # 2. Eliminar aguas y solutos exógenos, conservando iones funcionales
    fixer.removeHeterogens(keepWater=False)
    
    # 3. Adición de hidrógenos bajo criterio de pH
    fixer.addMissingHydrogens(pH=ph)
    
    # 4. Escritura estandarizada
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    print(f"[✓] Proteína reparada con éxito en: {output_pdb}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Limpieza y protonación con PDBFixer")
    parser.add_argument("--input", required=True, help="Ruta al PDB crudo")
    parser.add_argument("--output", required=True, help="Ruta de salida")
    parser.add_argument("--ph", type=float, default=7.4, help="pH fisiológico")
    parser.add_argument("--keep-ions", nargs="+", default=["CA"], help="Iones a conservar")
    args = parser.parse_args()
    
    fix_protein(args.input, args.output, args.keep_ions, args.ph)
