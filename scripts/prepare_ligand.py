import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

def prepare_ligand(input_pdb, resname, output_pdbqt, output_sdf):
    print(f"[+] Extrayendo coordenadas nativas de '{resname}' desde {input_pdb}...")

    # 1. Filtrar lineas del ligando desde el PDB para no sobrecargar RDKit
    with open(input_pdb, "r") as f:
        ligand_lines = [
            line for line in f 
            if line.startswith(("HETATM", "ATOM")) and line[17:20].strip() == resname
        ]

    if not ligand_lines:
        raise ValueError(f"No se encontro el residuo '{resname}' en {input_pdb}")

    # 2. Cargar bloque PDB sin protones (manteniendo las coordenadas 3D del cristal)
    pdb_block = "".join(ligand_lines)
    raw_mol = Chem.MolFromPDBBlock(pdb_block, removeHs=True, sanitize=False)
    if raw_mol is None:
        raise ValueError("Error critico: RDKit no pudo construir la molecula desde el bloque PDB.")

    # 3. Asignar ordenes de enlace y carga cationica (+1) usando plantilla SMILES
    # Benzamidina protonada en el grupo amidinio (pH fisiologico)
    template_smiles = "c1ccccc1C(=[NH2+])N"
    template = Chem.MolFromSmiles(template_smiles)
    
    # Asignacion de enlaces respetando la geometria tridimensional
    assigned_mol = AllChem.AssignBondOrdersFromTemplate(template, raw_mol)
    
    # 4. Anadir hidrogenos con coordenadas 3D coherentes
    mol_with_h = Chem.AddHs(assigned_mol, addCoords=True)

    # 5. Exportar SDF cristalográfico de referencia (necesario para autobox y RMSD)
    with Chem.SDWriter(output_sdf) as writer:
        writer.write(mol_with_h)
    print(f"[✓] SDF de referencia cristalográfica generado: {output_sdf}")

    # 6. Generar PDBQT mediante Meeko (API v0.5+)
    preparator = MoleculePreparation()
    mol_setup_list = preparator.prepare(mol_with_h)
    
    # Extraer el primer setup y escribir el string
    writer = PDBQTWriterLegacy()
    pdbqt_str = writer.write_string(mol_setup_list[0])

    with open(output_pdbqt, "w") as f:
        f.write(pdbqt_str[0])  # write_string devuelve una tupla, tomamos el string
    print(f"[✓] PDBQT generado con exito: {output_pdbqt}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extraccion y preparacion rigurosa de ligando")
    parser.add_argument("--input-pdb", required=True, help="Ruta al archivo PDB crudo")
    parser.add_argument("--resname", default="BEN", help="Codigo de residuo (ej. BEN)")
    parser.add_argument("--output-pdbqt", required=True, help="Ruta de salida del PDBQT")
    parser.add_argument("--output-sdf", required=True, help="Ruta de salida del SDF de referencia")
    args = parser.parse_args()

    prepare_ligand(args.input_pdb, args.resname, args.output_pdbqt, args.output_sdf)
