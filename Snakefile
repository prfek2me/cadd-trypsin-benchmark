import yaml

configfile: "config/config.yaml"

PDB_ID = config["system"]["pdb_id"]
LIG_RESNAME = config["system"]["raw_ligand_resname"]

rule all:
    input:
        "results/analysis/md_stability_panel.png"

rule fix_protein:
    input:
        f"data/raw/{PDB_ID}.pdb"
    output:
        "results/protein_fixed.pdb"
    shell:
        "python scripts/fix_protein.py --input {input} --output {output} --keep-ions CA --ph 7.4"

rule prepare_ligand:
    input:
        f"data/raw/{PDB_ID}.pdb"
    output:
        pdbqt="results/ligand.pdbqt",
        sdf="results/ligand_ref.sdf"
    shell:
        "python scripts/prepare_ligand.py --input-pdb {input} --resname {LIG_RESNAME} --output-pdbqt {output.pdbqt} --output-sdf {output.sdf}"

rule run_docking:
    input:
        receptor="results/protein_fixed.pdb",
        lig_pdbqt="results/ligand.pdbqt",
        lig_sdf="results/ligand_ref.sdf"
    output:
        out_sdf="results/docking_poses.sdf",
        best_pose="results/best_pose.sdf",
        log="results/docking.log"
    shell:
        "python scripts/run_docking.py --receptor {input.receptor} --ligand {input.lig_sdf} --autobox {input.lig_sdf} --out-sdf {output.out_sdf} --out-log {output.log} --best-pose {output.best_pose}"

rule build_topology:
    input:
        ref_sdf="results/ligand_ref.sdf",
        pose_sdf="results/best_pose.sdf"
    output:
        directory("results/topology")
    shell:
        "python scripts/build_topology.py --ref-sdf {input.ref_sdf} --pose-sdf {input.pose_sdf} --basename {LIG_RESNAME} --net-charge 1 --out-dir {output}"

rule build_complex:
    input:
        protein="results/protein_fixed.pdb",
        ligand_dir="results/topology"
    output:
        directory("results/complex")
    shell:
        "python scripts/build_complex.py --protein {input.protein} --ligand-dir {input.ligand_dir} --out-dir {output}"

rule write_mdps:
    output:
        directory("config/mdp")
    shell:
        "python scripts/write_mdps.py --out-dir {output} --time-ns 20 --seed 123"

rule run_md:
    input:
        complex_dir="results/complex",
        mdp_dir="config/mdp"
    output:
        xtc="results/md_run/md.xtc",
        tpr="results/md_run/md.tpr"
    shell:
        "python scripts/run_md.py --complex-dir {input.complex_dir} --md-dir results/md_run --mdp-dir {input.mdp_dir}"

rule stabilize_trajectory:
    input:
        xtc="results/md_run/md.xtc",
        tpr="results/md_run/md.tpr"
    output:
        ref_pdb="results/md_run/md_ref.pdb",
        xtc_nopbc="results/md_run/md_nopbc.xtc"
    params:
        md_dir="results/md_run"
    shell:
        """
        echo "0" | gmx trjconv -s {input.tpr} -f {input.xtc} -o {output.ref_pdb} -dump 0
        echo "0" | gmx trjconv -s {input.tpr} -f {input.xtc} -o {params.md_dir}/md_whole.xtc -pbc whole
        echo -e "1\\n0" | gmx trjconv -s {input.tpr} -f {params.md_dir}/md_whole.xtc -o {output.xtc_nopbc} -pbc mol -center
        rm -f {params.md_dir}/md_whole.xtc
        """

rule analyze_trajectory:
    input:
        tpr="results/md_run/md_ref.pdb",
        xtc="results/md_run/md_nopbc.xtc"
    output:
        "results/analysis/md_stability_panel.png"
    shell:
        "python scripts/analyze_trajectory.py --tpr {input.tpr} --xtc {input.xtc} --out-dir results/analysis"
