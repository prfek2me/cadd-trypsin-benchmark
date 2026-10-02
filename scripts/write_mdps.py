import argparse
import os

def write_em(out_dir):
    content = """title       = Minimizacion (AMBER)
integrator  = steep
emtol       = 1000.0
nsteps      = 50000
nstlist     = 10
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
vdwtype     = Cut-off
rvdw        = 1.0
DispCorr    = EnerPres
pbc         = xyz
"""
    with open(os.path.join(out_dir, "em.mdp"), "w") as f: f.write(content)

def write_nvt(out_dir, seed=123):
    content = f"""title       = Equilibracion NVT (AMBER)
define      = -DPOSRES
integrator  = md
dt          = 0.002
nsteps      = 50000      ; 100 ps
nstxout     = 0
nstvout     = 0
nstenergy   = 5000
nstlog      = 5000
continuation = no
constraint_algorithm = lincs
constraints = h-bonds
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
vdwtype     = Cut-off
rvdw        = 1.0
DispCorr    = EnerPres
pbc         = xyz
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1     0.1
ref_t       = 300     300
pcoupl      = no
gen_vel     = yes
gen_temp    = 300
gen_seed    = {seed}
"""
    with open(os.path.join(out_dir, "nvt.mdp"), "w") as f: f.write(content)

def write_npt(out_dir):
    content = """title       = Equilibracion NPT (AMBER)
define      = -DPOSRES
integrator  = md
dt          = 0.002
nsteps      = 50000      ; 100 ps
nstxout     = 0
nstvout     = 0
nstenergy   = 5000
nstlog      = 5000
continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
vdwtype     = Cut-off
rvdw        = 1.0
DispCorr    = EnerPres
pbc         = xyz
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1     0.1
ref_t       = 300     300
pcoupl      = c-rescale
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
refcoord_scaling = com
"""
    with open(os.path.join(out_dir, "npt.mdp"), "w") as f: f.write(content)

def write_md(out_dir, time_ns=20):
    nsteps = int((time_ns * 1000) / 0.002)
    content = f"""title       = Produccion Dinamica Molecular {time_ns} ns
integrator  = md
dt          = 0.002
nsteps      = {nsteps}
nstxout     = 0
nstvout     = 0
nstenergy   = 5000
nstlog      = 5000
nstxout-compressed = 5000
compressed-x-grps  = System
continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
vdwtype     = Cut-off
rvdw        = 1.0
DispCorr    = EnerPres
pbc         = xyz
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1     0.1
ref_t       = 300     300
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
"""
    with open(os.path.join(out_dir, "md.mdp"), "w") as f: f.write(content)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default="config/mdp")
    parser.add_argument("--time-ns", type=float, default=20.0)
    parser.add_argument("--seed", type=int, default=123)
    args = parser.parse_args()
    
    os.makedirs(args.out_dir, exist_ok=True)
    write_em(args.out_dir)
    write_nvt(args.out_dir, args.seed)
    write_npt(args.out_dir)
    write_md(args.out_dir, args.time_ns)
    print(f"[✓] Archivos MDP generados en: {args.out_dir}/")
