#!/usr/bin/env python3
"""Step 3: CARTA ILP (tree mode) for k = 1..kmax progenitors.

usage: python steps/03_carta_sweep.py config/<dataset>.yaml            # run here
       python steps/03_carta_sweep.py config/<dataset>.yaml --submit   # as a SLURM CPU job

writes runs/<run_name>/carta/sweep/k<k>_{prog,assignment,labeling}.csv (larry.py output) and
k<k>.log; a k that already has its _prog.csv is skipped, so the step can be re-run after a timeout.
"""
import os
import subprocess
import sys
import time
from pathlib import Path


from common import load

cfg_path, repo, cfg, run = load(sys.argv[1])
tools, cc = cfg['tools'], cfg['carta']
carta = run / 'carta'
sweep = carta / 'sweep'
sweep.mkdir(parents=True, exist_ok=True)

if '--submit' in sys.argv:
    job = sweep / 'sweep.sbatch'
    job.write_text(f"""#!/bin/bash
#SBATCH -J clokin_carta
#SBATCH -A {tools['slurm_cpu_account']}
#SBATCH -p {tools['slurm_cpu_partition']}
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 12:00:00
#SBATCH -o {sweep}/sweep.%j.log
cd {repo}
{tools['python']} {Path(__file__).resolve()} {cfg_path}
""")
    print(subprocess.run(['sbatch', str(job)], capture_output=True, text=True, check=True).stdout.strip())
    sys.exit()

env = dict(os.environ, GRB_LICENSE_FILE=tools['gurobi_license'])
larry = Path(tools['carta_dir']) / 'scripts' / 'larry.py'
for k in range(1, cc['kmax'] + 1):
    prefix = sweep / f'k{k}'
    if Path(f'{prefix}_prog.csv').exists():
        print(f'k={k}: done before')
        continue
    t0 = time.time()
    cmd = [tools['carta_python'], str(larry), '-i', str(carta / 'clone_type_mat.csv'),
           '-k', str(k), '-t', str(cc['threshold']), '-r', '-o', str(prefix)]
    with open(f'{prefix}.log', 'w') as log:
        r = subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT)
    print(f'k={k}: {"ok" if r.returncode == 0 else "FAILED, see " + str(prefix) + ".log"} ({time.time() - t0:.0f} s)', flush=True)
    if r.returncode != 0:
        sys.exit(1)
