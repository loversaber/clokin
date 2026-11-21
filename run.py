#!/usr/bin/env python3
"""Run clokin in two stages, with a decision between them.

usage: python run.py tree config/<dataset>.yaml       # steps 1-5 as one SLURM CPU job (--here: run in this shell)
       python run.py kinetics config/<dataset>.yaml   # steps 6-7: graph, then CLADES fit -> Gillespie on the GPU,
                                                      # steps 8-9 after the fit and again after Gillespie
       python run.py status config/<dataset>.yaml     # jobs, STATUS.txt and which outputs exist

between the stages: look at runs/<run_name>/carta/ and support/ (edges.csv, trees.pdf, paga_umap.pdf), then write the
chosen edges in the `graph` section of the config (and carta.k if the tree should use another k).
`python` must be the environment in tools.python (it reads the YAML config).
"""
import subprocess
import sys
from pathlib import Path

repo = Path(__file__).resolve().parent
sys.path.insert(0, str(repo / 'steps'))
from common import load

if len(sys.argv) < 3 or sys.argv[1] not in ('tree', 'kinetics', 'status'):
    sys.exit(__doc__)
stage = sys.argv[1]
cfg_path, _, cfg, run = load(sys.argv[2])
run.mkdir(parents=True, exist_ok=True)
tools, py = cfg['tools'], cfg['tools']['python']
tree_steps = ['01_inputs', '02_carta_matrix', '03_carta_sweep', '04_carta_tree', '05a_paga', '05b_cospar',
              '05c_edges', '05d_trees']


def sh(*args):
    return subprocess.run(list(args), capture_output=True, text=True, check=True).stdout.strip()


if stage == 'tree':
    lines = [f'{py} {repo}/steps/{s}.py {cfg_path}' for s in tree_steps]
    if '--here' in sys.argv:
        for line in lines:
            print('>>', line.split('/')[-1], flush=True)
            subprocess.run(line.split(), check=True)
        sys.exit()
    job = run / 'tree.sbatch'
    job.write_text(f"""#!/bin/bash
#SBATCH -A {tools['slurm_cpu_account']}
#SBATCH -p {tools['slurm_cpu_partition']}
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -J clokin_tree
#SBATCH -t 12:00:00
#SBATCH -o {run}/tree.%j.log
set -e
cd {repo}
""" + '\n'.join(f'echo ">> {s}"; {line}' for s, line in zip(tree_steps, lines)) + f"""
echo "$(date '+%Y-%m-%d %H:%M') tree job $SLURM_JOB_ID: steps 1-5 done" >> {run}/STATUS.txt
""")
    jid = sh('sbatch', '--parsable', str(job))
    with open(run / 'jobs.txt', 'a') as f:
        f.write(f'tree {jid}\n')
    print(f'submitted tree job {jid}; log {run}/tree.{jid}.log')

elif stage == 'kinetics':
    assert (run / 'support' / 'edges.csv').exists(), 'run the tree stage first'
    assert cfg.get('graph'), 'write the chosen edges in the `graph` section of the config first'
    subprocess.run([py, f'{repo}/steps/06_graph.py', str(cfg_path)], check=True)
    subprocess.run([py, f'{repo}/steps/07_clades.py', str(cfg_path), '--submit', '--post'], check=True)

else:
    jobs = {}
    for f in [run / 'jobs.txt', run / 'clades' / 'jobs.txt']:
        if f.exists():
            jobs.update(dict(line.split()[:2] for line in f.read_text().splitlines() if line.strip()))
    if jobs:
        print(sh('sacct', '-X', '-j', ','.join(jobs.values()), '-o', 'JobID,JobName%18,State,Elapsed,End'))
    if (run / 'STATUS.txt').exists():
        print('\n' + (run / 'STATUS.txt').read_text())
    expect = {'carta': ['tree_k*.pdf', 'pareto.pdf'], 'support': ['edges.csv', 'trees.pdf', 'paga_umap.pdf',
              'cospar_coupling.pdf'], 'data': ['graph_table.csv'], 'diagnostics': ['fit_by_row.csv', 'fate_prob.csv'],
              'figures': ['fit_overview.pdf', 'fate_prob.pdf']}
    for d, pats in expect.items():
        print(f'{d + "/":13s}' + '  '.join(f'{p} {"ok" if list((run / d).glob(p)) else "-"}' for p in pats))
