#!/usr/bin/env python3
"""Step 7: CLADES fit and Gillespie simulations, as two chained SLURM GPU jobs.

usage: python steps/07_clades.py config/hdbi.yaml            # write configs + job scripts only
       python steps/07_clades.py config/hdbi.yaml --submit   # and submit (Gillespie starts after the fit succeeds)
       python steps/07_clades.py config/hdbi.yaml --post     # queue steps 8 + 9 (CPU) after the fit and after Gillespie
                                                             # (job ids from clades/jobs.txt; with --submit, both at once)

reads runs/<run_name>/data/ (steps 1 and 6) and writes runs/<run_name>/clades/:
  fit.json, gillespie.json, logger_config.json   CLADES configs (gillespie.json gets the fitted model path at run time)
  fit.sbatch, gillespie.sbatch, post.sbatch, *.log   job scripts and logs; post.sbatch runs steps 8 and 9 and
                                                      appends a line to runs/<run_name>/STATUS.txt
  checkpoints/<name>/<time>/model_last.pt        fitted model (CLADES main.py)
  gillespie/checkpoints/.../models/clone_<i>/    simulations (CLADES main_gillespie.py)
"""
import json
import subprocess
import sys
from pathlib import Path

import yaml

cfg_path = Path(sys.argv[1]).resolve()
repo = cfg_path.parents[1]
cfg = yaml.safe_load(open(cfg_path))
cl = cfg['clades']
run = repo / cfg['runs_dir'] / cfg['run_name']
data, out = run / 'data', run / 'clades'
out.mkdir(parents=True, exist_ok=True)
params = json.load(open(data / 'clades_params.json'))
assert (data / 'graph_table.csv').exists(), 'run step 6 first'
name = f"{cfg['run_name']}_{cl['K_type']}"

json.dump({
    'version': 1, 'disable_existing_loggers': False,
    'formatters': {'simple': {'format': '%(message)s'},
                   'datetime': {'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s'}},
    'handlers': {'console': {'class': 'logging.StreamHandler', 'level': 'DEBUG', 'formatter': 'simple',
                             'stream': 'ext://sys.stdout'},
                 'info_file_handler': {'class': 'logging.handlers.RotatingFileHandler', 'level': 'INFO',
                                       'formatter': 'datetime', 'filename': 'info.log', 'maxBytes': 10485760,
                                       'backupCount': 20, 'encoding': 'utf8'}},
    'root': {'level': 'INFO', 'handlers': ['console', 'info_file_handler']},
}, open(out / 'logger_config.json', 'w'), indent=1)

json.dump({
    'name': name,
    'logger_config_path': str(out / 'logger_config.json'),
    'system': {'gpu_id': 0, 'seed': cl['seed']},
    'arch': {'type': 'ODEBlock', 'args': {'activation': cl['activation'], 'hidden_dim': cl['hidden_dim'],
                                          'K_type': cl['K_type']}},
    'base_trainer': {'epochs': cl['epochs'], 'save_dir': str(out), 'verbosity': 2},
    'user_trainer': {'adjoint': False, 'alphas': cl['alphas'],
                     'no_proliferation_pops': params['no_proliferation_pops'],
                     'no_apoptosis_pops': params['no_apoptosis_pops'],
                     't_observed': params['t_observed'], 'scaling_factor': params['scaling_factor'],
                     'ub_for_prol': cl['ub_for_prol'], 'ub_for_diff': cl['ub_for_diff']},
    'optimizer': {'learning_rate': cl['learning_rate'], 'lrs_ms': cl['lr_milestones']},
    'data_loader': {'type': 'ClonalTransDataLoader',
                    'args': {'data_dir': str(data), 'num_populations': params['num_populations'],
                             'annots': 'annotations.csv', 'graphs': 'graph_table.csv',
                             'day_zero': 'initial_condition.csv', 'cell_counts': 'kinetics.txt'}},
}, open(out / 'fit.json', 'w'), indent=1)

gl = cl['gillespie']
json.dump({
    'name': f'{name}_gillespie', 'logger_config_path': str(out / 'logger_config.json'), 'model_path': '',
    'time_interval': gl['time_interval'], 'concurrent': gl['concurrent'], 'epoch': gl['epoch'],
    't_cutoff': gl['t_cutoff'], 'base_trainer': {'save_dir': str(out / 'gillespie')},
}, open(out / 'gillespie.json', 'w'), indent=1)

head = f"""#!/bin/bash
#SBATCH -A {cl['slurm_gpu_account']}
#SBATCH -p {cl['slurm_gpu_partition']}
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
"""
(out / 'fit.sbatch').write_text(head + f"""#SBATCH -J clokin_fit
#SBATCH -t {cl['fit_hours']}:00:00
#SBATCH -o {out}/fit.%j.log
cd {cl['code_dir']}
{cl['python']} main.py --config {out}/fit.json
""")
(out / 'gillespie.sbatch').write_text(head + f"""#SBATCH -J clokin_gillespie
#SBATCH -t {gl['hours']}:00:00
#SBATCH -o {out}/gillespie.%j.log
ckpt=$(ls -t {out}/checkpoints/{name}/*/model_last.pt | head -1)
[ -z "$ckpt" ] && {{ echo "no fitted model"; exit 1; }}
echo "model: $ckpt"
{cl['python']} -c "import json, sys; f = '{out}/gillespie.json'; d = json.load(open(f)); d['model_path'] = sys.argv[1]; json.dump(d, open(f, 'w'), indent=1)" "$ckpt"
cd {cl['code_dir']}
{cl['python']} main_gillespie.py --config {out}/gillespie.json
""")
tools = cfg['tools']
(out / 'post.sbatch').write_text(f"""#!/bin/bash
#SBATCH -A {tools['slurm_cpu_account']}
#SBATCH -p {tools['slurm_cpu_partition']}
#SBATCH -N 1
#SBATCH -c 32
#SBATCH -J clokin_post
#SBATCH -t 06:00:00
#SBATCH -o {out}/post.%j.log
cd {repo}
{tools['python']} steps/08_diagnostics.py {cfg_path} && {tools['python']} steps/09_figures.py {cfg_path}
echo "$(date '+%Y-%m-%d %H:%M') post job $SLURM_JOB_ID ($1): steps 8-9 exit $?" >> {run}/STATUS.txt
""")
print(f'wrote configs and job scripts in {out}')

if '--submit' in sys.argv:
    fit = subprocess.run(['sbatch', '--parsable', str(out / 'fit.sbatch')], capture_output=True, text=True, check=True).stdout.strip()
    gil = subprocess.run(['sbatch', '--parsable', f'--dependency=afterok:{fit}', str(out / 'gillespie.sbatch')],
                         capture_output=True, text=True, check=True).stdout.strip()
    (out / 'jobs.txt').write_text(f'fit {fit}\ngillespie {gil} (afterok:{fit})\n')
    print(f'submitted fit {fit}, gillespie {gil} (after the fit)')

if '--post' in sys.argv:
    jobs = dict(line.split()[:2] for line in (out / 'jobs.txt').read_text().splitlines() if line.strip())
    lines = []
    for after in ['fit', 'gillespie']:
        jid = subprocess.run(['sbatch', '--parsable', f'--dependency=afterok:{jobs[after]}', str(out / 'post.sbatch'),
                              f'after {after}'], capture_output=True, text=True, check=True).stdout.strip()
        lines.append(f'post_{after} {jid} (afterok:{jobs[after]})')
    with open(out / 'jobs.txt', 'a') as f:
        f.write('\n'.join(lines) + '\n')
    print('submitted', '; '.join(lines))
