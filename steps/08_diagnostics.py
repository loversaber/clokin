#!/usr/bin/env python3
"""Step 8: diagnostics of the CLADES fit (and of the Gillespie simulations when they are complete).

usage: python steps/08_diagnostics.py config/hdbi.yaml   (re-runs itself with clades.python, on CPU)

writes runs/<run_name>/diagnostics/:
  fit_by_row.csv        per row (meta-clone, BG): observed / predicted totals per timepoint, log2(pred/obs),
                        relative error, Spearman, cells and populations at the last timepoint, trust label
  trajectories.csv      predicted cells per row and population every half day
  composition_last.csv  observed and predicted share of each population at the last timepoint
  rates.csv             fitted rates per row and timepoint: differentiation (source -> target) and net growth
  checks.csv            net growth of no-proliferation populations > 0, net loss of early populations, rates near the bound
  convergence.csv       loss terms over the epochs (from the fit log)
  fate_prob.csv         Gillespie: P(a lineage from one root cell reaches each population) per row, with trust label
                        (simulations counted once per seed, as in CLADES get_fate_prob; res_div.txt can hold duplicates)
  divisions.csv         Gillespie: mean number of divisions from the root to reach each population
trust: a row is 'use' when it has at least diagnostics.min_cells_last cells in at least diagnostics.min_pops_last
populations at the last timepoint; otherwise 'do not quote' (its rates and fate probabilities follow the penalty
terms and the graph more than the data).
"""
import json
import os
import sys
from pathlib import Path

cfg_path = Path(sys.argv[1]).resolve()
repo = cfg_path.parents[1]
if '--json' not in sys.argv:
    from common import load
    _, _, cfg, run = load(cfg_path)
    js = run / 'diagnostics' / 'config.json'
    js.parent.mkdir(parents=True, exist_ok=True)
    json.dump(cfg, open(js, 'w'))
    py = cfg['clades']['python']
    os.execv(py, [py, '-u', os.path.abspath(__file__), str(cfg_path), '--json', str(js)])

import ast
import glob
import re
import warnings

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr

warnings.filterwarnings('ignore')
cfg = json.load(open(sys.argv[sys.argv.index('--json') + 1]))
cl, dg = cfg['clades'], cfg['diagnostics']
run = repo / cfg['runs_dir'] / cfg['run_name']
out = run / 'diagnostics'
sys.path.insert(0, cl['code_dir'])
os.chdir(cl['code_dir'])

pops = list(cfg['populations'])
rows = [r for r in pd.read_csv(run / 'data' / 'annotations.csv').clones.dropna()]
days = [0] + list(cfg['timepoints'].values())
tps = ['D0'] + list(cfg['timepoints'])
name = f"{cfg['run_name']}_{cl['K_type']}"
ckpt = sorted(glob.glob(str(run / 'clades' / 'checkpoints' / name / '*' / 'model_last.pt')))[-1]
print('model:', ckpt, flush=True)

m = torch.load(ckpt, map_location='cpu', weights_only=False)
m.config['system']['gpu_id'] = 'cpu'
for k in dir(m):
    a = getattr(m, k, None)
    if hasattr(a, 'gpu_id'):
        try:
            a.gpu_id = 'cpu'
        except Exception:
            pass
m.gpu_id = 'cpu'
with torch.no_grad():
    obs = m.N.numpy()
    pred = m.eval_model(m.t_observed).numpy()

grid = np.round(np.arange(0, days[-1] + 0.01, 0.5), 2)
with torch.no_grad():
    curve = m.eval_model(torch.tensor(grid, dtype=torch.float32)).numpy()
pd.DataFrame([{'row': r, 'day': t, 'population': pops[j], 'predicted': round(float(curve[ti, i, j]), 2)}
              for ti, t in enumerate(grid) for i, r in enumerate(rows) for j in range(len(pops))]
             ).to_csv(out / 'trajectories.csv', index=False)

last = len(days) - 1
fit = []
for i, r in enumerate(rows):
    o, p = obs[1:, i].ravel(), pred[1:, i].ravel()
    nz = o > 0
    d = {'row': r}
    for t in range(1, len(days)):
        d[f'obs_{tps[t]}'] = round(float(obs[t, i].sum()))
        d[f'pred_{tps[t]}'] = round(float(pred[t, i].sum()))
        d[f'log2_pred_obs_{tps[t]}'] = round(float(np.log2((pred[t, i].sum() + 1) / (obs[t, i].sum() + 1))), 2)
    d['rel_error'] = round(float(np.abs(p[nz] - o[nz]).sum() / o[nz].sum()), 3) if nz.any() else np.nan
    d['spearman'] = round(float(spearmanr(o, p)[0]), 3)
    d['cells_last'] = round(float(obs[last, i].sum()))
    d['pops_last'] = int((obs[last, i] > 0).sum())
    d['trust'] = 'use' if d['cells_last'] >= dg['min_cells_last'] and d['pops_last'] >= dg['min_pops_last'] else 'do not quote'
    fit.append(d)
fit = pd.DataFrame(fit)
fit.to_csv(out / 'fit_by_row.csv', index=False)
trust = dict(zip(fit.row, fit.trust))

comp = []
for i, r in enumerate(rows):
    o = obs[last, i] / max(obs[last, i].sum(), 1)
    p = pred[last, i] / max(pred[last, i].sum(), 1e-9)
    comp += [{'row': r, 'population': pops[j], 'observed': round(float(o[j]), 4), 'predicted': round(float(p[j]), 4)}
             for j in range(len(pops))]
pd.DataFrame(comp).to_csv(out / 'composition_last.csv', index=False)

graph = pd.read_csv(run / 'data' / 'graph_table.csv', index_col=0).loc[pops, pops].values
rates, checks = [], []
term, early = set(cfg.get('no_proliferation', cfg['terminal'])), set(cfg['early'])
with torch.no_grad():
    for t in days[1:]:
        k = m.get_matrix_K(K_type=cl['K_type'], eval=True, tpoint=float(t)).numpy()
        for i, r in enumerate(rows):
            for a in range(len(pops)):
                rates.append({'row': r, 'day': t, 'source': pops[a], 'target': pops[a], 'kind': 'net growth',
                              'rate': round(float(k[i, a, a]), 4)})
                for b in np.where(graph[a] == 1)[0]:
                    rates.append({'row': r, 'day': t, 'source': pops[a], 'target': pops[b], 'kind': 'differentiation',
                                  'rate': round(float(k[i, a, b]), 4)})
rates = pd.DataFrame(rates)
rates.to_csv(out / 'rates.csv', index=False)
g = rates[rates.kind == 'net growth']
for _, x in g[g.source.isin(term) & (g.rate > 0.1)].iterrows():
    checks.append({**x[['row', 'day', 'source', 'rate']], 'check': 'no-proliferation population grows'})
for _, x in g[g.source.isin(early) & (g.rate < -0.1)].iterrows():
    checks.append({**x[['row', 'day', 'source', 'rate']], 'check': 'early population shrinks'})
dif = rates[rates.kind == 'differentiation']
for _, x in dif[dif.rate >= 0.9 * cl['ub_for_diff']].iterrows():
    checks.append({**x[['row', 'day', 'source', 'rate']], 'target': x.target, 'check': 'differentiation rate near bound'})
for _, x in g[g.rate >= 0.9 * cl['ub_for_prol']].iterrows():
    checks.append({**x[['row', 'day', 'source', 'rate']], 'check': 'net growth near bound'})
checks = pd.DataFrame(checks)
if len(checks):
    checks['trust'] = checks.row.map(trust)
checks.to_csv(out / 'checks.csv', index=False)

logs = sorted(glob.glob(str(run / 'clades' / 'fit.*.log')))
conv = []
if logs:
    txt = open(logs[-1], errors='ignore').read().replace('\r', '\n')
    for line in txt.split('\n'):
        ep = re.search(r'\|\s*(\d+)/(\d+) \[', line)
        if ep and 'L/' in line:
            conv.append({'epoch': int(ep.group(1)), **{a: float(b) for a, b in re.findall(r'L/(\w+) ([-\d.]+)', line)}})
conv = pd.DataFrame(conv).drop_duplicates('epoch', keep='last')
conv.to_csv(out / 'convergence.csv', index=False)

pd.set_option('display.width', 250)
print(fit.to_string(index=False))
if len(conv):
    obs_cols = [c for c in conv.columns if re.fullmatch(r'L\d+', c)]
    end, ref = conv.iloc[-1], conv[conv.epoch <= 0.8 * conv.epoch.max()].iloc[-1]
    print('\nloss at the end vs at 80 % of the epochs:',
          ', '.join(f'{c} {end[c]:.1f} ({100 * (end[c] - ref[c]) / ref[c]:+.1f} %)' for c in obs_cols))
print(f'\nchecks: {len(checks)} flags')
if len(checks):
    print(checks.groupby(['check', 'trust']).size().to_string())

gdir = glob.glob(str(run / 'clades' / 'gillespie' / 'checkpoints' / '*' / '*' / 'models'))
dirs = sorted(glob.glob(f'{gdir[0]}/clone_*'), key=lambda x: int(x.split('_')[-1])) if gdir else []
if len(dirs) != len(rows):
    print(f'\nGillespie not complete ({len(dirs)}/{len(rows)} rows): fate probabilities skipped')
    sys.exit()

import pl
fp, dv = [], []
for i, d in enumerate(dirs):
    if not os.path.exists(f'{d}/res_div.txt'):
        pl.get_num_div(pops, gillespie_dir=d)
    sims = {}
    for line in open(f'{d}/res_div.txt'):
        if line.strip():
            x = ast.literal_eval(ast.literal_eval(line))
            sims[x.pop('seed')] = x
    sims = list(sims.values())
    n = len(sims)
    reach = {p: [min(s[p]) for s in sims if p in s] for p in pops}
    fp.append({'row': rows[i], 'simulations': n, **{p: round(len(reach[p]) / n, 4) for p in pops}, 'trust': trust[rows[i]]})
    dv.append({'row': rows[i], **{p: round(float(np.mean(reach[p])), 2) if reach[p] else np.nan for p in pops},
               'trust': trust[rows[i]]})
pd.DataFrame(fp).to_csv(out / 'fate_prob.csv', index=False)
pd.DataFrame(dv).to_csv(out / 'divisions.csv', index=False)
print('\nP(reach) from one root cell, terminal populations:')
print(pd.DataFrame(fp).set_index('row')[cfg['terminal'] + ['trust']].to_string())
