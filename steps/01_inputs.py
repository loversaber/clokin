#!/usr/bin/env python3
"""Step 1: CLADES input files from the cell table, the meta-clone labels and the scaling factors.

usage: python steps/01_inputs.py config/<dataset>.yaml

writes runs/<run_name>/data/:
  kinetics.txt            (groups + BG) x (timepoints x populations), scaled cell counts; with clades.version '1.2'
                          the day-0 block comes first (CLADES v1.2.0 reads day 0 from here)
  initial_condition.csv   day-0 cells per row (CLADES before v1.2.0 puts them in the root population)
  day0.csv                day-0 cells per row and population
  annotations.csv         row names (groups, BG) and population names
  counts_long.csv         the same counts as a table, plus raw cell numbers
  clades_params.json      t_observed, scaling_factor, no_proliferation_pops, no_apoptosis_pops
and a copy of the config in runs/<run_name>/
"""
import json
import sys

import numpy as np
import pandas as pd

from common import load, population_of, read_table, save_config, selected

cfg_path, repo, cfg, run = load(sys.argv[1])
col = cfg['columns']
out = run / 'data'
out.mkdir(parents=True, exist_ok=True)
save_config(cfg, run / 'config.yaml')

pops = list(cfg['populations'])
pop_of = population_of(cfg)
days = cfg['timepoints']

cells = pd.read_parquet(repo / cfg['inputs']['cells'])
labels = pd.read_csv(repo / cfg['inputs']['labels'])
scf = read_table(repo / cfg['inputs']['scaling'])

obs = pd.DataFrame({
    'clone': cells[col['clone']].astype(str).replace({'nan': None, 'None': None, '': None}),
    'pop': cells[col['celltype']].astype(str).map(pop_of),
    'tp': cells[col['timepoint']].astype(str),
    'exp': cells[col['experiment']].astype(str),
})
keep = selected(cells, cfg) & obs['pop'].notna() & obs.exp.isin(scf.index)
unmapped = sorted(set(cells.loc[obs['pop'].isna(), col['celltype']].astype(str)))
obs = obs[keep].copy()
obs['sf'] = scf.stack().reindex(list(zip(obs.exp, obs.tp))).values
assert obs.sf.notna().all()

group_of = labels.set_index(col['label_clone'])[col['label_group']].astype(str)
obs['group'] = obs.clone.map(group_of)
groups = sorted(group_of.unique(), key=lambda g: (len(g), g))
if cfg.get('keep_groups'):
    groups = [g for g in groups if g in set(map(str, cfg['keep_groups']))]
names = [f"{cfg['group_prefix']}{g}" for g in groups] + ['BG']

print(f'cells kept: {len(obs):,} of {len(cells):,} ({obs.exp.nunique()} experiments); cell types not used: {unmapped}')
print(f'labelled clones: {len(group_of):,}, with kept cells: {obs.clone.isin(group_of.index).groupby(obs.clone).any().sum():,}')

labelled = obs[obs.group.notna()]
mc = labelled[labelled.group.isin(groups)]
if cfg.get('clades_rounding'):
    # as the CLADES cord blood notebook: each clone's scaled count rounded (float32), background truncated
    per_clone = mc.groupby(['group', 'clone', 'tp', 'pop']).agg(n_raw=('sf', 'size'), sf=('sf', 'first')).reset_index()
    per_clone['n_scaled'] = np.round(per_clone.n_raw.astype(np.float32) * per_clone.sf.astype(np.float32))
    mc_long = per_clone.groupby(['group', 'tp', 'pop'])[['n_scaled', 'n_raw']].sum().reset_index()
    bg_long = obs.groupby(['tp', 'pop']).agg(n_raw=('sf', 'size'), sf=('sf', 'first')).reset_index()
    bg_long['n_scaled'] = np.floor(bg_long.n_raw * bg_long.sf)
    bg_long = bg_long.drop(columns='sf')
else:
    mc_long = mc.groupby(['group', 'tp', 'pop']).agg(n_scaled=('sf', 'sum'), n_raw=('sf', 'size')).reset_index()
    bg_long = obs.groupby(['tp', 'pop']).agg(n_scaled=('sf', 'sum'), n_raw=('sf', 'size')).reset_index()
long = pd.concat([mc_long, bg_long.assign(group='BG')])
long['row'] = long.group.map(dict(zip(groups + ['BG'], names)))
long.to_csv(out / 'counts_long.csv', index=False)

tps = list(days)
kin = (long.pivot_table(index='row', columns=['tp', 'pop'], values='n_scaled', aggfunc='sum')
       .reindex(index=names, columns=pd.MultiIndex.from_product([tps, pops])).fillna(0).round())
np.savetxt(out / 'kinetics.txt', kin.values, fmt='%d')

# day 0: founders = one cell per clone (background: all clones); all in the root population, or with
# founders_split (CLADES v1.2.0) split over day0_populations as the row's first-timepoint cells in them
n_founders = group_of.value_counts().reindex(groups).tolist() + [cfg.get('day0_background_founders', obs.clone.dropna().nunique())]
day0 = pd.DataFrame(0.0, index=names, columns=pops)
if cfg['day0'] == 'founders':
    day0[pops[0]] = n_founders
elif cfg['day0'] == 'founders_split':
    d0_pops = cfg['day0_populations']
    base = obs[obs.clone.notna()]
    if cfg.get('day0_split_clones', 'all') == 'seen_at_all_timepoints':
        # CLADES v1.2.0: only clones with cells at every timepoint
        seen = base.groupby('clone').tp.nunique()
        base = base[base.clone.isin(seen[seen == len(tps)].index)]
    first = base[base.tp == tps[0]]
    for name, cells_g, n in zip(names, [first[first.group == g] for g in groups] + [first], n_founders):
        share = cells_g['pop'].value_counts().reindex(d0_pops).fillna(0)
        if share.sum():
            day0.loc[name, d0_pops] = np.round(share / share.sum() * n).values
        else:
            day0.loc[name, pops[0]] = n
elif cfg['day0'] == 'd3_total':
    day0[pops[0]] = kin[tps[0]].sum(axis=1).values
else:
    raise ValueError(f"day0 must be 'founders', 'founders_split' or 'd3_total', not {cfg['day0']}")
day0 = day0.round().astype(int)
d0 = day0.sum(axis=1).tolist()
day0.to_csv(out / 'day0.csv')
pd.DataFrame({'leiden': d0}).to_csv(out / 'initial_condition.csv')
if str(cfg.get('clades', {}).get('version', 'main')) == '1.2':
    np.savetxt(out / 'kinetics.txt', np.hstack([day0.values, kin.values]), fmt='%d')

n = max(len(names), len(pops))
pd.DataFrame({'clones': names + [''] * (n - len(names)),
              'populations': pops + [''] * (n - len(pops))}).to_csv(out / 'annotations.csv', index=False)

raw = obs.groupby('tp').size()
scaled = obs.groupby('tp').sf.sum()
params = {
    'num_populations': len(pops),
    't_observed': [0] + [days[t] for t in tps],
    'scaling_factor': [1.0] + [round(float(scaled[t] / raw[t]), 4) for t in tps],
    'no_proliferation_pops': [int(p in cfg.get('no_proliferation', cfg['terminal'])) for p in pops],
    'no_apoptosis_pops': [int(p in cfg['early']) for p in pops],
}
json.dump(params, open(out / 'clades_params.json', 'w'), indent=1)

summary = pd.DataFrame({t: kin[t].sum(axis=1) for t in tps}).round().astype(int)
summary.insert(0, 'day0', np.round(d0).astype(int))
summary[f'pops_{tps[-1]}'] = (kin[tps[-1]] > 0).sum(axis=1)
pd.set_option('display.width', 200)
print(f'\nkinetics {kin.shape} -> {out}')
print(json.dumps({k: params[k] for k in ['t_observed', 'scaling_factor']}))
print(summary.to_string())
