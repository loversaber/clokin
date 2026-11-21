#!/usr/bin/env python3
"""Step 5c: candidate edges between populations, with clonal, PAGA and CoSpar evidence.

usage: python steps/05c_edges.py config/hdbi.yaml

candidates: edges of edges.previous_graph (optional) + PAGA links >= edges.paga_candidate, directed from the earlier to the later
population in the `populations` order. Evidence for each edge A -> B:
  clonal_share   share of A's first-timepoint clones with cells in B at a later timepoint (barcoded cells);
                 low when B is rare at later timepoints (see target_cells_after_<first timepoint>)
  paga           PAGA connectivity (step 5a)
  cospar_partner fate that B's fate couples to most (CoSpar, all timepoints; step 5b)
status:
  agree          clonal_share >= clonal_support and paga >= paga_support
  clonal only    clonal support, weak PAGA link
  paga only      PAGA support, clonal share below clonal_support
  no clonal data A has fewer than min_d3_clones first-timepoint clones; PAGA only
  no support     neither
  mature-mature  both populations terminal: not a differentiation step
writes runs/<run_name>/support/edges.csv (column `use` is filled in by step 6) and potency.csv (fate output of each
population's first-timepoint clones)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from common import load, population_of, selected

cfg_path, repo, cfg, run = load(sys.argv[1])
col, ec, cc = cfg['columns'], cfg['edges'], cfg['carta']
out = run / 'support'

pops = list(cfg['populations'])
rank = {p: i for i, p in enumerate(pops)}
pop_of = population_of(cfg)
fate_of_pop = {p: {cc['fates'][ct] for ct in cts} - {'Undifferentiated'} for p, cts in cfg['populations'].items()}
terminal = set(cfg['terminal'])
tps = list(cfg['timepoints'])

cells = pd.read_parquet(repo / cfg['inputs']['cells'])
keep = selected(cells, cfg, barcoded=True)
c = pd.DataFrame({'clone': cells.loc[keep, col['clone']].astype(str),
                  'pop': cells.loc[keep, col['celltype']].astype(str).map(pop_of),
                  'tp': cells.loc[keep, col['timepoint']].astype(str)})
first = c[c.tp == tps[0]].groupby('pop').clone.unique()
later = pd.crosstab(c.loc[c.tp != tps[0], 'clone'], c.loc[c.tp != tps[0], 'pop']).reindex(columns=pops, fill_value=0) > 0
n_first = first.map(len).reindex(pops).fillna(0).astype(int)
n_later = c[c.tp != tps[0]]['pop'].value_counts().reindex(pops).fillna(0).astype(int)


def share(a, b):
    cl = [x for x in first.get(a, []) if x in later.index]
    return later.loc[cl, b].mean() if n_first[a] >= ec['min_d3_clones'] and cl else np.nan


fates = sorted({f for s in fate_of_pop.values() for f in s})
fate_later = pd.DataFrame({f: later[[p for p in pops if f in fate_of_pop[p]]].any(axis=1) for f in fates})
pot = pd.DataFrame({p: fate_later.loc[[x for x in first.get(p, []) if x in fate_later.index]].mean()
                    for p in pops if n_first[p] > 0}).T
pot.insert(0, f'clones_{tps[0]}', n_first[pot.index])
pot.round(3).to_csv(out / 'potency.csv', index_label='population')

paga = pd.read_csv(out / 'paga_connectivity.csv', index_col=0)
cos = pd.read_csv(out / 'cospar_coupling_all.csv', index_col=0)
prev = (pd.read_csv(repo / ec['previous_graph'], index_col=0) if ec.get('previous_graph')
        else pd.DataFrame(0, index=pops, columns=pops))

cand = {(a, b) for a in prev.index for b in prev.columns if prev.loc[a, b] == 1}
for i, a in enumerate(pops):
    for b in pops[i + 1:]:
        if paga.loc[a, b] >= ec['paga_candidate']:
            cand.add((a, b))
cand = {(a, b) if rank[a] < rank[b] else (b, a) for a, b in cand}

rows = []
for a, b in sorted(cand, key=lambda e: (rank[e[0]], rank[e[1]])):
    s, g = share(a, b), paga.loc[a, b]
    fb = sorted(fate_of_pop[b])
    partner = ''
    if len(fb) == 1 and fb[0] in cos.index:
        cp = cos.loc[fb[0]].drop(fb[0])
        partner = f'{cp.idxmax()} {cp.max():.3f}'
    if a in terminal and b in terminal:
        status = 'mature-mature'
    elif np.isnan(s):
        status = 'no clonal data' if g >= ec['paga_support'] else 'no support'
    elif s >= ec['clonal_support'] and g >= ec['paga_support']:
        status = 'agree'
    elif s >= ec['clonal_support']:
        status = 'clonal only'
    elif g >= ec['paga_support']:
        status = 'paga only'
    else:
        status = 'no support'
    rows.append({'source': a, 'target': b, 'in_previous_graph': int(prev.loc[a, b] == 1 if a in prev.index and b in prev.columns else 0),
                 f'source_clones_{tps[0]}': n_first[a], 'clonal_share': round(s, 3) if not np.isnan(s) else np.nan,
                 f'target_cells_after_{tps[0]}': n_later[b], 'paga': round(g, 3), 'target_fate': '+'.join(fb), 'cospar_partner': partner, 'status': status,
                 'use': ''})
edges = pd.DataFrame(rows)
edges.to_csv(out / 'edges.csv', index=False)

pd.set_option('display.width', 250)
pd.set_option('display.max_rows', 200)
print(pot.round(2).to_string(), '\n')
print(edges.drop(columns='use').to_string(index=False))
print('\n' + edges.status.value_counts().to_string())
