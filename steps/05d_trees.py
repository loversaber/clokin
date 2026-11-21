#!/usr/bin/env python3
"""Step 5d: fate similarity and fate trees from CARTA, CoSpar and PAGA, side by side.

usage: python steps/05d_trees.py config/hdbi.yaml

similarity between the CARTA fates, by three methods:
  CARTA   share of the progenitors (tree at carta.k) holding one of the two fates that hold both (Jaccard)
  CoSpar  clonal fate coupling, all timepoints (step 5b), divided by its maximum
  PAGA    exp(-d), d = shortest path in the PAGA population graph with edge length -log(connectivity), between the
          populations of the two fates (closest pair)
tree: average linkage on 1 - similarity.
writes runs/<run_name>/support/trees.pdf, fate_similarity_<method>.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.sparse.csgraph import shortest_path
from scipy.spatial.distance import squareform

from common import load

cfg_path, repo, cfg, run = load(sys.argv[1])
cc = cfg['carta']
out = run / 'support'
short = cc.get('short', {})

pot = pd.read_csv(run / 'carta' / 'potencies_by_k.csv')
k = pd.read_csv(run / 'carta' / 'pareto.csv').query('used').k.iloc[0]
progs = [set(p.split('+')) for p in pot[pot.k == k].potency.unique()]
fates = list(pd.read_csv(run / 'carta' / 'clone_type_mat.csv', index_col=0).columns[1:-1])

carta = pd.DataFrame([[sum(f in s and g in s for s in progs) / max(sum(f in s or g in s for s in progs), 1)
                       for g in fates] for f in fates], index=fates, columns=fates)

cos = pd.read_csv(out / 'cospar_coupling_all.csv', index_col=0).loc[fates, fates]
c = cos.values.copy()
np.fill_diagonal(c, np.nan)
cos = pd.DataFrame(np.clip(cos.values / np.nanmax(c), 0, 1), index=fates, columns=fates)

paga = pd.read_csv(out / 'paga_connectivity.csv', index_col=0)
pops = list(paga.index)
with np.errstate(divide='ignore'):
    length = np.where(paga.values > 0, -np.log(paga.values), 0)
dist = shortest_path(np.ascontiguousarray(length, dtype=float), directed=False)
pops_of = {f: [p for p, cts in cfg['populations'].items() if p in pops and {cc['fates'].get(ct) for ct in cts} == {f}]
           for f in fates}
pg = pd.DataFrame([[np.exp(-min(dist[pops.index(a), pops.index(b)] for a in pops_of[f] for b in pops_of[g]))
                    for g in fates] for f in fates], index=fates, columns=fates)

col = {}
for f in fates:
    term = [p for p in cfg['terminal'] if p in pops_of[f]]
    col[f] = cfg.get('colors', {}).get(term[0] if term else pops_of[f][-1], 'k')

sims = {'CARTA (progenitor co-membership, k = %d)' % k: carta, 'CoSpar (clonal coupling / max)': cos,
        'PAGA (exp(-path length))': pg}
fig, axes = plt.subplots(2, 3, figsize=(15, 9.5), gridspec_kw={'height_ratios': [1, 0.8]})
for j, (title, s) in enumerate(sims.items()):
    tag = title.split()[0].lower()
    s.round(4).to_csv(out / f'fate_similarity_{tag}.csv')
    d = 1 - s.values
    np.fill_diagonal(d, 0)
    z = linkage(squareform((d + d.T) / 2, checks=False), method='average')
    ax = axes[1, j]
    dd = dendrogram(z, labels=[short.get(f, f) for f in fates], orientation='left', ax=ax, color_threshold=0,
                    above_threshold_color='0.25')
    for t in ax.get_yticklabels():
        t.set_color(col[{short.get(f, f): f for f in fates}[t.get_text()]])
        t.set_fontweight('bold')
    ax.set_xlabel('1 - similarity')
    ax.spines[['top', 'right', 'left']].set_visible(False)
    order = [{short.get(f, f): f for f in fates}[x] for x in dd['ivl']][::-1]
    ax = axes[0, j]
    v = s.loc[order, order].values.copy()
    np.fill_diagonal(v, np.nan)
    im = ax.imshow(v, cmap='magma_r', vmin=0, vmax=1)
    ax.set_xticks(range(len(order)), [short.get(f, f) for f in order], rotation=45, ha='right')
    ax.set_yticks(range(len(order)), [short.get(f, f) for f in order])
    for i in range(len(order)):
        for jj in range(len(order)):
            if np.isfinite(v[i, jj]):
                ax.text(jj, i, f'{v[i, jj]:.2f}', ha='center', va='center', fontsize=6.5,
                        color='white' if v[i, jj] > 0.55 else 'black')
    ax.set_title(title, fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.75)
fig.suptitle('Which fates are close: CARTA (clonal, tree), CoSpar (clonal), PAGA (expression)', fontsize=12)
fig.tight_layout()
fig.savefig(out / 'trees.pdf', metadata={'CreationDate': None})
print('wrote', out / 'trees.pdf')
for title, s in sims.items():
    v = s.where(np.triu(np.ones(s.shape, bool), 1)).stack().sort_values(ascending=False)
    print(f'{title}: closest pairs ' + ', '.join(f'{short.get(a, a)}-{short.get(b, b)} {x:.2f}' for (a, b), x in v.head(5).items()))
