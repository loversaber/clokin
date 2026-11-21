#!/usr/bin/env python3
"""Step 5a: PAGA connectivity between the populations (expression only, no clonal information).

usage: python steps/05a_paga.py config/hdbi.yaml

cells: those kept by `cell_filter` with a population, all timepoints and experiments. Neighbours on `support.pca`.
writes runs/<run_name>/support/paga_connectivity.csv (populations x populations) and paga_connectivity.pdf
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc

from common import load, population_of, selected

cfg_path, repo, cfg, run = load(sys.argv[1])
col, sp = cfg['columns'], cfg['support']
out = run / 'support'
out.mkdir(parents=True, exist_ok=True)

pops = list(cfg['populations'])
pop_of = population_of(cfg)
cells = pd.read_parquet(repo / cfg['inputs']['cells'])
pca = pd.read_parquet(repo / sp['pca'])
assert pca.index.equals(cells.index), 'pca and cell table must have the same cells in the same order'

keep = selected(cells, cfg, timepoints=False) & cells[col['celltype']].astype(str).map(pop_of).notna()
a = ad.AnnData(obs=pd.DataFrame({'pop': pd.Categorical(cells.loc[keep, col['celltype']].astype(str).map(pop_of),
                                                       categories=pops)}, index=cells.index[keep]))
a.obsm['X_pca'] = pca.loc[keep].values
print(f'{a.n_obs:,} cells', flush=True)

sc.pp.neighbors(a, n_neighbors=sp['paga_neighbors'], use_rep='X_pca', random_state=0)
sc.tl.paga(a, groups='pop')
conn = pd.DataFrame(a.uns['paga']['connectivities'].toarray(), index=pops, columns=pops)
conn.to_csv(out / 'paga_connectivity.csv')

c = conn.values.copy()
np.fill_diagonal(c, np.nan)
fig, ax = plt.subplots(figsize=(9, 8))
im = ax.imshow(c, cmap='magma_r', vmin=0)
ax.set_xticks(range(len(pops)), pops, rotation=45, ha='right', fontsize=7)
ax.set_yticks(range(len(pops)), pops, fontsize=7)
for i in range(len(pops)):
    for j in range(len(pops)):
        if np.isfinite(c[i, j]) and c[i, j] >= 0.05:
            ax.text(j, i, f'{c[i, j]:.2f}', ha='center', va='center', fontsize=5,
                    color='white' if c[i, j] > 0.55 * np.nanmax(c) else 'black')
fig.colorbar(im, ax=ax, shrink=0.7, label='PAGA connectivity')
ax.set_title(f'PAGA connectivity ({a.n_obs:,} cells, {sp["paga_neighbors"]} neighbours)', fontsize=9)
fig.tight_layout()
fig.savefig(out / 'paga_connectivity.pdf', metadata={'CreationDate': None})
plt.close(fig)

umap = cells.loc[keep, ['UMAP1', 'UMAP2']].values
pos = pd.DataFrame(umap, index=a.obs_names).groupby(a.obs['pop'].values, observed=True).median().reindex(pops)
pos.to_csv(out / 'paga_node_positions.csv')
size = a.obs['pop'].value_counts().reindex(pops)
colors = [cfg.get('colors', {}).get(p, '0.5') for p in pops]
cell_col = a.obs['pop'].map(dict(zip(pops, colors))).astype(object).values
order = np.random.default_rng(0).permutation(a.n_obs)
thr = sp.get('paga_threshold', 0.1)
fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
for ax, show_cells in zip(axes, [True, False]):
    if show_cells:
        ax.scatter(umap[order, 0], umap[order, 1], c=cell_col[order], s=0.2, lw=0, rasterized=True)
    for i in range(len(pops)):
        for j in range(i + 1, len(pops)):
            w = conn.values[i, j]
            if w >= thr:
                ax.plot(pos.iloc[[i, j], 0], pos.iloc[[i, j], 1], color='0.15', lw=0.5 + 5 * w, alpha=0.75, zorder=2,
                        solid_capstyle='round')
    node_s = 40 + 600 * np.sqrt(size / size.max())
    ax.scatter(pos.iloc[:, 0], pos.iloc[:, 1], s=node_s, c=colors, edgecolors='0.15', lw=0.8, zorder=3)
    if not show_cells:
        for p, (x, y) in pos.iterrows():
            ax.annotate(p, (x, y), xytext=(0, 3 + np.sqrt(node_s[p]) / 2), textcoords='offset points', ha='center',
                        fontsize=7, zorder=4)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.spines[['top', 'right']].set_visible(False)
axes[0].set_title('cells on the UMAP, PAGA graph on top', fontsize=9)
axes[1].set_title(f'PAGA graph (edges: connectivity >= {thr}; width = connectivity; node size = cells)', fontsize=9)
for w in [0.25, 0.5, 1.0]:
    axes[1].plot([], [], color='0.15', lw=0.5 + 5 * w, alpha=0.75, label=f'{w}')
axes[1].legend(title='connectivity', frameon=False, fontsize=7, title_fontsize=7, loc='lower right')
fig.tight_layout()
fig.savefig(out / 'paga_umap.pdf', metadata={'CreationDate': None}, dpi=300)
plt.close(fig)

top = conn.where(np.triu(np.ones(conn.shape, bool), 1)).stack().sort_values(ascending=False)
print('strongest links:')
print(top[top >= 0.15].round(3).to_string())
