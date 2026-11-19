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
import yaml

cfg_path = Path(sys.argv[1]).resolve()
repo = cfg_path.parents[1]
cfg = yaml.safe_load(open(cfg_path))
col, sp = cfg['columns'], cfg['support']
out = repo / cfg['runs_dir'] / cfg['run_name'] / 'support'
out.mkdir(parents=True, exist_ok=True)

pops = list(cfg['populations'])
pop_of = {ct: p for p, cts in cfg['populations'].items() for ct in cts}
cells = pd.read_parquet(repo / cfg['inputs']['cells'])
pca = pd.read_parquet(repo / sp['pca'])
assert pca.index.equals(cells.index), 'pca and cell table must have the same cells in the same order'

keep = cells[col['celltype']].astype(str).map(pop_of).notna()
for c, v in cfg.get('cell_filter', {}).items():
    keep &= cells[c].astype(str) == v
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

top = conn.where(np.triu(np.ones(conn.shape, bool), 1)).stack().sort_values(ascending=False)
print('strongest links:')
print(top[top >= 0.15].round(3).to_string())
