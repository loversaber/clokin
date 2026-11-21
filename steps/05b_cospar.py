#!/usr/bin/env python3
"""Step 5b: CoSpar clonal fate coupling between the CARTA fates.

usage: python steps/05b_cospar.py config/hdbi.yaml

part 1 (this python): barcoded cells kept by `cell_filter`; state = CARTA fate (carta.fates), or the population for
Undifferentiated cells -> runs/<run_name>/support/cospar_input.npz
part 2 (support.cospar_python, started automatically): coupling from shared clones (source X_clone), for all
timepoints and each of support.cospar_times -> support/cospar_coupling_<time>.csv (fates x fates), cospar_coupling.pdf
"""
import json
import os
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings('ignore')
cfg_path = Path(sys.argv[1]).resolve()
repo = cfg_path.parents[1]

if '--json' not in sys.argv:
    import pandas as pd
    from common import load, population_of, selected
    _, _, cfg, run = load(cfg_path)
    col, sp, cc = cfg['columns'], cfg['support'], cfg['carta']
    out = run / 'support'
    out.mkdir(parents=True, exist_ok=True)

    pop_of = population_of(cfg)
    cells = pd.read_parquet(repo / cfg['inputs']['cells'])
    pca = pd.read_parquet(repo / sp['pca'])
    assert pca.index.equals(cells.index), 'pca and cell table must have the same cells in the same order'
    ct = cells[col['celltype']].astype(str)
    keep = selected(cells, cfg, barcoded=True) & ct.isin(cc['fates'])
    fate = ct[keep].map(cc['fates'])
    clones = pd.Categorical(cells.loc[keep, col['clone']].astype(str))
    np.savez(out / 'cospar_input.npz',
             cell=cells.index[keep].values.astype(str),
             pca=pca.loc[keep].values.astype(np.float32),
             emb=np.asarray(cells.loc[keep, ['UMAP1', 'UMAP2']].values if {'UMAP1', 'UMAP2'} <= set(cells.columns)
                            else pca.loc[keep].values[:, :2], dtype=np.float32),
             clone=clones.codes, n_clones=len(clones.categories),
             time=np.asarray(cells.loc[keep, col['timepoint']].astype(str).map(cfg['timepoints']).astype(str), dtype=str),
             state=np.asarray(np.where(fate == 'Undifferentiated', ct[keep].map(pop_of), fate), dtype=str),
             fates=np.array([f for f in fate.value_counts().index if f != 'Undifferentiated'], dtype=str))
    js = out / 'config.json'
    json.dump(cfg, open(js, 'w'))
    print(f'{keep.sum():,} cells, {len(clones.categories):,} clones -> {out / "cospar_input.npz"}', flush=True)
    py = sp['cospar_python']
    os.execv(py, [py, '-u', os.path.abspath(__file__), str(cfg_path), '--json', str(js)])

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
import cospar as cs

cfg = json.load(open(sys.argv[sys.argv.index('--json') + 1]))
sp, cc = cfg['support'], cfg['carta']
out = repo / cfg['runs_dir'] / cfg['run_name'] / 'support'
cs.settings.data_path = str(out / 'cospar_tmp')
cs.settings.figure_path = str(out / 'cospar_tmp')
cs.settings.verbosity = 0

d = np.load(out / 'cospar_input.npz', allow_pickle=False)
fates = list(d['fates'])
n = len(d['cell'])
x_clone = csr_matrix((np.ones(n), (np.arange(n), d['clone'])), shape=(n, int(d['n_clones'])))
adata = cs.pp.initialize_adata_object(
    X_state=d['pca'], cell_names=d['cell'], gene_names=np.array([f'pc{i + 1}' for i in range(d['pca'].shape[1])]),
    time_info=d['time'], X_clone=x_clone, X_pca=d['pca'], X_emb=d['emb'], state_info=d['state'], data_des='clokin')
cs.hf.update_time_ordering(adata, updated_ordering=[str(v) for v in cfg['timepoints'].values()])

res = {}
for t in sp['cospar_times']:
    sel = None if t == 'all' else [str(cfg['timepoints'][t])]
    cs.tl.fate_coupling(adata, selected_fates=fates, source='X_clone', selected_times=sel)
    m = pd.DataFrame(adata.uns['fate_coupling_X_clone']['X_coupling'], index=fates, columns=fates)
    m.to_csv(out / f'cospar_coupling_{t}.csv')
    res[t] = m
    print(f'\n{t}:\n{m.round(3).to_string()}', flush=True)

short = cc.get('short', {})
lab = [short.get(f, f) for f in fates]
fig, axes = plt.subplots(1, len(res), figsize=(3.6 * len(res), 3.4))
for ax, (t, m) in zip(np.atleast_1d(axes), res.items()):
    v = m.values.copy()
    np.fill_diagonal(v, np.nan)
    im = ax.imshow(v, cmap='magma_r', vmin=0)
    ax.set_xticks(range(len(fates)))
    ax.set_xticklabels(lab, rotation=45, ha='right', fontsize=7)
    ax.set_yticks(range(len(fates)))
    ax.set_yticklabels(lab, fontsize=7)
    for i in range(len(fates)):
        for j in range(len(fates)):
            if np.isfinite(v[i, j]):
                ax.text(j, i, f'{v[i, j]:.2f}', ha='center', va='center', fontsize=5,
                        color='white' if v[i, j] > 0.5 * np.nanmax(v) else 'black')
    ax.set_title(f'CoSpar clonal coupling, {t}', fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.7)
fig.tight_layout()
fig.savefig(out / 'cospar_coupling.pdf', metadata={'CreationDate': None})
