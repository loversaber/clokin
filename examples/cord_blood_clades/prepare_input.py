"""Make the clokin inputs from the CLADES cord blood data (figshare 10.6084/m9.figshare.27908142, CordBlood_Refine).

    python prepare_input.py <folder with the cord blood AnnData and CordBlood_Obs_Clones.pkl> <clonaltrans folder>

writes inputs/ next to this script:
  cells.parquet         one row per cell: clone (empty if no barcode), cell type (def_lab), timepoint, experiment,
                        GFP, UMAP1, UMAP2
  pca.parquet           X_pca of the same cells (50 components), for PAGA and CoSpar
  metaclone_labels.csv  CLADES meta-clone (0-11) of each clone: GFP+ barcoded clones with at least 3 cells, as in
                        CLADES (demo/CordBlood_Refine/prepare_input.ipynb)
  scaling_factors.csv   sequenced -> absolute cell number per timepoint, the factors of that notebook
  clades_graph_table.csv  the CLADES transition graph (demo/CordBlood_Refine/graph_table.csv, MIT licence)
"""
import shutil
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

src, clonaltrans = Path(sys.argv[1]), Path(sys.argv[2])
out = Path(__file__).parent / 'inputs'
out.mkdir(exist_ok=True)

a = ad.read_h5ad(src / 'adata_cord_blood_as_submitted.h5ad', backed='r')
obs = a.obs
clone = obs.clones.astype(str).where(obs.clones.astype(str) != 'Clone_')
cells = pd.DataFrame({'clone': clone, 'celltype': obs.def_lab.astype(str), 'timepoint': obs.Timepoint.astype(str),
                      'experiment': obs.Tissue.astype(str), 'GFP': obs.GFP.astype(str),
                      'UMAP1': a.obsm['X_umap'][:, 0], 'UMAP2': a.obsm['X_umap'][:, 1]}, index=obs.index)
cells.to_parquet(out / 'cells.parquet')
pd.DataFrame(np.asarray(a.obsm['X_pca'], dtype=np.float32), index=obs.index,
             columns=[f'pc{i + 1}' for i in range(a.obsm['X_pca'].shape[1])]).to_parquet(out / 'pca.parquet')

size = cells[(cells.GFP == 'GFP+') & cells.clone.notna()].groupby('clone').size()
kept = size[size >= 3].index
mc = obs.loc[obs['Meta clones'].astype(str) != '-1.0'].groupby('clones', observed=True)['Meta clones'].first()
mc = mc.astype(str).str.replace(r'\.0$', '', regex=True).reindex(kept)
assert mc.notna().all()
mc.rename('meta_clone').rename_axis('clone_id').to_csv(out / 'metaclone_labels.csv')

d3, d10 = 13284 / 2708, 500000 * (13284 / 4428) / 2827
d17 = 1764225 * (500000 / 460000) * (13284 / 4428) / 11310
pd.DataFrame({'Day3': [d3], 'Day10': [d10], 'Day17': [d17]},
             index=pd.Index([cells.experiment.iloc[0]], name='experiment')).to_csv(out / 'scaling_factors.csv')

shutil.copy(clonaltrans / 'demo' / 'CordBlood_Refine' / 'graph_table.csv', out / 'clades_graph_table.csv')
print(f'{len(cells):,} cells ({cells.clone.notna().sum():,} barcoded), {len(mc)} clones in {mc.nunique()} meta-clones')
