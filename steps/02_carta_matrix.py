#!/usr/bin/env python3
"""Step 2: clone x fate matrix for CARTA.

usage: python steps/02_carta_matrix.py config/hdbi.yaml

writes runs/<run_name>/carta/:
  clone_fate_counts.csv   one row per clone: cells of each fate (and Undifferentiated)
  clone_type_mat.csv      unique fate-presence profiles + number of clones (larry.py -i)
  profiles.csv            the multi-fate profiles that enter the ILP at the threshold
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from common import load, selected

cfg_path, repo, cfg, run = load(sys.argv[1])
col, cc = cfg['columns'], cfg['carta']
out = run / 'carta'
out.mkdir(parents=True, exist_ok=True)

cells = pd.read_parquet(repo / cfg['inputs']['cells'])
keep = selected(cells, cfg, barcoded=True)
df = pd.DataFrame({'clone': cells.loc[keep, col['clone']].astype(str),
                   'fate': cells.loc[keep, col['celltype']].astype(str).map(cc['fates'])})
missing = sorted(set(cells.loc[keep, col['celltype']].astype(str)) - set(cc['fates']))
assert not missing, f'cell types without a fate: {missing}'
if cc['clones'] == 'labelled':
    labelled = set(pd.read_csv(repo / cfg['inputs']['labels'])[col['label_clone']].astype(str))
    df = df[df.clone.isin(labelled)]

undiff = 'Undifferentiated'
fates = [f for f in df.fate.value_counts().index if f != undiff]
counts = df.pivot_table(index='clone', columns='fate', aggfunc='size', fill_value=0)
counts = counts.reindex(columns=[undiff] + fates, fill_value=0)
counts.to_csv(out / 'clone_fate_counts.csv')

pres = (counts[fates] >= cc['min_cells_per_fate']).astype(int)
pres.insert(0, undiff, (counts[undiff] > 0).astype(int))
uniq, inv = np.unique(pres.values, axis=0, return_inverse=True)
mat = pd.DataFrame(uniq, columns=pres.columns)
mat['counts'] = np.bincount(inv.ravel(), minlength=len(uniq))
mat.to_csv(out / 'clone_type_mat.csv')

sel = mat[(mat.counts >= cc['threshold']) & (mat[fates].sum(axis=1) > 1)]
prof = pd.DataFrame({'profile': ['+'.join(f for f in fates if r[f]) for _, r in sel.iterrows()],
                     'clones': sel.counts.values}, index=sel.index).sort_values('clones', ascending=False)
prof.to_csv(out / 'profiles.csv')

print(f'{len(df):,} cells, {len(counts):,} clones -> {len(mat)} fate profiles over {len(fates)} fates: {fates}')
print(f'at threshold {cc["threshold"]}: {len(sel)} multi-fate profiles, {sel.counts.sum():,} clones enter the ILP')
print(prof.head(10).to_string())
