#!/usr/bin/env python3
"""Step 9: figures from the tables of steps 1, 5, 6 and 8.

usage: python steps/09_figures.py config/hdbi.yaml

writes runs/<run_name>/figures/ (PDF):
  graph_on_umap        transition graph used by CLADES at the PAGA node positions; edge colour = evidence status
  fit_overview         log2(predicted / observed) cells per row and timepoint, trust label
  fit_scatter          predicted vs observed cells, every row x population x timepoint
  trajectories_<row>   observed points and predicted curves per population (one page per row)
  composition_last     observed vs predicted composition at the last timepoint
  rates_<day>          differentiation rate of every graph edge and net growth of every population, per row
  fate_prob            Gillespie P(reach) of the terminal populations, per row ('do not quote' rows greyed)
  divisions            Gillespie mean divisions from the root to each terminal population
  convergence          loss terms over the epochs
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from common import load, population_of

cfg_path, repo, cfg, run = load(sys.argv[1])
dg, out = run / 'diagnostics', run / 'figures'
out.mkdir(parents=True, exist_ok=True)
no_date = {'CreationDate': None}
pops = list(cfg['populations'])
colors = cfg.get('colors', {})
term = cfg['terminal']
tps = list(cfg['timepoints'])
days = list(cfg['timepoints'].values())
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})


def save(fig, name):
    fig.savefig(out / f'{name}.pdf', metadata=no_date, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('wrote', name, flush=True)


fit = pd.read_csv(dg / 'fit_by_row.csv')
rows = list(fit.row)
trust = dict(zip(fit.row, fit.trust))
lab = {r: f"{r}{'' if trust[r] == 'use' else ' *'}" for r in rows}

# transition graph on the PAGA layout
pos = pd.read_csv(run / 'support' / 'paga_node_positions.csv', index_col=0)
edges = pd.read_csv(run / 'support' / 'edges_used.csv')
cells = pd.read_parquet(repo / cfg['inputs']['cells'], columns=['UMAP1', 'UMAP2', cfg['columns']['celltype']])
pop_of = population_of(cfg)
cp = cells[cfg['columns']['celltype']].astype(str).map(pop_of)
cells = cells[cp.notna()]
sc_col = cp[cp.notna()].map(lambda p: colors.get(p, '0.8'))
stat_col = {'agree': '#1b7837', 'paga only': '#e08214', 'no clonal data': '#4575b4', 'not a candidate': '0.4'}
fig, ax = plt.subplots(figsize=(8, 7))
ax.scatter(cells.UMAP1, cells.UMAP2, c=sc_col.values, s=0.1, lw=0, alpha=0.25, rasterized=True)
for e in edges[edges.use == 1].itertuples():
    (x0, y0), (x1, y1) = pos.loc[e.source].values, pos.loc[e.target].values
    ax.annotate('', (x1, y1), (x0, y0), arrowprops=dict(arrowstyle='-|>', lw=1.6, mutation_scale=12, shrinkA=7, shrinkB=7,
                                                         color=stat_col.get(e.status, 'k')), zorder=3)
ax.scatter(pos.UMAP1 if 'UMAP1' in pos else pos.iloc[:, 0], pos.iloc[:, 1], s=60,
           c=[colors.get(p, '0.5') for p in pos.index], edgecolors='k', lw=0.6, zorder=4)
for p, (x, y) in pos.iterrows():
    ax.annotate(p, (x, y), xytext=(0, 7), textcoords='offset points', ha='center', fontsize=6.5, zorder=5)
for s, c in stat_col.items():
    if s in set(edges[edges.use == 1].status):
        ax.plot([], [], color=c, lw=2, label=s)
ax.legend(title='edge evidence', frameon=False, loc='lower right')
ax.set_xticks([])
ax.set_yticks([])
ax.set_title(f'transition graph used by CLADES ({int(edges.use.sum())} edges)')
save(fig, 'graph_on_umap')

# fit overview
cols = [f'log2_pred_obs_{t}' for t in tps]
v = fit.set_index('row')[cols].values
fig, ax = plt.subplots(figsize=(3.2, 0.25 * len(rows) + 1))
lim = max(2, np.nanmax(np.abs(v)))
im = ax.imshow(v, cmap='RdBu_r', vmin=-lim, vmax=lim, aspect='auto')
for i in range(v.shape[0]):
    for j in range(v.shape[1]):
        ax.text(j, i, f'{v[i, j]:+.1f}', ha='center', va='center', fontsize=6)
ax.set_xticks(range(len(tps)), tps)
ax.set_yticks(range(len(rows)), [lab[r] for r in rows])
fig.colorbar(im, ax=ax, label='log2(predicted / observed cells)', shrink=0.6)
ax.set_title('fit per row (* = do not quote)')
save(fig, 'fit_overview')

# trajectories + scatter
data = pd.read_csv(run / 'data' / 'counts_long.csv')
traj = pd.read_csv(dg / 'trajectories.csv')
d0 = pd.read_csv(run / 'data' / 'initial_condition.csv', index_col=0)['leiden'].values
fig, axes = plt.subplots(1, len(tps), figsize=(3.2 * len(tps), 3.2))
for ax, t, day in zip(axes, tps, days):
    o = data[data.tp == t].pivot_table(index='row', columns='pop', values='n_scaled', aggfunc='sum').reindex(index=rows, columns=pops).fillna(0)
    p = traj[traj.day == day].pivot_table(index='row', columns='population', values='predicted').reindex(index=rows, columns=pops)
    ok = (o.values > 0) | (p.values > 1)
    tr = np.array([[trust[r] == 'use'] * len(pops) for r in rows])
    for sel, c, l in [(tr & ok, 'k', 'use'), (~tr & ok, '0.7', 'do not quote')]:
        ax.scatter(o.values[sel] + 1, p.values[sel] + 1, s=4, c=c, lw=0, label=l)
    hi = max(o.values.max(), np.nanmax(p.values)) * 2
    ax.plot([1, hi], [1, hi], lw=0.6, color='tab:red')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('observed cells + 1')
    ax.set_ylabel('predicted cells + 1')
    ax.set_title(t)
axes[0].legend(frameon=False)
save(fig, 'fit_scatter')

with PdfPages(out / 'trajectories.pdf', metadata=no_date) as pdf:
    for i, r in enumerate(rows):
        fig, axes = plt.subplots(3, 7, figsize=(14, 6), sharex=True)
        for ax, p in zip(axes.ravel(), pops):
            tt = traj[(traj.row == r) & (traj.population == p)]
            ax.plot(tt.day, tt.predicted, color=colors.get(p, 'k'), lw=1.2)
            ob = [d0[i] if p == pops[0] else 0] + [data[(data.row == r) & (data.tp == t) & (data['pop'] == p)].n_scaled.sum() for t in tps]
            ax.plot([0] + days, ob, 'o', color='k', ms=3)
            ax.set_title(p, fontsize=7)
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            ax.set_xticks([0] + days)
        fig.suptitle(f'{r} ({trust[r]}): observed (points) and predicted (line) cells', fontsize=10)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
print('wrote trajectories (one page per row)')

# composition at the last timepoint
comp = pd.read_csv(dg / 'composition_last.csv')
fig, ax = plt.subplots(figsize=(10, 0.32 * len(rows) * 2 + 1))
y, ticks, tl = 0, [], []
for r in rows:
    for kind in ['observed', 'predicted']:
        left = 0
        for p in pops:
            w = comp[(comp.row == r) & (comp.population == p)][kind].iloc[0]
            ax.barh(y, w, left=left, color=colors.get(p, '0.5'), height=0.8, lw=0)
            left += w
        ticks.append(y)
        tl.append(f'{lab[r]} {kind[:4]}.')
        y += 1
    y += 0.6
ax.set_yticks(ticks, tl, fontsize=6)
ax.invert_yaxis()
ax.set_xlabel(f'share of cells at {tps[-1]}')
for p in pops:
    ax.barh([], [], color=colors.get(p, '0.5'), label=p)
ax.legend(ncol=1, fontsize=6, frameon=False, bbox_to_anchor=(1.01, 1), loc='upper left')
save(fig, 'composition_last')

# rates
rates = pd.read_csv(dg / 'rates.csv')
for day in days:
    rd = rates[rates.day == day]
    dif = rd[rd.kind == 'differentiation'].assign(edge=lambda x: x.source + ' > ' + x.target)
    grw = rd[rd.kind == 'net growth']
    a = dif.pivot_table(index='row', columns='edge', values='rate').reindex(index=rows)
    b = grw.pivot_table(index='row', columns='source', values='rate').reindex(index=rows, columns=pops)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 0.28 * len(rows) + 3), gridspec_kw={'width_ratios': [a.shape[1], b.shape[1]]})
    im1 = ax1.imshow(a.values, cmap='viridis', aspect='auto', vmin=0)
    ax1.set_xticks(range(a.shape[1]), a.columns, rotation=60, ha='right', fontsize=6)
    ax1.set_yticks(range(len(rows)), [lab[r] for r in rows])
    ax1.set_title(f'differentiation rate (per day), day {day}')
    fig.colorbar(im1, ax=ax1, shrink=0.6)
    lim = np.nanmax(np.abs(b.values))
    im2 = ax2.imshow(b.values, cmap='RdBu_r', aspect='auto', vmin=-lim, vmax=lim)
    ax2.set_xticks(range(len(pops)), pops, rotation=60, ha='right', fontsize=6)
    ax2.set_yticks(range(len(rows)), [lab[r] for r in rows])
    ax2.set_title(f'net growth (per day), day {day}')
    fig.colorbar(im2, ax=ax2, shrink=0.6)
    save(fig, f'rates_day{day}')

# Gillespie
for name, f, fmt, cmap, vmax in [('fate_prob', 'fate_prob.csv', '{:.2f}', 'magma_r', 1), ('divisions', 'divisions.csv', '{:.1f}', 'viridis', None)]:
    if not (dg / f).exists():
        print(f'no {f} yet (Gillespie not complete)')
        continue
    t = pd.read_csv(dg / f).set_index('row').reindex(rows)
    v = t[term].values.astype(float)
    fig, ax = plt.subplots(figsize=(0.75 * len(term) + 2.5, 0.3 * len(rows) + 1.5))
    im = ax.imshow(v, cmap=cmap, aspect='auto', vmin=0, vmax=vmax)
    for i, r in enumerate(rows):
        for j in range(len(term)):
            if np.isfinite(v[i, j]):
                ax.text(j, i, fmt.format(v[i, j]), ha='center', va='center', fontsize=6,
                        color=('white' if v[i, j] > 0.55 * np.nanmax(v) else 'black'))
        if trust[r] != 'use':
            ax.add_patch(plt.Rectangle((-0.5, i - 0.5), len(term), 1, color='white', alpha=0.6, lw=0))
    ax.set_xticks(range(len(term)), term, rotation=40, ha='right')
    ax.set_yticks(range(len(rows)), [lab[r] for r in rows])
    ax.set_title('P(a lineage from one root cell reaches the population)' if name == 'fate_prob'
                 else 'mean divisions from the root', fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.6)
    save(fig, name)

conv = pd.read_csv(dg / 'convergence.csv') if (dg / 'convergence.csv').stat().st_size > 1 else pd.DataFrame()
if len(conv):
    obs = [c for c in conv.columns if c.startswith('L') and c[1:].isdigit()]
    fig, ax = plt.subplots(figsize=(4, 3))
    for c in obs:
        ax.plot(conv.epoch, conv[c], label=f'day {c[1:]}')
    ax.set_xlabel('epoch')
    ax.set_ylabel('loss')
    ax.legend(frameon=False)
    save(fig, 'convergence')
