#!/usr/bin/env python3
"""Step 4: choose k and build the CARTA differentiation tree.

usage: python steps/04_carta_tree.py config/hdbi.yaml

reads runs/<run_name>/carta/sweep/ (step 3) and writes to runs/<run_name>/carta/:
  pareto.csv / pareto.pdf     discrepancy D(k) = sum over profiles of clones x fates predicted but not seen;
                              elbow (kneedle) marked, the k used for the tree circled
  potencies_by_k.csv          progenitor potencies and their clones at every k (are they stable as k grows?)
  tree_k<k>_progenitors.csv   progenitors of the tree: potency, parent, clones assigned
  tree_k<k>.pdf / .nwk        the tree: progenitors nested by potency, each fate under its smallest progenitor
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yaml

cfg_path = Path(sys.argv[1]).resolve()
repo = cfg_path.parents[1]
cfg = yaml.safe_load(open(cfg_path))
cc = cfg['carta']
carta = repo / cfg['runs_dir'] / cfg['run_name'] / 'carta'
sweep = carta / 'sweep'
no_date = {'CreationDate': None}

mat = pd.read_csv(carta / 'clone_type_mat.csv', index_col=0)
fates = list(mat.columns[1:-1])
inp = mat[(mat.counts >= cc['threshold']) & (mat[fates].sum(axis=1) > 1)]
ks = [k for k in range(1, cc['kmax'] + 1) if (sweep / f'k{k}_prog.csv').exists()]
assert ks == list(range(1, len(ks) + 1)) and ks, f'sweep incomplete: k found {ks}'


def short(p):
    return '+'.join(cc.get('short', {}).get(f, f) for f in p.split('+'))


def potency(row):
    return '+'.join(f for f in fates if row[f])


rows, dk = [], {}
for k in ks:
    lab = pd.read_csv(sweep / f'k{k}_labeling.csv', index_col=0).reindex(inp.index)
    extra = ((lab[fates].values == 1) & (inp[fates].values == 0)).sum(axis=1)
    dk[k] = float((inp.counts.values * extra).sum())
    prog = pd.read_csv(sweep / f'k{k}_prog.csv', index_col=0)
    asg = pd.read_csv(sweep / f'k{k}_assignment.csv', index_col=0).reindex(inp.index)
    for p, r in prog.iterrows():
        rows.append({'k': k, 'progenitor': p, 'potency': potency(r), 'n_fates': int(r[fates].sum()),
                     'clones': int(inp.counts[asg[str(p)] == 1].sum())})
pot = pd.DataFrame(rows)
pot.to_csv(carta / 'potencies_by_k.csv', index=False)

x = np.array(ks, float)
y = np.array([dk[k] for k in ks])
if len(ks) > 2:
    xn, yn = (x - x.min()) / np.ptp(x), (y - y.min()) / max(np.ptp(y), 1e-9)
    elbow = int(x[np.argmax((1 - xn) - yn)])
else:
    elbow = ks[-1]
k_use = elbow if cc['k'] == 'elbow' else int(cc['k'])
assert k_use in ks, f'k = {k_use} not in the sweep'
pd.DataFrame({'k': ks, 'discrepancy': y, 'elbow': [k == elbow for k in ks],
              'used': [k == k_use for k in ks]}).to_csv(carta / 'pareto.csv', index=False)

fig, ax = plt.subplots(figsize=(4.5, 3.2))
ax.plot(x, y, '-o', color='0.2', ms=4)
ax.plot(elbow, dk[elbow], 's', mfc='none', mec='tab:blue', ms=11, mew=1.5, label=f'elbow k = {elbow}')
ax.plot(k_use, dk[k_use], 'o', mfc='none', mec='tab:red', ms=15, mew=1.5, label=f'used k = {k_use}')
ax.set_xlabel('number of progenitors k')
ax.set_ylabel('discrepancy D(k)')
ax.set_xticks(ks)
ax.legend(frameon=False)
ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout()
fig.savefig(carta / 'pareto.pdf', metadata=no_date)
plt.close(fig)

# tree at k_use: progenitors with identical potency merged, parent = smallest strict superset
tk = pot[pot.k == k_use].groupby('potency', as_index=False).agg(clones=('clones', 'sum'), n_fates=('n_fates', 'first'))
sets = {r.potency: set(r.potency.split('+')) for r in tk.itertuples()}
parent = {}
for a, sa in sets.items():
    sup = [b for b, sb in sets.items() if sb > sa]
    parent[a] = min(sup, key=lambda b: len(sets[b])) if sup else None
tk['parent'] = tk.potency.map(parent)
leaf_of = {f: min([p for p, s in sets.items() if f in s], key=lambda p: len(sets[p])) for f in fates
           if any(f in s for s in sets.values())}
tk = tk.sort_values('n_fates', ascending=False).reset_index(drop=True)
tk.insert(0, 'node', [f'P{i}' for i in range(len(tk))])
node = dict(zip(tk.potency, tk.node))
tk['parent'] = tk.parent.map(node)
tk['fates_attached'] = tk.potency.map(lambda p: '+'.join(f for f, q in leaf_of.items() if q == p))
tk.to_csv(carta / f'tree_k{k_use}_progenitors.csv', index=False)

kids = {n: [] for n in tk.node}
for r in tk.itertuples():
    if isinstance(r.parent, str):
        kids[r.parent].append(r.node)
for f, p in leaf_of.items():
    kids[node[p]].append(f)
roots = tk.node[tk.parent.isna()].tolist()


def newick(n):
    return n if n not in kids or not kids[n] else '(' + ','.join(newick(c) for c in kids[n]) + ')' + n


(carta / f'tree_k{k_use}.nwk').write_text(';\n'.join(newick(r) for r in roots) + ';\n')

order, depth = [], {}


def walk(n, d):
    depth[n] = d
    for c in sorted(kids.get(n, []), key=lambda c: (c in kids, c)):
        walk(c, d + 1)
    if n not in kids:
        order.append(n)


for r in roots:
    walk(r, 0)
ypos = {f: i for i, f in enumerate(order)}


def place(n):
    if n in ypos:
        return ypos[n]
    ypos[n] = float(np.mean([place(c) for c in kids[n]]))
    return ypos[n]


for r in roots:
    place(r)
xmax = max(depth.values())
info = tk.set_index('node')
fig, ax = plt.subplots(figsize=(7, 0.45 * len(order) + 1))
for n, cs in kids.items():
    for c in cs:
        xc = xmax if c not in kids else depth[c]
        ax.plot([depth[n], depth[n], xc], [ypos[n], ypos[c], ypos[c]], color='0.3', lw=1)
for n in kids:
    ax.plot(depth[n], ypos[n], 'o', color='0.2', ms=5, zorder=3)
    ax.annotate(f"{short(info.potency[n])}\n{n}, {info.clones[n]:,} clones", (depth[n], ypos[n]),
                xytext=(4, 4), textcoords='offset points', fontsize=6.5, va='bottom')
for f in order:
    ax.text(xmax + 0.08, ypos[f], f, va='center', fontsize=9)
ax.set_xlim(-0.2, xmax + 1.2)
ax.set_axis_off()
ax.set_title(f'CARTA tree, k = {k_use} (threshold {cc["threshold"]}, {int(inp.counts.sum()):,} clones)', fontsize=9)
fig.tight_layout()
fig.savefig(carta / f'tree_k{k_use}.pdf', metadata=no_date)
plt.close(fig)

pd.set_option('display.width', 200)
pd.set_option('display.max_colwidth', None)
print(pd.Series(dk, name='D(k)').to_string())
print(f'\nelbow k = {elbow}; tree built at k = {k_use} (config carta.k = {cc["k"]})\n')
print(tk.to_string(index=False))
print('\npotencies by k:')
for k, s in pot.groupby('k').potency:
    print(f'k = {k}: ' + ' ; '.join(short(p) for p in sorted(s, key=len)))
