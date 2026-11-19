#!/usr/bin/env python3
"""Step 6: transition graph for CLADES from the `graph` section of the config.

usage: python steps/06_graph.py config/hdbi.yaml

writes runs/<run_name>/data/graph_table.csv (populations x populations, 1 = source -> target) and
support/edges_used.csv (support/edges.csv with `use` filled in, plus graph edges that were not candidates).
Every population other than the root must be reachable from the root.
"""
import sys
from pathlib import Path

import pandas as pd
import yaml

cfg_path = Path(sys.argv[1]).resolve()
repo = cfg_path.parents[1]
cfg = yaml.safe_load(open(cfg_path))
run = repo / cfg['runs_dir'] / cfg['run_name']
pops = list(cfg['populations'])

g = pd.DataFrame(0, index=pops, columns=pops, dtype=int)
for a, bs in cfg['graph'].items():
    for b in bs:
        assert a in pops and b in pops, f'unknown population in graph: {a} -> {b}'
        assert a != b, f'self edge {a}'
        g.loc[a, b] = 1

seen, todo = {pops[0]}, [pops[0]]
while todo:
    a = todo.pop()
    for b in g.columns[g.loc[a] == 1]:
        if b not in seen:
            seen.add(b)
            todo.append(b)
unreached = [p for p in pops if p not in seen]
assert not unreached, f'not reachable from {pops[0]}: {unreached}'
g.to_csv(run / 'data' / 'graph_table.csv')

edges = pd.read_csv(run / 'support' / 'edges.csv')
used = {(a, b) for a in pops for b in pops if g.loc[a, b]}
edges['use'] = [int((a, b) in used) for a, b in zip(edges.source, edges.target)]
extra = used - set(zip(edges.source, edges.target))
if extra:
    edges = pd.concat([edges, pd.DataFrame([{'source': a, 'target': b, 'status': 'not a candidate', 'use': 1}
                                            for a, b in sorted(extra)])], ignore_index=True)
edges.to_csv(run / 'support' / 'edges_used.csv', index=False)

pd.set_option('display.width', 200)
print(f'{g.values.sum()} edges, all {len(pops)} populations reachable from {pops[0]}')
print(edges.groupby(['status', 'use']).size().unstack(fill_value=0).to_string())
weak = edges[(edges.use == 1) & ~edges.status.isin(['agree'])]
print('\nused edges without clonal + PAGA agreement:')
print(weak[['source', 'target', 'status', 'clonal_share', 'paga']].to_string(index=False))
