"""Shared by the steps: config loading, run folder, cell selection.

A config may start with `base: <other config>` (path relative to it); its keys then override those of the base,
section by section (lists and values are replaced, dictionaries merged).
"""
from pathlib import Path

import yaml

missing_clone = ['nan', 'None', '']


def merge(base, over):
    out = dict(base)
    for k, v in over.items():
        out[k] = merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def read_config(path):
    path = Path(path).resolve()
    cfg = yaml.safe_load(open(path))
    if 'base' in cfg:
        cfg = merge(read_config(path.parent / cfg.pop('base')), cfg)
    return cfg


def load(path):
    """config path, repo folder, merged config, run folder"""
    path = Path(path).resolve()
    cfg = read_config(path)
    repo = path.parents[1]
    return path, repo, cfg, Path(repo / cfg['runs_dir'] / cfg['run_name'])


def save_config(cfg, path):
    with open(path, 'w') as f:
        yaml.safe_dump(cfg, f, sort_keys=False, allow_unicode=True)


def population_of(cfg):
    return {ct: p for p, cts in cfg['populations'].items() for ct in cts}


def selected(cells, cfg, barcoded=False, timepoints=True):
    """boolean mask: cells passing cell_filter, optionally with a clone and at a configured timepoint"""
    col = cfg['columns']
    keep = cells[col['celltype']].notna()
    for c, v in cfg.get('cell_filter', {}).items():
        keep &= cells[c].astype(str) == v
    if barcoded:
        keep &= ~cells[col['clone']].astype(str).isin(missing_clone)
    if timepoints:
        keep &= cells[col['timepoint']].astype(str).isin(cfg['timepoints'])
    return keep
