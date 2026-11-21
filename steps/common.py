"""Shared by the steps: config loading, run folder, cell selection.

A config may start with `base: <other config>` (path relative to it); its keys then override those of the base,
section by section (lists and values are replaced, dictionaries merged).
"""
from pathlib import Path

import pandas as pd
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
    """config path, root folder for the relative paths in the config, merged config, run folder

    root: the folder above the config's folder (config/<dataset>.yaml -> the clokin folder), or, with `root: <path>`
    in the config, that path relative to the config's folder (`root: .` for a self-contained example folder)
    """
    path = Path(path).resolve()
    cfg = read_config(path)
    root = (path.parent / cfg['root']).resolve() if 'root' in cfg else path.parents[1]
    return path, root, cfg, Path(root / cfg['runs_dir'] / cfg['run_name'])


def save_config(cfg, path):
    with open(path, 'w') as f:
        yaml.safe_dump(cfg, f, sort_keys=False, allow_unicode=True)


def read_table(path):
    path = Path(path)
    return pd.read_pickle(path) if path.suffix == '.pkl' else pd.read_csv(path, index_col=0)


def umap_columns(cfg):
    return list(cfg['columns'].get('umap', ['UMAP1', 'UMAP2']))


def clades_code(cfg):
    """CLADES code folder for clades.version ('main': the code before v1.2.0, '1.2': release v1.2.0); a site config
    can give one `code_dir` or one per version in `code_dirs`"""
    cl = cfg['clades']
    return cl['code_dirs'][clades_version(cfg)] if 'code_dirs' in cl else cl['code_dir']


def clades_version(cfg):
    return str(cfg.get('clades', {}).get('version', 'main'))


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
