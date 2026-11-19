# clokin

**CLOnes to KINetics.** From lineage-traced clones (e.g. LARRY barcodes) grouped into meta-clones (with
[clomet](https://github.com/loversaber/clomet)) to a differentiation tree and the proliferation and differentiation
rates of each meta-clone on that tree (CLADES, Gao, Barile et al., Nat Commun 2025).

All settings of a dataset are in one file (`config/<dataset>.yaml`); the steps read only that file. Each run writes to
`runs/<run_name>/` with a copy of the config used.

## Steps

| step | what | environment | output in `runs/<run_name>/` |
|---|---|---|---|
| `01_inputs.py` | cell counts per meta-clone x timepoint x population, scaled to absolute numbers; day 0; background row | python | `data/` |
| `02_carta_matrix.py` | clone x fate presence, unique fate profiles | python | `carta/clone_type_mat.csv` |
| `03_carta_sweep.py` | CARTA ILP, tree mode, k = 1..kmax (`--submit`: SLURM CPU job) | carta + Gurobi | `carta/sweep/` |
| `04_carta_tree.py` | discrepancy D(k), elbow, potencies at every k, tree at the chosen k | python | `carta/pareto.*`, `carta/tree_k*` |
| `05a_paga.py` | PAGA connectivity between populations (expression only); PAGA graph on the UMAP | python | `support/paga_connectivity.*`, `support/paga_umap.pdf` |
| `05b_cospar.py` | CoSpar clonal fate coupling, all cells and per timepoint (starts `support.cospar_python` itself) | python + cospar_env | `support/cospar_coupling_*` |
| 05c (to do) | one table of candidate edges with CARTA, CoSpar and PAGA evidence | python | `support/edges.csv` |
| 06 (to do) | transition graph between populations from the curated edge table | python | `data/graph_table.csv` |
| 07 (to do) | CLADES fit, Gillespie, bootstrap (chained SLURM GPU jobs) | CLADES (torch) | `clades/` |
| 08 (to do) | per-meta-clone diagnostics: observed vs predicted, trust label, constraint checks | CLADES (torch, CPU) | `diagnostics/` |
| 09 (to do) | figures | python | `figures/` |

```bash
python steps/01_inputs.py config/hdbi.yaml
python steps/02_carta_matrix.py config/hdbi.yaml
python steps/03_carta_sweep.py config/hdbi.yaml --submit
python steps/04_carta_tree.py config/hdbi.yaml        # after the sweep
python steps/05a_paga.py config/hdbi.yaml
python steps/05b_cospar.py config/hdbi.yaml
```

## Choices in the config

- **day 0** (`day0`): `founders` = number of clones per meta-clone, all in the root population (as in the CLADES cord
  blood example; the background row starts with the number of clones); `d3_total` = the scaled D3 total.
- **CARTA**: cell type → fate map, which clones (`all` or `labelled`), threshold (fate profiles seen in fewer clones are
  dropped), kmax, and the k of the tree (`elbow` or a number).

## Inputs (`inputs/`, not in git)

| file | content |
|---|---|
| `obs_umap.parquet` | one row per cell: clone, cell type, timepoint, experiment, sort labels |
| `metaclone_labels.csv` | one row per clone: `clone_id`, `meta_clone` (from clomet) |
| `pca_harmony.parquet` | harmony PCA (30 components) of the same cells, same order (PAGA, CoSpar) |
| `scaling_factors.pkl` | experiment x timepoint factor from sequenced to absolute cell number |

For the HDBI data these are `inputs/obs_umap.parquet` and `outputs/fig4_metaclones/Fig4_metaclone_labels.csv` of
`hdbi_larry_atlas` and `figure1_analysis_251028/df_scf_all_neg_251116.pkl`.

## External code

- CARTA: github.com/raphael-group/CARTA, commit cd5733c9 (`scripts/larry.py`); needs Gurobi
  (`tools.gurobi_license`).
- CLADES: github.com/StatBiomed/clonaltrans, commit a2d6013.
