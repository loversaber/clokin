# clokin

**CLOnes to KINetics.** From lineage-traced clones (e.g. LARRY barcodes) grouped into meta-clones (with
[clomet](https://github.com/loversaber/clomet)) to a differentiation tree and the proliferation and differentiation
rates of each meta-clone on that tree (CLADES, Gao, Barile et al., Nat Commun 2025).

All settings of a dataset are in one file (`config/<dataset>.yaml`); the steps read only that file. Each run writes to
`runs/<run_name>/` with a copy of the (merged) config used.

## Quick start

Everything a dataset needs is in one config file; `python` below is the environment in `tools.python`.

```bash
python run.py tree config/hdbi.yaml       # 1. steps 1-5 as one SLURM CPU job (about 1 h)
                                          # 2. decide: look at runs/<run_name>/carta and support (see below), write
                                          #    the chosen edges in the `graph` section of the config
python run.py kinetics config/hdbi.yaml   # 3. steps 6-7: graph, CLADES fit -> Gillespie (GPU); steps 8-9 run after each
python run.py status config/hdbi.yaml     #    jobs, STATUS.txt and which outputs exist
```

To check that everything runs (for example after changing the code or the environments), use the quick test: same
data, small CARTA sweep, 20 epochs, 2 simulations per row (`config/hdbi_quicktest.yaml`, run `quicktest`).

### Decision between the stages

| look at | to decide |
|---|---|
| `carta/pareto.pdf`, `carta/potencies_by_k.csv` | the k of the tree (`carta.k`): the elbow is marked; potencies that stay as k grows are robust |
| `carta/tree_k<k>.pdf` | the progenitors (sets of fates) and how they nest |
| `support/trees.pdf` | whether CARTA, CoSpar (both clonal) and PAGA (expression) group the fates the same way |
| `support/paga_umap.pdf` | the PAGA graph on the UMAP |
| `support/edges.csv`, `support/potency.csv` | each candidate edge with its evidence and status (`agree` / `paga only` / `no clonal data` / ...) |

Edges with status `agree` are usually kept; `no clonal data` edges (source too rare at the first timepoint) and
`paga only` edges need a judgement; `mature-mature` edges are not differentiation steps. Write the edges you keep in
`graph` (source: [targets]), with a comment on what changed and why. Step 6 checks every population is reachable from
the root.

### Results to read (`runs/<run_name>/`)

| folder | content |
|---|---|
| `diagnostics/fit_by_row.csv` | per row: observed vs predicted cells per timepoint, error, trust label (`use` / `do not quote`) |
| `diagnostics/checks.csv`, `convergence.csv` | rates that break the assumptions; loss over the epochs |
| `diagnostics/rates.csv` | fitted differentiation and net growth rates per row and day |
| `diagnostics/fate_prob.csv`, `divisions.csv` | Gillespie: chance to reach each population from one root cell; divisions to get there |
| `figures/` | the same as PDFs, plus the graph on the UMAP and the trajectories |

Notes on reading the results of a dataset (caveats specific to its data) can go in `notes/`, which is not in git.

## A new dataset

Copy `config/hdbi.yaml` and change: `inputs` (cell table, meta-clone labels, scaling factors, PCA) and `columns`;
`cell_filter`; `timepoints`; `populations` (cell type -> population, root first) with `terminal`,
`no_proliferation` and `early`; `carta.fates`; `colors`; `tools` (environments, SLURM accounts); `run_name`. Start
`graph` empty, run the tree stage, then fill it in. A variant of an existing config can start with
`base: <config>.yaml` and list only what differs (as `hdbi_quicktest.yaml` does).

## Steps

| step | what | environment | output in `runs/<run_name>/` |
|---|---|---|---|
| `01_inputs.py` | cell counts per meta-clone x timepoint x population, scaled to absolute numbers; day 0; background row | python | `data/` |
| `02_carta_matrix.py` | clone x fate presence, unique fate profiles | python | `carta/clone_type_mat.csv` |
| `03_carta_sweep.py` | CARTA ILP, tree mode, k = 1..kmax (`--submit`: on its own as a SLURM CPU job) | carta + Gurobi | `carta/sweep/` |
| `04_carta_tree.py` | discrepancy D(k), elbow, potencies at every k, tree at the chosen k | python | `carta/pareto.*`, `carta/tree_k*` |
| `05a_paga.py` | PAGA connectivity between populations (expression only); PAGA graph on the UMAP | python | `support/paga_connectivity.*`, `support/paga_umap.pdf` |
| `05b_cospar.py` | CoSpar clonal fate coupling, all cells and per timepoint (starts `support.cospar_python` itself) | python + cospar_env | `support/cospar_coupling_*` |
| `05c_edges.py` | candidate edges (previous graph + PAGA links) with clonal share, PAGA and CoSpar evidence and a status | python | `support/edges.csv`, `support/potency.csv` |
| `05d_trees.py` | fate similarity and fate trees from CARTA (progenitor co-membership), CoSpar (clonal coupling) and PAGA (graph distance), side by side | python | `support/trees.pdf`, `support/fate_similarity_*.csv` |
| `06_graph.py` | transition graph from the `graph` section of the config; checks every population is reachable | python | `data/graph_table.csv`, `support/edges_used.csv` |
| `07_clades.py` | CLADES configs and job scripts; `--submit`: fit, then Gillespie after the fit; `--post`: steps 8-9 after each (SLURM) | python | `clades/` |
| `08_diagnostics.py` | observed vs predicted per row, trust label, rates, constraint checks, convergence; Gillespie fate probabilities and divisions (starts `clades.python` itself, CPU) | python + CLADES env | `diagnostics/` |
| `09_figures.py` | graph on the UMAP, fit, trajectories, composition, rates, fate probabilities, divisions, convergence | python | `figures/` |

Each step can also be run on its own: `python steps/<step>.py config/<dataset>.yaml`. Shared code (config loading
with `base:`, cell selection) is in `steps/common.py`.

## Choices in the config

- **day 0** (`day0`): `founders` = number of clones per meta-clone, all in the root population (as in the CLADES cord
  blood example; the background row starts with the number of clones); `d3_total` = the scaled D3 total.
- **no_proliferation**: populations that may not divide in CLADES (soft penalty). As in the CLADES examples, the
  mature end populations only; precursors such as Pre-cDC (`DC`) are free.
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

## Methods used (please cite)

clokin organises the steps; the methods are those of the original authors.

| step | method | reference |
|---|---|---|
| meta-clones (input) | clomet | github.com/loversaber/clomet |
| 2–4 | CARTA, cell differentiation maps from lineage tracing | Sashittal et al., Raphael lab (full reference: TODO) |
| 5a | PAGA | Wolf et al., Genome Biology 2019 |
| 5b | CoSpar | Wang et al., Nature Biotechnology 2022 |
| 7–8 | CLADES | Gao, Barile et al., Nature Communications 2025 (full reference: TODO) |

## External code

clokin does not contain these programs; it runs them from a local copy. The commit is the version used, so that
results can be reproduced if the projects change later.

- CARTA: github.com/raphael-group/CARTA, commit cd5733c9. Step 3 runs `scripts/larry.py`, an integer linear
  program solved with Gurobi (licence: `tools.gurobi_license`).
- CLADES: github.com/StatBiomed/clonaltrans (the CLADES code), commit a2d6013. Steps 7 and 8 run `main.py`
  (fit) and `main_gillespie.py` (Gillespie simulations).
- CoSpar 0.5.0 and scanpy (PAGA), installed in the environments named in the config.
