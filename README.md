# clokin

**CLOnes to KINetics.** From lineage-traced clones (e.g. LARRY barcodes) grouped into meta-clones (with
[clomet](https://github.com/loversaber/clomet)) to a differentiation tree and the proliferation and differentiation
rates of each meta-clone on that tree (CLADES, Gao, Barile et al., Nat Commun 2025).

All settings of a dataset are in one file (`config/<dataset>.yaml`); the steps read only that file. Each run writes to
`runs/<run_name>/` with a copy of the (merged) config used.

## Setup

1. Three Python environments: one for clokin itself (pandas, numpy, scipy, matplotlib, scanpy, anndata, pyarrow,
   pyyaml), one for CARTA (its requirements and gurobipy, with a Gurobi licence, free for academic use), one for
   CLADES (its requirements: PyTorch with CUDA, torchdiffeq, natsort), plus CoSpar 0.5.0 in the first or a fourth one.
2. Copy `config/site_example.yaml` to `config/site.yaml` (not in git) and fill in the environments, the CARTA and
   CLADES folders, the Gurobi licence and the SLURM accounts and partitions. Dataset configs start with
   `base: <path>/config/site.yaml`.

## Example: CLADES cord blood data

`examples/cord_blood_clades/`: the human cord blood LARRY data published with CLADES (Gao, Barile et al., Nat Commun
2025; figshare 10.6084/m9.figshare.27908142, CC BY 4.0), with the 12 CLADES meta-clones and the CLADES settings.
`inputs/` (made by `prepare_input.py` from the figshare files) and `results/` are included.

```bash
python run.py tree examples/cord_blood_clades/config.yaml       # steps 1-5
python run.py kinetics examples/cord_blood_clades/config.yaml   # graph (the CLADES one, in the config), CLADES, steps 8-9
python run.py status examples/cord_blood_clades/config.yaml
```

Step 1 rebuilds the CLADES cord blood input (`demo/CordBlood_Refine` in the CLADES code) exactly: all 468 cell counts,
the day-0 numbers and the CLADES parameters (`clades_rounding: true` rounds as the CLADES notebook does).
`examples/cord_blood_clades/quicktest.yaml` runs the whole pipeline on the example in minutes, to check an installation.

## Quick start (your data)

```bash
python run.py tree config/<dataset>.yaml       # 1. steps 1-5 as one SLURM CPU job
                                               # 2. decide: look at runs/<run_name>/carta and support (below), write
                                               #    the chosen edges in the `graph` section of the config
python run.py kinetics config/<dataset>.yaml   # 3. steps 6-7: graph, CLADES fit -> Gillespie (GPU); 8-9 after each
python run.py status config/<dataset>.yaml     #    jobs, STATUS.txt and which outputs exist
```

`python` is the environment in `tools.python`. Configs and inputs of unpublished data stay out of git (`config/`
keeps only `site_example.yaml` in git, `inputs/` only its README).

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

The Gillespie simulations (`clades/gillespie/`) take several GB per run (22 GB for 20 rows x 1,000 simulations).

Notes on reading the results of a dataset (caveats specific to its data) can go in `notes/`, which is not in git.

## A new dataset

Make the four input files described in `inputs/README.md`, copy `examples/cord_blood_clades/config.yaml` to
`config/<dataset>.yaml` and change: `base` and `root` (or remove `root` for a config in `config/`), `run_name`,
`inputs` and `columns`, `cell_filter`, `timepoints`, `populations` (cell type -> population, root first), `terminal`,
`no_proliferation` and `early`, `day0`, `colors`, `carta.fates`, `edges` (`previous_graph` is optional) and
`diagnostics`. Leave `graph` empty, run the tree stage, then fill it in. A variant of a config can start with
`base: <config>.yaml` and list only what differs (as `quicktest.yaml` does).

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

- **CLADES version** (`clades.version`): `'1.2'` = the published release (branch v1.2.0: day 0 inside the counts,
  adaptive learning rate, rate clipping, size-weighted background penalty, Gillespie lineages started as in day 0) or
  `main` = the code before it (commit a2d6013: separate day-0 file, fixed learning-rate milestones).
- **day 0** (`day0`): one founder cell per clone (the background row: all clones). `founders` puts them all in the
  root population (CLADES before v1.2.0); `founders_split` splits them over `day0_populations` as the row's
  first-timepoint cells in those populations (CLADES v1.2.0; `day0_split_clones: seen_at_all_timepoints` uses only
  clones seen at every timepoint, as v1.2.0 does); `d3_total` = the scaled D3 total.
- **keep_groups**: the meta-clones to model (default all); the background row always uses all cells.
- **no_proliferation**: populations that may not divide in CLADES (soft penalty). As in the CLADES examples, the
  mature end populations only; precursors such as Pre-cDC (`DC`) are free.
- **CARTA**: cell type → fate map, which clones (`all` or `labelled`), threshold (fate profiles seen in fewer clones are
  dropped), kmax, and the k of the tree (`elbow` or a number).

## Inputs

Formats: `inputs/README.md`. Example: `examples/cord_blood_clades/inputs/`.

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
