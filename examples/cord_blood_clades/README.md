# Example: CLADES cord blood LARRY data

Human cord blood LARRY lineage tracing (days 3, 10 and 17; 24,885 cells, 3,246 clones), published with

> Gao M, Barile M, Chabra S, Haltalli M, Calderbank EF, Chao Y, Zheng W, Wilson NK, Laurenti E, Göttgens B, Huang Y.
> CLADES: a hybrid NeuralODE-Gillespie approach for unveiling clonal cell fate and differentiation dynamics.
> Nat Commun 16 (2025). https://doi.org/10.1038/s41467-025-63150-6

Raw data: GEO GSE276896. Processed data: figshare https://doi.org/10.6084/m9.figshare.27908142 (`CordBlood_Refine`),
licence CC BY 4.0. The files in `inputs/` are derived from the processed data. The CLADES graph is from
`demo/CordBlood_Refine/graph_table.csv` of the CLADES code (github.com/StatBiomed/clonaltrans, MIT licence).

The example uses the 12 CLADES meta-clones and the settings of the CLADES cord blood run, so that it reproduces the
CLADES analysis.

## Run

From the clokin folder, after setting up `config/site.yaml` (see the main README):

```bash
python run.py tree examples/cord_blood_clades/config.yaml       # steps 1-5 (minutes; add --here to run without SLURM)
python run.py kinetics examples/cord_blood_clades/config.yaml   # CLADES fit -> Gillespie (GPU), steps 8-9 after each
python run.py status examples/cord_blood_clades/config.yaml
```

Outputs go to `runs/cord_blood/` (not in git). `quicktest.yaml` runs the same in minutes (20 epochs, 2 simulations per
row) to check an installation.

## Files

| file | content |
|---|---|
| `prepare_input.py` | makes `inputs/` from the figshare AnnData, `CordBlood_Obs_Clones.pkl` and the CLADES code |
| `inputs/cells.parquet` | one row per cell: clone, cell type (`def_lab`), timepoint, experiment, GFP, UMAP |
| `inputs/pca.parquet` | PCA of the same cells (50 components) |
| `inputs/metaclone_labels.csv` | CLADES meta-clone (0-11) of the 606 clones with at least 3 GFP+ cells |
| `inputs/scaling_factors.csv` | sequenced -> absolute cells per timepoint (factors of the CLADES notebook) |
| `inputs/clades_graph_table.csv` | the CLADES transition graph |
| `inputs/clades/` | the CLADES input clokin makes from these files (identical to `demo/CordBlood_Refine`); see its README |
| `config.yaml`, `quicktest.yaml` | settings of the example and of the quick test |
| `results/` | small outputs of the example run: CARTA (`carta/`), PAGA / CoSpar / candidate edges (`support/`), and `graph_comparison.md` |

## Checks

- Step 1 rebuilds the CLADES cord blood input (`demo/CordBlood_Refine` in the CLADES code) exactly: all 468 cell
  counts of `kinetics_array_correction_factor.txt`, `initial_condition.csv`, the scaling factors and the
  proliferation / apoptosis settings (`clades_rounding: true` rounds as the CLADES notebook).
- The tree stage (CARTA, PAGA, CoSpar) supports or does not contradict every edge of the CLADES graph; the CLADES
  fit uses that graph unchanged. Details: `results/graph_comparison.md`.
