# Input files

A dataset needs four files (paths in the `inputs` and `support` sections of its config). `examples/cord_blood_clades/inputs/`
has a complete set, made by `examples/cord_blood_clades/prepare_input.py`. Files of unpublished data go in this folder,
which is not in git.

| file (config key) | one row per | columns |
|---|---|---|
| cells (`inputs.cells`), parquet | cell (index: cell name) | the columns named in `columns`: clone (empty if the cell has no barcode), cell type, timepoint (names as in `timepoints`), experiment; the UMAP coordinates (`columns.umap`); any column used in `cell_filter` |
| meta-clone labels (`inputs.labels`), csv | clone | `columns.label_clone` (the clone names used in the cell table), `columns.label_group` (meta-clone) |
| scaling factors (`inputs.scaling`), csv or pkl | experiment (index: as in the experiment column) | one column per timepoint (names as in `timepoints`): factor from a sequenced cell to the absolute number of cells it stands for |
| PCA (`support.pca`), parquet | cell, same cells in the same order as the cell table | the components (any names) |

Each cell carries at most one clone. Only experiments in the scaling table are used for the CLADES counts.

Example (cord blood): `cells.parquet` has `clone`, `celltype`, `timepoint` (`Day3`, `Day10`, `Day17`), `experiment`,
`GFP`, `UMAP1`, `UMAP2`; `metaclone_labels.csv` has `clone_id`, `meta_clone`; `scaling_factors.csv`:

```
experiment,Day3,Day10,Day17
CB_HDBI58,4.905,530.598,508.657
```
