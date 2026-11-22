# CLADES v1.2.0 fit of the example (`config.yaml`)

Run with CLADES v1.2.0 (github.com/StatBiomed/clonaltrans, release v1.2.0, commit a0442270), the published version.

| folder / file | content |
|---|---|
| `input/` | the CLADES input made by step 1 and step 6: `kinetics.txt` (11 rows: Clone 0-8, Clone 10, BG; day 0 block first, then Day3, Day10, Day17 x 12 populations), `day0.csv`, `annotations.csv`, `graph_table.csv` (15 edges), `edges_used.csv`, `clades_params.json`, `counts_long.csv` |
| `fit.json`, `gillespie.json` | the CLADES configs clokin wrote |
| `diagnostics/` | step 8: `fit_by_row.csv` (observed vs predicted cells per row and day, trust label), `checks.csv`, `convergence.csv`, `composition_last.csv`, `rates.csv`, `trajectories.csv`; Gillespie: `fate_prob.csv`, `divisions.csv`, `gillespie_starts.csv` |
| `figures/` | step 9 |

## Input

`kinetics.txt` is identical to `demo/Human_CordBlood/kinetics_corrected.txt` of CLADES v1.2.0 (all 528 numbers) and
`graph_table.csv` to its `paga_transitions.csv`. Day 0 as in v1.2.0: the founders of each meta-clone are split over
HSC/MPP 1, HSC/MPP 2, MEMP, NMP and Mono precursor in proportion to the day-3 cells of the clones seen at all three
timepoints. Meta-clones 9 and 11 have no such clone, get no day-0 cells and are left out (as in v1.2.0).

## Fit

1500 epochs. The loss levels off after about 500 epochs (change over the last 20 %: -0.1 %).

- Day 10 and day 17 are fitted well: log2(predicted / observed) is within 0.07 at day 17 and within 0.62 at day 10 for the
  7 meta-clones labelled "use" (more than 3 populations at day 17) and for BG.
- Day 3 is over-predicted in the larger meta-clones (up to 14x in Clone 7, 40 observed cells); its counts are small
  compared with day 10 and 17, so it weighs little in the loss. The old CLADES version shows the same.
- Clone 5, 8 and 10 have at most 3 populations at day 17: labelled "do not quote".
- `checks.csv`: 7 flags, all Monocyte (and Mast cell of BG at day 3) growing although set as no-proliferation; the
  setting is a penalty in CLADES, not a hard constraint.
- In 592 of the 1500 epochs the loss is NaN (one or two of the 132 observed values). This is CLADES v1.2.0 behaviour:
  every 100 epochs it goes back to the best model when the loss is NaN or higher. The model saved at the end
  (`model_last.pt`, used for Gillespie) is the one of the last epoch; all its predictions are finite.

## Gillespie

1,000 simulations per row, each from one founder cell that starts in a population drawn from the row's day 0
(`gillespie_starts.csv`), run to day 17 (`t_cutoff` 1e-3, as the v1.2.0 config).

- `fate_prob.csv`: share of the simulations that reach each population. Of the rows labelled "use", Clone 2 and
  Clone 7 reach all five end populations (Mast cell, Late Erythroid, NMP, Monocyte, DC) in more than half of the
  simulations; Clone 1 and Clone 3 reach Mast cell, Late Erythroid and Monocyte but not NMP or DC; Clone 4 mostly
  NMP and Monocyte (about 0.9 each); Clone 6 every end population except Late Erythroid.
- `divisions.csv`: mean number of divisions from the founder to the first cell of each population (0 for the
  population it starts in).
