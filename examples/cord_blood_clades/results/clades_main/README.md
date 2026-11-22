# Older CLADES code (main branch) fit of the example (`config_main.yaml`)

Run with the CLADES code of the main branch (`demo/CordBlood_Refine`, as used before release v1.2.0). For the
published release see `../clades_v1.2/`.

| folder / file | content |
|---|---|
| `input/` | the CLADES input made by step 1 and step 6: `kinetics.txt` (13 rows: Clone 0-11, BG; Day3, Day10, Day17 x 12 populations), `initial_condition.csv` (day 0), `annotations.csv`, `graph_table.csv` (13 edges), `edges_used.csv`, `clades_params.json`, `counts_long.csv` |
| `fit.json`, `gillespie.json` | the CLADES configs clokin wrote |
| `diagnostics/` | step 8, the same tables as `../clades_v1.2/diagnostics/` |
| `figures/` | step 9 |

## Input

Identical to `demo/CordBlood_Refine` (all 468 cell counts, `initial_condition.csv`, scaling factors). Day 0: all
founders of each meta-clone in HSC/MPP 1.

## Fit

1500 epochs, MultiStep learning rate. The loss falls from 230 to 75 by epoch 1200 and changes by -2.4 % over the last
20 % of the epochs.

- Day 17 is fitted well: |log2(predicted / observed)| <= 0.09 for the 7 meta-clones labelled "use" (Clone 0-4, 6, 7)
  and BG; day 10 within 0.91 (Clone 2 under-predicted); day 3 over-predicted (up to 3.4 log2), as in v1.2.0.
- Clone 5, 8, 9, 10 and 11 have at most 3 populations at day 17: labelled "do not quote".
- `checks.csv`: 11 flags, Mast cell and Monocyte growing although set as no-proliferation (a penalty in CLADES, not a
  hard constraint).

## Gillespie

1,000 simulations per row from one founder in HSC/MPP 1; 855-996 per row reach day 17 with cells (`simulations`
column of `fate_prob.csv`).

- Clone 7 and BG reach every population in more than 0.55 of the simulations; Clone 0 and Clone 2 reach Mast cell,
  Late Erythroid, NMP and Monocyte in about half or more; Clone 3 never reaches NMP or DC; Clone 6 never reaches the
  erythroid populations.
- Compared with v1.2.0 (rows labelled "use" in both): the same rows never reach DC (Clone 0, 1, 3) or Late
  Erythroid (Clone 6), but the shares differ by up to 0.6 (Clone 4 Mast cell 0.80 here vs 0.28 in v1.2.0; Clone 3
  Monocyte 0.19 vs 0.81), as the two versions start the founders in different populations and fit differently.
