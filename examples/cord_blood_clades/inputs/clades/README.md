# CLADES input of the example

Made by clokin from the files one folder up (`python run.py tree` / `kinetics`, steps 1 and 6), the same files that go
into `runs/cord_blood/data/`. They are what CLADES reads (`data_loader` in the CLADES config):

| file | content |
|---|---|
| `kinetics.txt` | 13 rows (Clone 0 .. Clone 11, BG) x 36 columns (Day3, Day10, Day17 x 12 populations): absolute cell numbers |
| `initial_condition.csv` | day-0 cells per row, all in the first population (HSC/MPP 1): the number of clones; BG: all clones |
| `annotations.csv` | row names and population names |
| `graph_table.csv` | transition graph (1 = row population -> column population) |
| `clades_params.json` | `t_observed`, `scaling_factor`, `no_proliferation_pops`, `no_apoptosis_pops` for the CLADES config |
| `counts_long.csv` | the counts of `kinetics.txt` as a table, with the sequenced (raw) cell numbers |

Compared with `demo/CordBlood_Refine` of the CLADES code: `kinetics.txt` = `kinetics_array_correction_factor.txt`,
`initial_condition.csv` and `graph_table.csv` identical, the same populations, and the same parameters as
`config/main_dynamic_cordblood.json`. Only the background row is named `BG` here (`Clone BG` there).
