# clokin tree stage vs the CLADES cord blood graph

The CLADES fit in this example uses the published CLADES graph unchanged (`graph` in `config.yaml`, from
`demo/CordBlood_Refine/graph_table.csv`). The tree stage of clokin was run on the same data to test it and to compare
its evidence with that graph. Files: `support/edges.csv`, `support/potency.csv`, `support/trees.pdf`,
`carta/tree_k3.pdf` (in this folder).

## CARTA

Discrepancy D(k) for k = 1-4: 480, 218, 128, 128 (262 clones in 10 multi-fate profiles enter the ILP). D(k) stops
falling at k = 3 (the elbow is k = 2). Progenitors at k = 3: {Ery, Mast} (131 clones), {Mono, DC} (45 clones) and
{Mono, Ery, Mast, DC} (86 clones).

- {Ery, Mast} is MEMP -> {Mast cell, Early Erythroid} in the CLADES graph.
- {Mono, DC} has no single node in the CLADES graph, where Mono precursor and DC precursor come from HSC/MPP 1 and 2
  separately. The clones of HSC/MPP 1 make both (D3 clones seen again: Monocyte 81 %, DC precursor 41 %, DC 31 %),
  those of HSC/MPP 2 almost only monocytes (91 %; DC 5 %), so HSC/MPP 1 is closest to the {Mono, DC} progenitor.
- CoSpar (clonal coupling) groups the fates the same way: Ery-Mast 0.044 is the strongest pair, then Mono-DC 0.020.
  PAGA (expression) puts Mono and DC closest.

## Edges of the CLADES graph (13)

| evidence | edges |
|---|---|
| clonal and PAGA agree (7) | HSC/MPP 1 -> HSC/MPP 2, HSC/MPP 1 -> DC precursor, HSC/MPP 2 -> NMP, HSC/MPP 2 -> Mono precursor, MEMP -> Mast cell, MEMP -> Early Erythroid, Early -> Mid Erythroid |
| clonal only (1) | HSC/MPP 1 -> MEMP (6 % of the clones; PAGA 0.09) |
| PAGA only (2) | HSC/MPP 2 -> MEMP, HSC/MPP 2 -> DC precursor (each 4.5 % of 22 clones) |
| no clonal data (3) | Mid -> Late Erythroid, Mono precursor -> Monocyte, DC precursor -> DC: their D3 clones are not seen again later |

No edge of the CLADES graph is contradicted. Clonal evidence is thin in this data set: of the D3 clones of each
progenitor, only some are seen again at day 10 or 17 (HSC/MPP 1: 32 of 177, HSC/MPP 2: 22 of 107, NMP: 46 of 169,
MEMP: 98 of 559, Mono precursor 0 of 45, DC precursor 1 of 48).

## Candidates not in the CLADES graph

| edge | evidence | note |
|---|---|---|
| HSC/MPP 1 -> DC | clonal 31 %, PAGA 0.50 | probably the route HSC/MPP 1 -> DC precursor -> DC that the graph already has |
| NMP -> Mono precursor | PAGA 0.19, clonal 2 % | not supported by the clones |
| Mono precursor -> DC precursor | PAGA 0.48, no clonal data | expression similarity only |
| Late Erythroid - DC | PAGA 0.15 | two end populations: not a differentiation step |

NMP is an end point in the CLADES graph, but 91 % of the NMP clones seen again later (42 of 46) have monocytes; the PAGA
link NMP - Monocyte is 0.148, just below the candidate cutoff (0.15). The clones suggest that NMP also feeds the
monocytes.

## Used for the CLADES fit

The CLADES graph, unchanged, so that the example reproduces the CLADES analysis.

## CLADES v1.2.0 graph

The graph of CLADES v1.2.0 (`inputs/clades_graph_table_v1.2.csv`, used by `config_v1.2.yaml`) has the 13 edges above
plus MEMP -> Mono precursor and NMP -> Mono precursor.

- NMP -> Mono precursor: PAGA 0.19 (a candidate edge here), and 91 % of the NMP clones seen again have monocytes
  (only 2 % have Mono precursor cells, which are few after day 3). It fits the NMP -> monocyte route the clones show.
- MEMP -> Mono precursor: no support in this data (PAGA 0.009; not a candidate edge).
