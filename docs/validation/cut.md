# The cut line on hand-cut trees

The automatic cut against 28 Rhizophora trees cut by hand. Only the cut is tested, not the multi-stem split, the volumes or the DBH.

| trees | 28 |
|---|---:|
| cut found | 27 |
| mean difference | -0.05 m |
| mean absolute difference | 0.38 m |
| within 0.5 m | 22/27 |
| flagged `short run` | 3 |
| no cut found | 1 |

- **Worst cases:** `SU14_psp3_4` is 1.71 m too high (flagged `short run`). `SU15_psp1_6` is 1.34 m too low and not flagged: two stems
  were detected, and 26 of the 27 trees with a cut got one stem. `SU14_psp3_13` got no cut.
- **Script:** `validate_cut.py` in this folder, with the final code. Run it in the conda env: `python validate_cut.py <data folder> <out dir>`.
- **Data, not in this repo:** hand-cut trees of the Suriname TLS study, Q-ForestLab, UGent. Per tree, a top and a roots file cut by hand,
  in `tops/` and `roots/`, and the height of each cut in `roots/z_refdb_part2.csv`. Each tree is put back together from the two
  parts, then the cut is searched for. The clouds hold 36 k to 1.5 M points.

`h_hand` and `h_auto` are heights above `z0`, the 0.1th percentile of the tree's z. `diff_m` is automatic minus hand.

| tree | h_hand | h_auto | diff_m | n_stems | clean_run | flag |
|---|---|---|---|---|---|---|
| Suriname_SU14_psp3_12_a | 2.3 | 2.54 | 0.24 | 1 | 0.66 |  |
| Suriname_SU14_psp3_13_a | 1.91 |  |  | 0 | 0.0 | no cut found |
| Suriname_SU14_psp3_4_a | 3.47 | 5.18 | 1.71 | 1 | 0.4 | short run |
| Suriname_SU14_psp3_7_a | 2.1 | 2.48 | 0.39 | 1 | 0.72 |  |
| Suriname_SU15_psp1_11_a | 1.1 | 1.01 | -0.09 | 1 | 0.42 | short run |
| Suriname_SU15_psp1_12_a | 1.78 | 1.82 | 0.04 | 1 | 1.74 |  |
| Suriname_SU15_psp1_13_a | 1.14 | 1.36 | 0.22 | 1 | 1.8 |  |
| Suriname_SU15_psp1_1_a | 2.97 | 3.22 | 0.25 | 1 | 7.36 |  |
| Suriname_SU15_psp1_2_a | 2.99 | 2.28 | -0.7 | 1 | 8.3 |  |
| Suriname_SU15_psp1_3_a | 2.88 | 2.42 | -0.45 | 1 | 5.9 |  |
| Suriname_SU15_psp1_4_a | 3.14 | 2.74 | -0.39 | 1 | 3.84 |  |
| Suriname_SU15_psp1_5_a | 2.2 | 2.46 | 0.26 | 1 | 8.12 |  |
| Suriname_SU15_psp1_6_a | 4.37 | 3.03 | -1.34 | 2 | 0.62 |  |
| Suriname_SU15_psp1_7_a | 2.66 | 2.98 | 0.32 | 1 | 1.58 |  |
| Suriname_SU15_psp1_8_a | 2.58 | 2.81 | 0.24 | 1 | 0.56 | short run |
| Suriname_SU15_psp1_9_a | 2.54 | 2.7 | 0.16 | 1 | 3.02 |  |
| Suriname_SU6_psp2_1_a | 2.77 | 2.6 | -0.17 | 1 | 7.98 |  |
| Suriname_SU6_psp2_3_a | 3.43 | 3.14 | -0.29 | 1 | 7.44 |  |
| Suriname_SU6_psp2_4_a | 4.63 | 4.26 | -0.37 | 1 | 3.7 |  |
| Suriname_SU6_psp2_5_a | 3.77 | 3.08 | -0.69 | 1 | 7.5 |  |
| Suriname_SU6_psp2_6_a | 4.64 | 4.5 | -0.13 | 1 | 6.08 |  |
| Suriname_SU6_psp2_7_a | 4.57 | 4.06 | -0.51 | 1 | 6.52 |  |
| Suriname_SU6_psp4_2_a | 3.29 | 3.02 | -0.27 | 1 | 4.0 |  |
| Suriname_SU6_psp4_3_a | 2.06 | 2.26 | 0.2 | 1 | 6.14 |  |
| Suriname_SU6_psp4_4_a | 4.02 | 4.14 | 0.12 | 1 | 2.34 |  |
| Suriname_SU6_psp4_5_a | 2.68 | 2.4 | -0.28 | 1 | 8.18 |  |
| Suriname_SU6_psp4_6_a | 3.26 | 3.44 | 0.19 | 1 | 7.14 |  |
| Suriname_SU6_psp4_7_a | 4.0 | 3.88 | -0.12 | 1 | 6.7 |  |
