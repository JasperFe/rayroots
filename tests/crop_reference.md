# Expected result on the crop, against two hand-corrected trees
From `tests/compare_reference.py`, ground grid, RCT 492c079. A known result to compare a run against, not a validation:
SU7_psp3_364_136 was used in the method's development. Both trees are cut by the crop. Crop RCT segments differ from
the full plot's, and the scored trees are not whole. Precision near 1 is what this tree gives: it is a clean, isolated
root system (development report: RCT 0.21, v1 0.87); do not expect it elsewhere, see Validation in the README.

## SU7_psp3_364_136, Rhizophora, hand roots (crop segments 22, 56; split stem)
| method | precision | recall | F1 | TP | FP | FN |
|---|---:|---:|---:|---:|---:|---:|
| RCT | 1.000 | 0.120 | 0.214 | 61,131 | 14 | 449,507 |
| v1 | 0.999 | 0.762 | 0.865 | 389,055 | 242 | 121,583 |
| v2 | 1.000 | 0.423 | 0.594 | 215,851 | 52 | 294,787 |
| union | 0.999 | 0.795 | 0.886 | 405,956 | 275 | 104,682 |

## SU7_psp3_119, Avicennia (crop segment 8)
Hand-corrected tree, no hand roots: an overlap check only, not a score of root finding on Avicennia.

| method | precision | recall | F1 | TP | FP | FN |
|---|---:|---:|---:|---:|---:|---:|
| RCT | 0.994 | 0.667 | 0.798 | 29,624 | 171 | 14,791 |
| v1 | 0.939 | 0.925 | 0.932 | 41,090 | 2,670 | 3,325 |
| v2 | 0.989 | 0.958 | 0.973 | 42,555 | 467 | 1,860 |
| union | 0.941 | 0.978 | 0.959 | 43,456 | 2,721 | 959 |
