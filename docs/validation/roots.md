# Root proposals on unseen plots

**Where these numbers come from.** The untouched test of the Suriname TLS study: the frozen v1 and v2 were run once on plots
not used in their development, and nothing was tuned on these trees.

- **Plots:** SU6_psp2 (7 trees) and SU14_psp3 (1 tree).
- **Reference:** hand-cleaned trees with their roots.
- **Domain:** within 5 m of the stem, below 4 m. Counts are points. RCT = the RCT segment without root proposals.

These scores measure the **root proposals, which are meant to be corrected by hand**. They are not the accuracy of a finished
tree. The scoring script and the reference trees are not in this repo.

## Per colour
The share of the points in each colour of the review files that is hand-cleaned wood, same domain.

| colour | found by | share that is hand-cleaned wood, SU6_psp2 | SU14_psp3 | meaning |
|---|---|---:|---:|---|
| green | both | 0.92 | 0.98 | mostly keep |
| blue | v2 only | 0.77 | 0.76 | mostly keep |
| orange | v1 only | 0.35 | 0.81 | check |
| grey | RCT | 0.33 | 0.61 | |

SU6_psp2 is pooled over 6 trees, without SU6_psp2_2 (zero by construction); SU14_psp3 is 1 tree, so its orange value carries little weight.

## SU6_psp2, 7 hand-root trees, pooled
| method | band | precision | recall | F1 | TP | FP | FN |
|---|---|---:|---:|---:|---:|---:|---:|
| RCT | all | 0.314 | 0.331 | 0.322 | 518,650 | 1,131,110 | 1,049,214 |
| v1 | all | 0.394 | 0.717 | 0.508 | 1,124,279 | 1,730,015 | 443,585 |
| v2 | all | 0.437 | 0.589 | 0.502 | 923,061 | 1,186,927 | 644,803 |
| RCT | below 2 m | 0.480 | 0.236 | 0.317 | 299,530 | 324,478 | 968,322 |
| v1 | below 2 m | 0.526 | 0.691 | 0.597 | 876,190 | 789,235 | 391,662 |
| v2 | below 2 m | 0.645 | 0.532 | 0.583 | 675,002 | 371,867 | 592,850 |

## SU14_psp3, 1 hand-root tree (SU14_psp3_4)
| method | band | precision | recall | F1 | TP | FP | FN |
|---|---|---:|---:|---:|---:|---:|---:|
| RCT | all | 0.612 | 0.225 | 0.329 | 29,096 | 18,466 | 100,460 |
| v1 | all | 0.824 | 0.941 | 0.878 | 121,862 | 26,096 | 7,694 |
| v2 | all | 0.823 | 0.768 | 0.795 | 99,544 | 21,405 | 30,012 |

Notes:
- SU6_psp2_2 scores zero by construction: _2 and _3 are scored on the same predicted points (shared segments), and all
  the matching hand wood belongs to _3. Pooled without _2: F1 0.534 for v1, 0.522 for v2.
- Precision above 2 m is about 0.2 in every method, RCT included (0.21-0.25); most of those false positives are already
  in the RCT segments. The below-2 m rows are the fairer view of the root step.

**Not validated:**
- The union of v1 and v2, and the v2-wins rule for contested points: v1 and v2 were scored separately.
- Avicennia roots: the test trees are Rhizophora.
- `ground = mesh`.
- Any other site, scanner or scan setting.
- The whole-cloud mode (no boundary), where edge trees have one-sided data.
- Only 8 hand-root trees in 2 plots carry these numbers.
