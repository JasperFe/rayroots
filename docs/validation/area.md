# vol_area checks

`vol_area` is an experimental extra; `vol_meas` stays the default. Method: [METHOD.md](../METHOD.md#part-b-single-tree).
Script: `tests/area_check.py <single --out folder> <out dir>`, in the conda env. Data not in this repo.

**Synthetic.** Noisy rings (5 mm), one slab: fully covered rings within 4 % of the true area at 4-10 cm. Half-covered rings
(one-sided scan): within 5 % at 10 cm, 8-22 % too small at 4-5 cm with ratio 1.5-2, since rejected slabs fall back to the
circle. A 1 m oval root with ratio 2: the circle alone is 50 % too small, `vol_area` 1-9 %.

**Arch sections** of `SU6_psp4_17` (cm²; outline = convex hull of the points in a 8 cm slab):

| section | outline | rct | meas | area |
|---|---:|---:|---:|---:|
| a | 361 | 104 | 240 | 240 (circle) |
| b | 274 | 62 | 112 | 175 |
| c | 247 | 61 | 121 | 177 |

Section a keeps its circle: two of its three slabs fit half-covered ellipses of 518-533 cm², larger than the whole outline.

**Root totals**, `vol_roots_area` / `vol_roots_meas`: 1.011 (`SU6_psp2_6`), 1.035 (`SU6_psp4_17`), 1.015 (`SU7_psp3_364_136`).

**Slice check.** Horizontal 2 cm slices every 10 cm of the root points; blobs as above, holes filled. Only closed rings
(filled holes add >= 20 %) and only blobs whose cylinders are all steeper than 45 deg. Points = the middle of the point
shell (filled blob minus half its shell) x cos(tilt), against the cylinders crossing the blob:

| tree | rings | rct / points | meas / points | area / points |
|---|---:|---:|---:|---:|
| SU6_psp2_6 | 146 | 0.83 | 0.85 | 0.86 |
| SU6_psp4_17 | 52 | 1.15 | 0.85 | 0.90 |
| SU7_psp3_364_136 | 64 | 0.96 | 0.90 | 0.87 |

`vol_meas` is consistently 10-15 % below the points; `vol_rct` swings between -17 and +15 %. `vol_area` is closer than
`vol_meas` in 95 of 220 blobs. Open arcs (one-sided scans) fill only to the shell and are not used. On three small
hand-cut trees (`SU15_psp1_5`, `SU6_psp2_3`, `SU6_psp4_5`, rebuilt as in [cut.md](cut.md)) RCT models the roots as
many thin cylinders: 3-9 % of the closed-ring area on two of them.

**Rejected: a convex-hull bound** (`area_v2`, an ellipse no larger than the hull of its blob's points, instead of
4 pi r_meas^2). The hull of a half-covered ring is about half the ellipse, so one-sided roots fell back to the circle
(33-50 % too small on half-covered ovals) and section a stayed a circle. Not kept.
