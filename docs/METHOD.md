# rayroots: details

Back to the [README](../README.md).

## Inputs
- A plot LAS (`<stem>.las`), or a `<stem>_segmented.ply` from `rayextract trees`.
- Optional boundary, csv with columns `x,y`: scan positions or a polygon, **in the segmented frame**. The plot is
  the convex hull of these points (a concave polygon gets filled in). Without a boundary the whole cloud is used,
  and trees at the edge of the scan get proposals from one-sided data.

**Segmented frame** = LAS coordinates minus the LAS header offset (translation only).
`rayimport --remove_start_pos` subtracts the header offset, not the first point as its help says. Take the offset
from `las_offset_exact` in `<stem>.rct.log`: the value rayimport prints has 6 significant digits, which is up to
0.5 m off. rayroots stops when the boundary misses the cloud entirely, but cannot detect a small shift.

## Part A: plot
All steps write in `<project>/A_plot/`. A step whose output exists is skipped; `--redo` runs it again. A path given to one
step is stored in `work/<plot>/run.json`, so later steps need only `--plot` and `--project`.

| step | needs | writes |
|---|---|---|
| `rayroots-plot segment <las>` | plot LAS | `work/<plot>/`: `<plot>_segmented.ply`, `<plot>_mesh.ply`, `<plot>.rct.log` |
| `rayroots-plot roots` | segmented cloud; `--mesh`, `--boundary` optional | `work/<plot>/`: `<stem>_roots.npz`, `<stem>_roots.csv` |
| `rayroots-plot select` | segmented cloud | `work/<plot>/trees_in_plot.csv` |
| `rayroots-plot review` | the above; roots optional | `review/<plot>/` |
| `rayroots-plot context [trees]` | review | `review/<plot>/context/` |
| `rayroots-plot tree <seg ids>` | `roots --keep` | `work/<plot>/tree/`: trees with `root_src` and `contested` |
| `rayroots-plot finish [--apply]` | `checklist.csv` | `finished_trees/<genus>/`, `tree_register.csv` |
| `rayroots-plot overview` | reviews | `overview.html` |

`scripts/run_plot.sh` chains segment, roots, select and review. Run it detached: `rayextract trees` peaked at 7.4 GB on the
22 M point crop, so a full plot needs a server (`scripts/server_template.sh`). The root step ran on 8 GB: 30 min on a 172 M
point plot with 4 CPUs.
```
setsid nohup bash scripts/run_plot.sh <plot>.las --boundary b.csv --project P > run.log 2>&1 &
```

**Your own rayextract output.** Skip `segment` and start from your file; add the terrain mesh if you have one:
```
rayroots-plot roots --plot X --project P --segmented own_segmented.ply --mesh own_mesh.ply --boundary b.csv
rayroots-plot select --plot X --project P
rayroots-plot review --plot X --project P
```
`run_plot.sh` does the same when it is given a `_segmented.ply`.

**RCT settings.** The `rct` section of `rayroots/plot/plot.yaml` is passed through as flags: `rayimport`, `terrain` and `trees`.
Empty = RCT's defaults, `true` = a bare flag. Copy what you change into your own yaml and pass it with `--settings`:
```yaml
rct:
  trees: {max_diameter: 0.9}
```
Each command goes to `<plot>.rct.log`. RCT inside the root method is separate: it is part of `roots.yaml`, and changing it marks
the run as not the frozen method.

**Without roots:** `--no-roots` for plots without prop roots. **Several plots:** `plots.csv`, columns `plot`, `cloud`, `boundary`,
`roots` yes/no, optional `plot_genus`; paths relative to the csv. A failed plot is logged and skipped, the overview comes last:
```
setsid nohup bash scripts/run_plot.sh --plots plots.csv --project P > run.log 2>&1 &
```

**Other settings.** `rayroots/plot/plot.yaml` (`--settings`). H, D and `hag` always use the ground grid. Review options:
`--plot-genus G`, `--genus-hints h.csv` (columns `tree` or `seg_id`, and `genus`), `--existing-trees e.csv` (columns `genus`,
`H_m`, `D_cm`: count as already covered in the priority). The genus suggestion from root volume is off by default
(`suggest_genus`); its threshold is calibrated on Suriname data.

**Root step.** Ground: `grid` (default, the frozen method) computes the ground from the cloud's lowest points; `mesh` uses the RCT
terrain mesh and needs `--mesh`. A mesh run is logged as not the frozen method. `--config` takes a yaml with changed method settings.

| exit code | meaning |
|---|---|
| 0 | done |
| 2 | bad input (the message says what) |
| 3 | memory guard: rerun the same command, it resumes from the tile cache (`run_plot.sh` does this up to 3 times) |
| 5 | outputs written, but flip tiles failed: those areas have no v1 roots. See the WARNING lines and the `tile_failed` flag |

While running, `tmp_rayroots_<stem>/` holds the caches; it is deleted at the end. A resume with other inputs or settings is refused.

## Outputs of part A
```
<project>/A_plot/
  overview.html
  tree_register.csv
  review/<plot>/        00_gallery.html, checklist.csv, trees/, previews/, context/
  finished_trees/<genus>/<tree>.ply
  work/<plot>/          everything else
```
| file | content |
|---|---|
| `review/<plot>/00_gallery.html` | previews of trees with H >= 3 m, per genus by priority |
| `review/<plot>/checklist.csv` | fill in `status` (ok / fix / drop), `genus`, `notes` |
| `review/<plot>/trees/<priority>_<tree>.ply` | own segment + roots: RGB as below, `root_src`, `contested`, `hag`. Without roots: x, y, z, `hag` |
| `review/<plot>/context/<tree>_context.ply` | from `context`: surroundings, field `seg_id` |
| `finished_trees/<genus>/<tree>.ply`, `tree_register.csv` | from `finish`: the input of part B |
| `work/<plot>/<stem>_roots.npz` | `fi` (point index), `seg_id`, `src`: one row per (point, segment) a method attaches |
| `work/<plot>/<stem>_roots.csv` | per segment: `n_rct`, `n_v1_only`, `n_v2_only`, `n_both`, `n_contested`, `H_m`, `base_x`, `base_y`, `dist_to_boundary`, `flags` |
| `work/<plot>/trees_in_plot.csv`, `trees.csv` | selected trees; everything per tree |
| `work/<plot>/trees_below_3m.csv`, `trees_edge_truncated.csv` | listed only, no tree file |
| `work/<plot>/review_log.csv`, `<plot>.rct.log`, `run.json` | `finish` actions; exact LAS offset, RCT commands and timings; inputs |
| `work/overview.csv`, `work/covered_trees.csv` | per plot; finished trees, usable as `--existing-trees` |

`--keep` also keeps the LAS copy, the ray cloud and RCT's tree mesh (`segment`), and the recoloured plot cloud
`<stem>_segmented_roots.ply`, `<stem>_unsegmented.ply` and one cloud per segment in `trees/` from `raysplit` (`roots`). These are
deleted by default: on the crop 2.7 of 3.6 GB.

`root_src`: 0 RCT, 1 v1 only, 2 v2 only, 3 both.

Edge-truncated = at least 200 points, and 0.2 % of the tree's points, in the cells along the cloud's outer edge. These
trees are cut off by the scan and get no tree file.

`flags`: `outside_boundary` (base more than 2 m outside the hull), `no_stem_track` (no stem found, often not a tree),
`absorbed` (no points left: all taken by other segments), `contested`, `tile_failed`.

**Contested points**: v1 gives the point to one tree, v2 to another (4 % of the root points of SU7_psp3). The npz
keeps both rows. A point in the plot PLY has one colour, so there it goes to v2's tree (v2 had the higher precision
below 2 m; contested points were not scored separately). `rayroots-plot tree` puts it in both trees.

A point can move away from its RCT segment only when that segment does not reach `h_cut` (or is unsegmented). So
small segments low down can lose points to a tree (`n_rct` counts them before that).

## Cleaning a tree in CloudCompare
Start from `review/<plot>/00_gallery.html` and `checklist.csv`. Trees marked `fix` need `rayroots-plot context --plot <plot>`,
then open `trees/<priority>_<tree>.ply` and `context/<tree>_context.ply` in CloudCompare, and save the result as
`trees/<tree>_clean.ply`. `finish` copies `ok` trees and `_clean` files. An `ok` tree with `n_contested` > 0 can share
points with a neighbour: `finish` warns, it does not block.

Outside a review, to make a file that shows where each point comes from (needs `roots --keep`):
```
rayroots-plot tree 364 136 --plot X --project P
```
This writes `<stem>_tree_<id>.ply` with its RCT points and all its roots, contested ones included, with fields
`root_src` and `contested` and RGB by source:

| colour | RGB | meaning |
|---|---|---|
| grey | 127,127,127 | RCT tree (`root_src` 0) |
| green | 0,150,60 | roots from both v1 and v2 (3) |
| orange | 235,104,52 | v1 only (1) |
| blue | 42,120,214 | v2 only (2) |
| magenta | 220,0,200 | contested (`contested` = 1); also in another tree |

1. Open the file; keep the default PLY field mapping. `root_src` and `contested` appear as scalar fields.
2. Check magenta points against the other tree they are in.
3. Remove false points with the Segment tool (scissors), working in small views around the stem.
4. To see only one class: Edit > Scalar fields > Filter by value on `root_src`.
5. If a tree has two segments (split stem), merge them with Edit > Merge after cleaning; drop one copy of a
   contested point that is in both.
6. Save as binary PLY.

A root that no method found: open `<stem>_unsegmented.ply` alongside, cut the root out and merge it in.

## Part B: single tree
Input: `A_plot/finished_trees/<genus>/<tree>.ply` and `A_plot/tree_register.csv`; trees with status ok or fix. A tree already
done is skipped; `--tree T` redoes T, `--redo` all. Not validated yet, see the [README](../README.md#known-limits-and-not-tested).

`split` removes noise (largest connected component), finds the ground `z_ref`, splits the tree and writes one QC figure;
`model` runs RCT and adds panel (d) to it.

```
<project>/B_trees/
  tree_volumes.csv      one row per tree, one per stem
  qc/<tree>.png         (a) slices, last root, cut; (b) blobs per height; (c) slices below the last root and at the cut; (d) cylinders
  overrides.csv         your corrections, optional
  work/                 standard/ (after noise removal), stems/, roots/, base/, profile/, split_log.csv, rct/<tree>/
```
`rct/<tree>/` holds RCT's `trees.txt` and the measured radii per part; roots are in the flipped frame. The ensemble's tree files are
deleted unless `model --keep`.

**RCT settings:** the `rct` section of `rayroots/single/single.yaml`, passed through as flags to `rayextract trees` (`stem`, `roots`);
`ray_vector` for rayimport. Pass your own yaml with `--settings`. The tests below hold for the defaults.

| column in `tree_volumes.csv` | meaning |
|---|---|
| `vol_rct_m3`, `vol_meas_m3` | central run; with `--ensemble` the median of the runs, with `_central_`, `_mean_`, `_sd_` columns |
| `vol_stems_*`, `vol_roots_*`, `vol_base_*` | the same per part of a tree |
| `DBH_h_field_m`, `DBH_h_m` | field height, and the clean height used. Rhizophora: 0.30 m above the last root; Avicennia: 1.3 m above the ground |
| `DBH_points_cm`, `DBH_rct_cm`, `DBH_meas_cm` | from the points, the RCT model, the measured radii. Several stems: `/` |
| `h_cut_m`, `height_m`, `z_ref`, `n_stems`, `n_runs` | cut and height above the ground, and counts |
| `flags` | see below |

| flag | meaning |
|---|---|
| `rct split`, `stems merged` | RCT made more or fewer trees than there are stems |
| `vol gap` | `vol_rct` and `vol_meas` of the tree differ by more than 25 % (`gap_max`) |
| `short run` | the clean stem run is shorter than 0.6 m: check the cut |
| `base off stem` | the base of the RCT tree is off the stem at the cut |
| `DBH moved ...`, `DBH points: ...`, `DBH: no clean height ...` | DBH rules |
| `small stem merged`, `fragment to roots`, `override` | split events |
| `fork near DBH`, `late start`, `sparse base` | Avicennia: fork within 0.15 m of 1.3 m, stem found high up, few points at the base |

`overrides.csv` in `B_trees/`: `tree`, `h_cut_above_z0` (Rhizophora, as the column in `split_log.csv`), `n_stems` (Avicennia),
`reason`. `z0` = the 0.1th percentile of the tree's z, in `split_log.csv`; it is not `z_ref`, the ground of `h_cut_m`.

**Known limits.** A central run can be well off where RCT splits a tree: on `SU8_PSP4_52` the central run gives 2.41 m3
and the ensemble median 3.43 m3. On 364_136 the central total is 2 % below the median. Without a hand-cleaned input the
cut can be far off: `crop_22` and `crop_56` from the crop are the uncleaned stems of 364_136 and were cut at 1.74 and
3.91 m, against 3.53 m on the cleaned tree. The ensemble sd covers precision only, not the shape bias of RCT.

**Tests.** Against reference logs, `tests/compare_single.py`: identical cut, stem count, point counts and central volumes on
`SU7_psp3_364_136` and `SU7_psp3_119`, and identical ensemble medians on 364_136 (10 min, 0.9 GB). The cut: [validation/cut.md](validation/cut.md).

## Method
The method is frozen: v1 and v2 are not changed, a change is a new version.

- **Ground grid**: the 20th lowest point per 0.5 m cell, outlier cells replaced, gaps filled, smoothed.
- **Region**: points 0.15 to 8.6 m above ground, within hull + 10 m. Segments within hull + 9 m and at least 2 m tall
  get a stem track and compete for roots, whether in the plot or not.
- **v1** (flipped RCT): below `h_cut` the cloud is mirrored in height and cut in tiles, so that `rayextract trees` grows
  from the stems at the cut towards the root tips. Points beyond 0.05 m of their owner's cylinder surface are dropped.
- **v2** (fragment graph): points are cut in 0.20 m block fragments with a measured width and axis. All stems grow at
  once over this graph with rules on thickness, direction, mats (fans) and splits; near-ties stay unassigned. Second pass:
  fragments that stopped growth on thickness are re-cut at 0.10 m, which separates touching roots.
- **Union**: a point is a root of segment s when v1 or v2 gives it s and RCT did not; `src` records which.

| parameter | default | key in roots.yaml |
|---|---:|---|
| height of the root zone | 4.0 m | `h_cut` |
| v1 cylinder tolerance | 0.05 m | `cyl_tol` |
| ground removed below | 0.15 m | `cuts.Z_MIN` |
| root reach | 7 m | `cuts.R_MAX` |
| trees of the plot: hull + | 2 m | `cuts.BUFFER` |
| smallest tree | 2 m | `cuts.MIN_HEIGHT` |
| v1 tile / overlap | 16 m / 7 m | `flip.TILE`, `flip.MARGIN` |
| v2 tile | 12 m | `v2.CORE` |
| v2 voxel / block | 0.025 m / 8 voxels | `connect.VOX`, `connect.BLOCK` |
| v2 thickness tolerance | 0.30 | `v2.TOL` |
| v2 largest bend / branch | 40 / 90 deg | `v2.TURN_MAX`, `v2.BRANCH_MAX` |

`rayroots/roots.yaml` lists all named constants with these values as defaults. Copy it, keep only what you change, and
pass it with `--config`. Any changed method value makes the run "not frozen" in the log and in the npz. A few inline
numbers inside the method code (for example the 0.3 m seed band below `h_cut`) are not in the yaml.

## Known limits
- **Constant `h_cut` of 4 m** for every tree, whatever its size. Roots attached higher are missed, and for small trees
  the zone includes stem and crown.
- **No species gate.** Every tree gets root proposals, Avicennia too. The method was developed for Rhizophora prop
  roots; in the original workflow the step only ran on plots with Rhizophora in the inventory. The quality of the
  proposals on Avicennia is not known.
- **Same scanner and settings.** Fixed lengths (voxels, raster) assume the point density of the Suriname scans.
- **Mesh ground** is an option on top of the frozen method. It was only compared with the grid on one crop, where
  it changed 18 % of the root rows (mesh on average 2 cm above the grid, `grid_vs_mesh.md`), and never validated against
  hand roots.
- **Segmentation errors stay**: rayroots does not merge or split RCT stems. A tree in two segments needs merging by hand.
- **Crown points**: precision above 2 m is low in every method, see [validation/roots.md](validation/roots.md).
- `raysplit` output has no `root_src`; use `rayroots-plot tree`.
- Results depend on the RCT build and the pinned numpy/scipy; another build gives other point sets, and the expected
  test values no longer apply.

## Validation
[Root proposals on unseen plots](validation/roots.md), [the cut on hand-cut trees](validation/cut.md).

**Not validated:**
- The union of v1 and v2, and the v2-wins rule for contested points: v1 and v2 were scored separately.
- Avicennia roots: the test trees are Rhizophora.
- `ground = mesh`.
- Any other site, scanner or scan setting.
- The whole-cloud mode (no boundary), where edge trees have one-sided data.
- Only 8 hand-root trees in 2 plots carry the root scores.

## Tests
`docs/make_figures.py` makes the figures in `docs/img/`. `pytest tests/test_config.py` checks that the yaml defaults equal the frozen constants. `tests/` also holds the
regression, the crop check and the comparisons; their data is not in git. Crop run:
```
python tests/make_crop.py <SU7_psp3>.las <crop>.las -13 7 -11 9
bash scripts/run_plot.sh <crop>.las --project P --plot crop --keep
python tests/crop_check.py P/A_plot/work/crop crop tests/crop_expected.json
```
The crop holds the hand-root tree SU7_psp3_364_136 (Rhizophora) and the reference tree SU7_psp3_119 (Avicennia).
`tests/compare_reference.py` scores the result against the hand-root tree; expected values in `tests/crop_reference.md`.
That tree was used in the method's development, so this is a known result to compare against, not a validation.
