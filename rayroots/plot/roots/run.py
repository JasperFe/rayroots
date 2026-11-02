"""The root step: rayroots-plot roots.
Exit codes: 0 done; 2 bad input; 3 memory guard (rerun: resumes from the tile caches); 5 finished, but tiles failed."""
import argparse
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

from rayroots.plot.roots import config, cuts, flip, output, v2
from rayroots.plot.roots import plotdata as cm

EXIT_TILES = 5


def use_tmp(d):
    """Point every cache of the modules at the run's tmp folder d (deleted when the run ends)."""
    d.mkdir(parents=True, exist_ok=True)
    cm.CACHE = cm.rw.CACHE = cuts.OUT = flip.CACHE = d
    flip.TILE_DIR, flip.TREES_DIR, flip.TMP = d / "flip_tiles", d / "flip_trees", d / "tmp_rct"


def check_tmp(tmp, stamp):
    """A resume must use the same inputs and settings as the run that made the tile caches."""
    f = tmp / "run_stamp.txt"
    if f.exists() and f.read_text() != stamp:
        raise SystemExit(f"{tmp} holds a killed run with other inputs or settings; delete it or use them again")
    f.write_text(stamp)


def check_inputs(plot):
    """Fail early when the boundary is not in the cloud's frame, or no tree reaches h_cut (rc 2)."""
    hull = cm.hull_xy(plot)
    g = cm.ground(plot)
    lo = np.array([g.x0, g.y0])
    hi = lo + np.array(g.g.shape) * cm.rw.GROUND_CELL
    if (hull.max(0) < lo).any() or (hull.min(0) > hi).any():
        print(f"boundary {hull.min(0).round(2)} .. {hull.max(0).round(2)} lies outside the cloud "
              f"{lo.round(2)} .. {hi.round(2)}: it must be in the segmented frame (LAS minus the exact header offset)")
        raise SystemExit(2)


def failed_tiles(plot, hcut):
    """Failed flip tiles: [(tile name, status, box lo, box hi)]; box = the tile with its margin, in x y."""
    p = flip.CACHE / f"{plot}_flip_h{hcut}_tiles.csv"
    t = pd.read_csv(p) if p.exists() and p.stat().st_size > 1 else pd.DataFrame()
    if t.empty:
        return []
    bad = t[t.status != "ok"]
    if bad.empty:
        return []
    r = cuts.assign_region(plot)
    pt = flip.patches(r, hcut, cuts.stems(plot))[0]
    lo = pt[["cx", "cy"]].values.min(0) - flip.MARGIN
    out = []
    for q in bad.itertuples():
        i, j = (int(v) for v in q.tile.split(","))
        c0 = lo + np.array([i, j]) * flip.TILE
        out.append((q.tile, q.status, c0 - flip.MARGIN, c0 + flip.TILE + flip.MARGIN))
    return out


def run(seg, mesh=None, boundary=None, out=Path("."), cfg_path=None, ground=None):
    t0 = time.time()
    cfg = config.load(cfg_path)
    if ground:
        cfg["ground"] = ground
    if cfg["ground"] not in ("grid", "mesh"):
        print("ground must be grid or mesh")
        return 2
    if cfg["ground"] == "mesh" and mesh is None:
        print("ground = mesh needs --mesh <stem>_mesh.ply")
        return 2
    for p in (seg, mesh, boundary):
        if p is not None and not p.is_file():
            print(f"not found: {p}")
            return 2
    if cfg["ground"] == "grid" and mesh is not None:
        print("ground = grid: --mesh not used")
    config.apply(cfg)
    frozen = config.is_frozen(cfg)
    h = config.digest(cfg)
    print(f"{seg.name}: settings {h}, " + ("frozen method (v1 + v2, ground grid)" if frozen else
          "NOT the frozen method (changed settings or ground = mesh)"), flush=True)
    stem = seg.stem.removesuffix("_segmented")
    plot = stem
    out.mkdir(parents=True, exist_ok=True)
    tmp = out / f"tmp_rayroots_{stem}"
    use_tmp(tmp)
    check_tmp(tmp, f"{seg.resolve()} {seg.stat().st_size} {boundary} {mesh} {h}\n")
    cm.setup(seg, boundary, mesh if cfg["ground"] == "mesh" else None)
    check_inputs(plot)
    cuts.assign_region(plot)                          # cached in tmp; flip and v2 load it themselves
    cuts.stems(plot)
    r = cuts.assign_region(plot)
    if flip.patches(r, flip.H_CUT, cuts.stems(plot))[0].empty:
        print(f"{seg.name}: no segment reaches {flip.H_CUT} m above the ground: nothing to attach roots to")
        return 2
    del r
    flip.run(plot, flip.H_CUT)
    if not v2.labels_path(plot).exists():             # resume: v2 finished before a kill
        v2.run(plot)
    r = cuts.assign_region(plot)
    z = np.load(tmp / f"{plot}_flip_h{flip.H_CUT}.npz")
    lab1 = flip.cyl_labels(r["ids"], z, flip.CYL_TOL)
    del z
    lab2 = np.load(v2.labels_path(plot))["lab"]
    assert len(lab1) == len(lab2) == len(r["ids"])
    fi, sid, src = output.roots_table(r["ids"], r["idx"], lab1, lab2)
    del r, lab1, lab2
    failed = failed_tiles(plot, flip.H_CUT)
    for name, status, _, _ in failed:
        print(f"WARNING: flip tile {name} failed, no v1 roots from it: {status}", flush=True)
    n_after = output.write_ply(out / f"{stem}_segmented_roots.ply", seg, fi, sid, src)
    output.write_csv(out / f"{stem}_roots.csv", plot, fi, sid, src, n_after, [(a, b) for _, _, a, b in failed])
    output.write_npz(out / f"{stem}_roots.npz", fi, sid, src, seg, h, frozen)
    shutil.rmtree(tmp)
    n = np.bincount(src, minlength=4)
    nc = len(np.unique(fi[output.contested_rows(fi)]))
    print(f"{stem}: {len(fi):,} root rows ({len(np.unique(sid))} segments): v1 only {n[1]:,}, v2 only {n[2]:,}, "
          f"both {n[3]:,}; contested points {nc:,}; failed tiles {len(failed)}; {time.time() - t0:.0f} s", flush=True)
    return EXIT_TILES if failed else 0


def step(plot, project=".", seg=None, mesh=None, boundary=None, cfg_path=None, ground=None, redo=False, keep=False):
    """rayroots-plot roots: run() in A_plot/work/<plot>/. Keeps roots.npz and roots.csv; the recoloured plot cloud only
    with --keep, then also split into one file per tree with raysplit (work/<plot>/trees/)."""
    from rayroots.plot import settings
    r = settings.inputs(project, plot, segmented=seg, mesh=mesh, boundary=boundary)
    seg = settings.need(r, "segmented", "roots --segmented")
    w = settings.work_dir(project, plot)
    stem = seg.stem.removesuffix("_segmented")
    npz = w / f"{stem}_roots.npz"
    if npz.exists() and not redo:
        print(f"{plot}: roots exist, skipped (--redo runs them again)")
        return 0
    rc = run(seg, Path(r["mesh"]) if "mesh" in r else None, Path(r["boundary"]) if "boundary" in r else None, w, cfg_path, ground)
    if rc in (0, EXIT_TILES):
        settings.inputs(project, plot, roots=npz)
        ply = w / f"{stem}_segmented_roots.ply"
        if not keep:
            ply.unlink(missing_ok=True)
        elif rc == 0:
            subprocess.run(["raysplit", ply.name, "seg_colour"], cwd=w, check=True, stdout=subprocess.DEVNULL)
            (w / "trees").mkdir(exist_ok=True)
            (w / f"{stem}_segmented_roots_-1.ply").rename(w / f"{stem}_unsegmented.ply")
            for p in w.glob(f"{stem}_segmented_roots_[0-9]*.ply"):
                p.rename(w / "trees" / p.name)
    return rc
