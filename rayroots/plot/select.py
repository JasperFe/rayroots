"""select: the trees of one plot from its segmented cloud -> work/<plot>/trees_in_plot.csv.
A tree = an RCT segment with a stem track (height >= 2 m) whose base lies within boundary + buffer; without a boundary
all segments with a stem track."""
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from rayroots.plot.roots import run as roots_run
from rayroots.plot.roots import config, cuts
from rayroots.plot.roots import plotdata as cm
from rayroots.plot import settings


def track_at(tracks, s, h, d_top, d_win):
    """(centre at height h, rough D in m) of segment s's stem track, or (nan, nan). D = the narrowest d_win running
    median of the track diameter between h and d_top: above the prop roots for Rhizophora (whose slices at 1.3 m
    merge stem and roots), about DBH for Avicennia (little taper over 1.3-4 m)."""
    if s not in tracks:
        return np.full(2, np.nan), np.nan
    c, r = tracks[s]
    ok = np.flatnonzero(np.isfinite(r))
    if not len(ok):
        return np.full(2, np.nan), np.nan
    k = ok[np.argmin(abs(cuts.HEIGHTS[ok] - h))]
    hs = cuts.HEIGHTS[ok]
    w = [np.median(r[ok[abs(hs - z) <= d_win / 2]]) for z in hs[(hs >= h) & (hs <= d_top)]]
    return c[k], 2 * float(min(w)) if w else np.nan


def run(plot, project=".", seg=None, boundary=None, settings_path=None, redo=False):
    r = settings.inputs(project, plot, segmented=seg, boundary=boundary)
    seg = settings.need(r, "segmented", "select --segmented")
    boundary = Path(r["boundary"]) if "boundary" in r else None
    out = settings.work_dir(project, plot)
    if (out / "trees_in_plot.csv").exists() and not redo:
        print(f"{plot}: trees_in_plot.csv exists, skipped (--redo runs it again)")
        return 0
    cfg = settings.load(settings_path)
    p = cfg["select"]
    root_cfg = config.load()
    root_cfg["cuts"]["BUFFER"] = p["buffer_m"]
    config.apply(root_cfg)
    tmp = out / f"tmp_select_{plot}"
    roots_run.use_tmp(tmp)
    cm.setup(seg, boundary)
    st = cuts.segments(plot)
    tracks = cuts.stems(plot, cuts.region(plot))
    g = cm.ground(plot)
    rows = []
    for q in st[st.seg_id.isin(list(tracks))].itertuples():
        s = int(q.seg_id)
        if not q.in_plot:
            continue
        cen, d = track_at(tracks, s, p["h_dbh_m"], p["d_top_m"], p["d_win_m"])
        usable = np.isfinite(d) and 100 * d <= cfg["d_max_cm"]
        base_z = float(g(np.array([q.base_x]), np.array([q.base_y]))[0])
        rows.append({"plot": plot, "seg_id": s, "tree": f"{plot}_{s}", "base_x": q.base_x, "base_y": q.base_y,
                     "base_z": base_z, "x": cen[0] if usable else q.base_x, "y": cen[1] if usable else q.base_y,
                     "H_m": q.z_max - base_z, "D_cm": 100 * d, "dist_to_boundary": q.dist_to_hull,
                     "n_points": q.n_points})
    t = pd.DataFrame(rows)
    part = out / "trees_in_plot.part"
    t.to_csv(part, index=False, float_format="%.3f")
    part.replace(out / "trees_in_plot.csv")
    shutil.rmtree(tmp)
    print(f"{plot}: {len(t)} trees in the plot (base within boundary + {p['buffer_m']:g} m)", flush=True)
    return 0
