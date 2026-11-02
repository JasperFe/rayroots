# Frozen method (v1 + v2). Do not edit here: a change is a new version.
"""Frozen v2: plot-wide labels in tiles with the fragment graph and growth of connect.py."""
import time

import numpy as np

from rayroots.plot.roots import connect as cn
from rayroots.plot.roots import cuts, flip
from rayroots.plot.roots import plotdata as cm

VERSION = "v2"
H_CUT = 4.0
CORE = 12.0                 # m: tile core (a 28 m tile below H_CUT holds <= ~15 M points in the densest plots)
PAD = cn.R_REACH + 1.0      # m: context around the core
COL_PAD = 0.05                 # m added to the stem track radius: stem column (seed fragments)
DEFAULT = (0.30, 40.0, 90.0)    # prototype values (TOL, TURN_MAX), BRANCH_MAX


def labels_path(plot):
    return cm.CACHE / f"{plot}_v2_labels.npz"



def run(plot):
    t0 = time.time()
    cn.TOL, cn.TURN_MAX, cn.BRANCH_MAX = DEFAULT
    r = cuts.assign_region(plot)
    tracks = cuts.stems(plot)
    pt, cls, _ = flip.patches(r, H_CUT, tracks)
    reaching = set(cls.loc[cls["class"] == "reaches", "seg_id"].tolist())
    lab = np.where(np.isin(r["ids"], list(reaching)), r["ids"], -1).astype(np.int32)
    low = np.flatnonzero(r["hag"] < H_CUT)
    xy = r["xyz"][low, :2].astype(np.float64)
    lo, hi = xy.min(0), xy.max(0)
    nx, ny = np.ceil((hi - lo) / CORE).astype(int)
    idx = cuts.CellIndex(xy)
    n_tiles = n_new = n_ties = 0
    tdir = cm.CACHE / "v2_tiles" / plot                  # per-tile results: a killed run resumes here
    tdir.mkdir(parents=True, exist_ok=True)
    for i in range(nx):
        for j in range(ny):
            tp = tdir / f"t_{i}_{j}.npz"
            if tp.exists():
                z = np.load(tp)
                lab[z["pts"]] = z["lab"]
                n_tiles += int(z["used"]); n_new += len(z["pts"]); n_ties += int(z["ties"])
                continue
            c_lo = lo + np.array([i, j]) * CORE
            c_hi = c_lo + CORE
            box = idx.box(c_lo - PAD, c_hi + PAD)
            empty = dict(pts=np.zeros(0, np.int64), lab=np.zeros(0, np.int32), used=0, ties=0)
            if not len(box):
                np.savez(tp, **empty)
                continue
            box.sort()
            p = low[box]
            core = ((xy[box] >= c_lo) & (xy[box] < c_hi)).all(1)
            if not core.any():
                np.savez(tp, **empty)
                continue
            pid = r["ids"][p]
            reass = ~np.isin(pid, list(reaching))
            if not reass[core].any():
                np.savez(tp, **empty)
                continue
            pin = pt[((pt[["cx", "cy"]].values >= c_lo - PAD) & (pt[["cx", "cy"]].values < c_hi + PAD)).all(1)]
            ctx = {"xyz": r["xyz"][p].astype(np.float64), "hag": r["hag"][p], "pid": pid, "reass": reass,
                   "pin": pin, "tracks": tracks}
            B1 = build(ctx)
            if not B1["seeds"]:
                np.savez(tp, **empty)
                continue
            B2, gr, own, _ = two_pass(ctx, B1)
            po = own[B2["frag"]]
            m = core & reass & (po >= 0)
            lab[p[m]] = po[m]
            np.savez(tp, pts=p[m], lab=po[m].astype(np.int32), used=1, ties=int(gr.tie_sub.sum()))
            n_tiles += 1
            n_new += int(m.sum())
            n_ties += int(gr.tie_sub.sum())
            print(f"{plot} tile {i},{j}: {len(p):,} pts, {B1['nf']:,} fragments, {len(B1['seeds'])} stems, "
                  f"{int(m.sum()):,} new, {time.time() - t0:.0f} s", flush=True)
    np.savez(labels_path(plot), lab=lab, idx=r["idx"])
    info = {"plot": plot, "version": VERSION, "region_points": len(lab), "tiles": n_tiles,
            "reaching_segments": len(reaching), "new_points": n_new, "tie_fragments": n_ties,
            "seconds": round(time.time() - t0)}
    print(info, flush=True)
    return info



def seeds_for(s, pt_s, tracks, fh, fc, rct_owner):
    """Seed fragments of stem s: fragments in its tracked stem column (track radius + COL_PAD at the fragment's
    height) and, in the 0.3 m below H_CUT, in its patches; not owned by another segment's RCT points."""
    free = (rct_owner == s) | (rct_owner < 0)
    m = np.zeros(len(fh), bool)
    for q in pt_s.itertuples():
        m |= (fh >= H_CUT - 0.3) & (np.linalg.norm(fc[:, :2] - [q.cx, q.cy], axis=1) <= q.radius)
    diam = 2 * pt_s.radius.max()
    if s in tracks:
        cen, rad = tracks[s]
        k = np.clip(np.round((fh - cuts.H_START) / cuts.STEP).astype(int), 0, len(rad) - 1)
        ok = np.isfinite(rad[k]) & (fh >= cuts.H_START - 0.05)
        d = np.linalg.norm(fc[:, :2] - np.nan_to_num(cen[k]), axis=1)
        m |= ok & (d <= np.nan_to_num(rad[k]) + COL_PAD)
        below = np.isfinite(rad) & (cuts.HEIGHTS < H_CUT)
        if below.any():
            diam = max(diam, 2 * np.nanmax(rad[below]))
    return np.flatnonzero(m & free), diam



def build(ctx, fine=None):
    """Fragments, measurements, RCT owner per fragment, adjacency, seeds of the competing stems."""
    xyz, reass, pid = ctx["xyz"], ctx["reass"], ctx["pid"]
    frag, nf, pairs, _ = cn.fragments(xyz, fine)
    f = cn.measure(xyz, frag, nf)
    fh = np.bincount(frag, ctx["hag"], nf) / f["n"]
    k = ~reass
    u, cnt = np.unique(np.c_[frag[k], pid[k]], axis=0, return_counts=True)
    rct_owner = np.full(nf, -1, np.int64)
    big = cnt >= 0.5 * f["n"][u[:, 0]]
    rct_owner[u[big, 0]] = u[big, 1]
    seeds, diam, sxy = {}, {}, {}
    for s, g in ctx["pin"].groupby("seg_id"):
        sd, dm = seeds_for(int(s), g, ctx["tracks"], fh, f["c"], rct_owner)
        if len(sd):
            seeds[int(s)], diam[int(s)] = sd, dm
            sxy[int(s)] = g[["cx", "cy"]].values.mean(0)
    return {"frag": frag, "nf": nf, "f": f, "fh": fh, "rct_owner": rct_owner, "nb": cn.adjacency(nf, pairs),
            "seeds": seeds, "diam": diam, "sxy": sxy}


def grow(B):
    gr = cn.Growth(B["f"], B["nb"], B["seeds"], B["diam"], B["sxy"], B["rct_owner"]).run()
    return gr, gr.final_owner()


def two_pass(ctx, B1):
    gr1, _ = grow(B1)
    th = cn.thick_frags(gr1)
    seed = np.zeros(B1["nf"], bool)
    for v in B1["seeds"].values():
        seed[v] = True
    th = th[~seed[th] & (B1["rct_owner"][th] < 0)]
    B2 = build(ctx, np.isin(B1["frag"], th))
    gr2, own2 = grow(B2)
    return B2, gr2, own2, len(th)

