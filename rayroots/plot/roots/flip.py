# Frozen method (v1 + v2). Do not edit here: a change is a new version.
"""Frozen v1 (method c_cyl): flipped RCT below H_CUT in tiles and the cylinder filter. v1 setting: H_CUT 4.0, CYL_TOL 0.05."""
import shutil
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree

from rayroots.plot.roots import rct_tools as rt
from rayroots.plot.roots import cuts
from rayroots.plot.roots import plotdata as cm

H_CUT, CYL_TOL = 4.0, 0.05      # frozen v1
Z_MIN = cuts.Z_MIN
SLAB = 0.10           # m below H_CUT: the slice whose own blobs become patches
PATCH_M = 0.05        # m added to a blob's radius (bases snap to 0.1 m grid vertices: steps 03/04)
GRID_STEP = 0.1       # m: patch triangle grid
MAP_TOL = 0.3         # m: an extracted tree's base may lie this far outside a patch
TILE, MARGIN = 16.0, 7.0   # m: tile core and overlap (roots reach < 7 m: root_window_check.py grid); 16 m keeps RCT inputs < ~10 M points
TREES = {"gravity_factor": 0, "global_taper_factor": 0, "girth_height_ratio": 0.001, "height_min": 0.5,
         "crop_length": 0.15, "max_diameter": 2}
CYL_STEP = 0.02       # m: axis sample spacing for point-to-cylinder distances (error <= 0.01 m)
CYL_K, CYL_MAX = 16, 0.5   # neighbours searched; distances beyond CYL_MAX m are stored as CYL_MAX
MEM_MIN_GB, MEM_WAIT = 2.0, 1800
CACHE = TILE_DIR = TREES_DIR = TMP = None   # set by lib.use_tmp (the step's tmp folder)


def mem_available_gb():
    for line in open("/proc/meminfo"):
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) / 2**20
    return float("nan")


def memory_guard(tag):
    """After a tile: pause while MemAvailable < MEM_MIN_GB; give up (exit 3) after MEM_WAIT s."""
    t0 = time.time()
    while mem_available_gb() < MEM_MIN_GB:
        if time.time() - t0 > MEM_WAIT:
            print(f"  {tag}: memory still < {MEM_MIN_GB} GB after {MEM_WAIT} s, exiting for a restart", flush=True)
            raise SystemExit(3)
        print(f"  {tag}: paused, {mem_available_gb():.2f} GB available", flush=True)
        time.sleep(60)


def blob_list(xy):
    """Own blobs in a slice: (centre, radius) per blob of >= MIN_CELLS cells (1 cell dilation)."""
    if len(xy) < 4:
        return []
    origin = xy.min(0) - 2 * cuts.RASTER
    shape = tuple(np.floor((xy.max(0) - origin) / cuts.RASTER).astype(int) + 3)
    lab = cuts.blobs(xy, origin, shape, 1)
    out = []
    for b in np.unique(lab[lab > 0]):
        c = (np.argwhere(lab == b) + 0.5) * cuts.RASTER + origin
        cen = c.mean(0)
        out.append((cen, float(np.linalg.norm(c - cen, axis=1).max() + cuts.RASTER / 2)))
    return out


def patches(r, hcut, tracks):
    """Patch table (seg_id, cx, cy, radius, kind) and per-segment class; continuity labels for 'above' segments."""
    ids, hag = r["ids"], r["hag"]
    seg = ids >= 0
    u, inv = np.unique(ids[seg], return_inverse=True)
    lo = np.full(len(u), np.inf); hi = np.full(len(u), -np.inf)
    np.minimum.at(lo, inv, hag[seg]); np.maximum.at(hi, inv, hag[seg])
    cls = pd.DataFrame({"seg_id": u, "hag_min": lo, "hag_max": hi})
    cls["class"] = np.where(cls.hag_max < hcut, "below", np.where(cls.hag_min >= hcut, "above", "reaches"))
    rows = []
    sl = seg & (hag >= hcut - SLAB) & (hag < hcut)
    o = np.argsort(ids[sl], kind="stable")
    sid_sl, xy_sl = ids[sl][o], r["xyz"][sl][o][:, :2].astype(np.float64)
    us, first, cnt = np.unique(sid_sl, return_index=True, return_counts=True)
    span = dict(zip(us.tolist(), zip(first.tolist(), (first + cnt).tolist())))
    reaches = set(cls.loc[cls["class"] == "reaches", "seg_id"].tolist())
    for s in reaches:
        a, b = span.get(s, (0, 0))
        bl = blob_list(xy_sl[a:b])
        for cen, rad in bl:
            rows.append({"seg_id": s, "cx": cen[0], "cy": cen[1], "radius": rad + PATCH_M, "kind": "slice blob"})
        if not bl and s in tracks:                       # reaches H_CUT but no own blob in the slab: track
            c, rr = tracks[s]
            ok = np.flatnonzero(np.isfinite(rr))
            if len(ok):
                k = ok[np.argmin(abs(cuts.HEIGHTS[ok] - hcut))]
                rows.append({"seg_id": s, "cx": c[k, 0], "cy": c[k, 1], "radius": rr[k] + PATCH_M, "kind": "track (gap)"})
    cont = np.full(len(ids), -2, np.int32)               # continuity labels for -1 stem points above H_CUT
    for s in cls.loc[cls["class"] == "above", "seg_id"]:
        if s not in tracks:
            continue
        c, rr = tracks[s]
        ok = np.flatnonzero(np.isfinite(rr))
        if not len(ok):
            continue
        k = ok.min()
        cen, rad = c[k], rr[k]
        rows.append({"seg_id": s, "cx": cen[0], "cy": cen[1], "radius": rad + PATCH_M, "kind": "track (starts above)"})
        bottom = cls.loc[cls.seg_id == s, "hag_min"].iloc[0]
        m = (ids < 0) & (hag >= hcut) & (hag < bottom)
        m[m] = np.linalg.norm(r["xyz"][m, :2].astype(np.float64) - cen, axis=1) <= rad + 0.1
        cont[m] = s
    return pd.DataFrame(rows), cls, cont


def write_ply_xyz(path, xyz):
    v = np.ascontiguousarray(xyz, dtype="<f8")
    head = ("ply\nformat binary_little_endian 1.0\nelement vertex %d\nproperty double x\nproperty double y\n"
            "property double z\nend_header\n" % len(v))
    with open(path, "wb") as f:
        f.write(head.encode()); f.write(v.tobytes())


def write_patches(path, centres, radii, z, step=GRID_STEP):
    """Triangle grid restricted to discs (centre, radius): a face is kept when all its corners are in a disc."""
    lo = centres.min(0) - radii.max() - step
    hi = centres.max(0) + radii.max() + step
    xs, ys = np.arange(lo[0], hi[0] + step, step), np.arange(lo[1], hi[1] + step, step)
    X, Y = np.meshgrid(xs, ys)
    P = np.c_[X.ravel(), Y.ravel()]
    ins = np.zeros(len(P), bool)
    kd = cKDTree(P)
    for c, rr in zip(centres, radii):
        ins[kd.query_ball_point(c, max(rr, step * 0.75))] = True
    nx = len(xs)
    j, i = np.meshgrid(np.arange(len(ys) - 1), np.arange(nx - 1), indexing="ij")
    a = (j * nx + i).ravel(); b = a + 1; c = a + nx + 1; e = a + nx
    tri = np.r_[np.c_[a, b, c], np.c_[a, c, e]]
    tri = tri[ins[tri].all(1)]
    used = np.unique(tri)
    remap = np.full(len(P), -1); remap[used] = np.arange(len(used))
    v = np.c_[P[used], np.full(len(used), z)].astype("<f4")
    faces = np.c_[np.full(len(tri), 3), remap[tri]].astype("<i4")
    head = (f"ply\nformat binary_little_endian 1.0\nelement vertex {len(v)}\nproperty float x\nproperty float y\n"
            f"property float z\nelement face {len(faces)}\nproperty list int int vertex_indices\nend_header\n")
    with open(path, "wb") as f:
        f.write(head.encode()); f.write(v.tobytes()); f.write(faces.tobytes())
    return len(tri)


def cylinder_samples(path):
    """Axis samples (every CYL_STEP m) of every segment in a rayextract trees.txt: xyz, radius, base node (xyz of
    the tree's first node) per sample."""
    pts, rad, base = [], [], []
    for line in open(path):
        if line.startswith("#") or line.startswith("x,"):
            continue
        a = np.array([[float(v) for v in t.strip(",").split(",")] for t in line.split()])
        if a.ndim != 2 or not len(a):
            continue
        xyz, r, par = a[:, :3], a[:, 3], a[:, 4].astype(int)
        for i in np.flatnonzero(par >= 0):
            p0, p1 = xyz[par[i]], xyz[i]
            n = max(2, int(np.ceil(np.linalg.norm(p1 - p0) / CYL_STEP)) + 1)
            t = np.linspace(0, 1, n)[:, None]
            pts.append(p0 + t * (p1 - p0)); rad.append(np.full(n, r[i])); base.append(np.repeat(xyz[:1], n, 0))
    if not pts:
        return np.zeros((0, 3)), np.zeros(0), np.zeros((0, 3))
    return np.vstack(pts), np.concatenate(rad), np.vstack(base)


def run_tile(xy, zf, pt, keep_trees=None):
    """rayextract on flipped points; returns owner patch index per point (-1 none), distance beyond the surface of
    the owner segment's cylinders (NaN: no owner), distance beyond the nearest cylinder of any tree, and info.
    keep_trees: path to copy the flipped trees.txt to."""
    work = Path(tempfile.mkdtemp(dir=TMP))
    try:
        pts = np.c_[xy, zf]
        write_ply_xyz(work / "w.ply", pts)
        n_tri = write_patches(work / "patches.ply", pt[["cx", "cy"]].values, pt.radius.values,
                              zf.min() + rt.GROUND_OFFSET)
        err = rt.rayimport(work / "w.ply", work / "w_rc.ply")
        if not err:
            err = rt.rayextract_trees(work / "w_rc.ply", work / "patches.ply", TREES, str(work / "out"),
                                      keep=("segmented.ply", "trees.txt"))
        if err:
            nan = np.full(len(xy), np.nan)
            return np.full(len(xy), -1), nan, nan, {"status": err, "faces": n_tri}
        s = cm.rw.ply_memmap(work / "out_segmented.ply")
        sp = np.c_[s["x"], s["y"], s["z"]].astype(np.float64)
        tid = cm.rw.decode_ids(s["red"], s["green"], s["blue"])
        # tree -> patch by its lowest points (not by colour order)
        zplane = zf.min()
        owner = {}
        pkd = cKDTree(pt[["cx", "cy"]].values)
        for t in np.unique(tid[tid >= 0]):
            q = sp[tid == t]
            low = q[q[:, 2] <= q[:, 2].min() + 0.1]
            c = low[:, :2].mean(0)
            d, k = pkd.query(c)
            if d <= pt.radius.values[k] + MAP_TOL and q[:, 2].min() <= zplane + 0.3:
                owner[t] = k
        lab = np.array([owner.get(t, -1) for t in tid])
        d, k = cKDTree(sp).query(pts, distance_upper_bound=1e-3)
        out = np.full(len(xy), -1)
        ok = np.isfinite(d)
        out[ok] = lab[k[ok]]
        # cylinders: samples along every segment; a tree's owner patch from its base node (same rule as points)
        cs, cr, cb = cylinder_samples(work / "out_trees.txt")
        d_own = np.full(len(xy), np.nan); d_any = np.full(len(xy), np.nan)
        if len(cs):
            bd, bk = pkd.query(cb[:, :2])
            cown = np.where(bd <= pt.radius.values[bk] + MAP_TOL, pt.seg_id.values[bk], -1)
            ckd = cKDTree(cs)
            pseg = np.where(out >= 0, pt.seg_id.values[np.maximum(out, 0)], -1)
            for s0 in range(0, len(pts), 1_000_000):
                sl = slice(s0, s0 + 1_000_000)
                dd, kk = ckd.query(pts[sl], k=CYL_K, distance_upper_bound=CYL_MAX)
                fin = np.isfinite(dd)
                kc = np.minimum(kk, len(cr) - 1)
                surf = np.where(fin, dd - cr[kc], np.inf)
                da = surf.min(1)
                same = fin & (cown[kc] == pseg[sl][:, None])
                do = np.where(same, surf, np.inf).min(1)
                d_any[sl] = np.where(np.isfinite(da), da, CYL_MAX)
                d_own[sl] = np.where(out[sl] >= 0, np.where(np.isfinite(do), do, CYL_MAX), np.nan)
        if keep_trees is not None:
            shutil.copy(work / "out_trees.txt", keep_trees)
        return out, d_own, d_any, {"status": "ok", "faces": n_tri, "trees": int((np.unique(tid) >= 0).sum()),
                                   "attached_trees": len(owner), "matched": float(ok.mean()), "cylinder_samples": len(cs)}
    finally:
        shutil.rmtree(work, ignore_errors=True)


def run(plot, hcut):
    out_npz = CACHE / f"{plot}_flip_h{hcut}.npz"
    if out_npz.exists():
        return
    TMP.mkdir(parents=True, exist_ok=True)
    TREES_DIR.mkdir(parents=True, exist_ok=True)
    TILE_DIR.mkdir(parents=True, exist_ok=True)
    print(f"{plot} h{hcut}: start, {mem_available_gb():.2f} GB available", flush=True)
    t_all = time.time()
    r = cuts.assign_region(plot)
    tracks = cuts.stems(plot)
    pt, cls, cont = patches(r, hcut, tracks)
    below_ids = set(cls.loc[cls["class"] == "below", "seg_id"].tolist())
    ids, hag = r["ids"], r["hag"]
    re_ok = (ids < 0) | np.isin(ids, list(below_ids))            # re-assignable
    inreg = (hag > Z_MIN) & (hag < hcut)
    xy = r["xyz"][:, :2]
    new = ids.astype(np.int32).copy()
    new[cont >= 0] = cont[cont >= 0]
    dist_own = np.full(len(ids), np.nan, np.float16)    # re-assigned points: beyond the owner's cylinder surface
    dist_any = np.full(len(ids), np.nan, np.float16)    # region points in a tile core: beyond any cylinder
    lo = pt[["cx", "cy"]].values.min(0) - MARGIN
    hi = pt[["cx", "cy"]].values.max(0) + MARGIN
    nx, ny = np.ceil((hi - lo) / TILE).astype(int)
    tiles = []
    for i in range(nx):
        for j in range(ny):
            c0 = lo + np.array([i, j]) * TILE
            core = inreg & np.all((xy >= c0) & (xy < c0 + TILE), axis=1)
            if not (core & re_ok).any():
                continue
            m = inreg & np.all((xy >= c0 - MARGIN) & (xy < c0 + TILE + MARGIN), axis=1)
            pin = pt[np.all((pt[["cx", "cy"]].values >= c0 - MARGIN) & (pt[["cx", "cy"]].values < c0 + TILE + MARGIN), axis=1)]
            if not len(pin):
                continue
            tfile = TILE_DIR / f"{plot}_h{hcut}_tile{i}_{j}.npz"
            if tfile.exists():                           # finished in an earlier (killed) run
                z = np.load(tfile, allow_pickle=True)
                new[z["take"]] = z["sid"]
                dist_own[z["take"]] = z["d_own"]
                dist_any[z["core"]] = z["d_any"]
                info = z["info"].item()
                info["resumed"] = True
                tiles.append(info)
                print(f"  {plot} h{hcut} tile {i},{j}: from tile cache", flush=True)
                continue
            t0 = time.time()
            w = np.flatnonzero(m)
            own, d_own, d_any, info = run_tile(xy[w].astype(np.float64), 2 * hcut - hag[w].astype(np.float64),
                                               pin.reset_index(drop=True),
                                               keep_trees=TREES_DIR / f"{plot}_h{hcut}_tile{i}_{j}_trees.txt")
            sid = np.where(own >= 0, pin.seg_id.values[np.maximum(own, 0)], -1)
            take = core[w] & re_ok[w] & (sid >= 0)
            new[w[take]] = sid[take]
            c_ = core[w]
            dist_own[w[take]] = d_own[take]
            dist_any[w[c_]] = d_any[c_]
            info.update({"plot": plot, "hcut": hcut, "tile": f"{i},{j}", "points": len(w), "patches": len(pin),
                         "assigned": int(take.sum()), "seconds": time.time() - t0,
                         "mem_avail_gb": round(mem_available_gb(), 2)})
            tmp = tfile.with_suffix(".part.npz")
            np.savez(tmp, take=w[take].astype(np.int64), sid=sid[take].astype(np.int32),
                     d_own=d_own[take].astype(np.float16), core=w[c_].astype(np.int64),
                     d_any=d_any[c_].astype(np.float16), info=np.array(info, dtype=object))
            tmp.replace(tfile)                           # atomic: a kill never leaves a half-written tile
            tiles.append(info)
            print(f"  {plot} h{hcut} tile {i},{j}: {info}", flush=True)
            memory_guard(f"{plot} h{hcut} tile {i},{j}")
    # root check per segment: newly attached points below H_CUT
    gained = (new != ids) & inreg & (new >= 0)
    rows = []
    g_ids = new[gained]
    gxy = xy[gained].astype(np.float64); ghag = hag[gained]
    for s in cls.seg_id:
        m = g_ids == s
        p = pt[pt.seg_id == s]
        row = {"seg_id": s, "class": cls.loc[cls.seg_id == s, "class"].iloc[0], "patches": len(p),
               "patch_kind": ";".join(sorted(set(p.kind))), "new_points": int(m.sum()),
               "continuity_points": int((cont == s).sum())}
        if m.any() and len(p):
            d = cKDTree(p[["cx", "cy"]].values).query(gxy[m])[0]
            row["spread_p95"] = float(np.percentile(d, 95))
            low = gxy[m][ghag[m] < 0.5]
            if len(low):
                ij = np.floor((low - low.min(0)) / 0.2).astype(int)
                grid = np.zeros(ij.max(0) + 1, bool); grid[ij[:, 0], ij[:, 1]] = True
                row["ground_contacts"] = int(ndimage.label(grid, np.ones((3, 3)))[1])
            else:
                row["ground_contacts"] = 0
        rows.append(row)
    seg = pd.DataFrame(rows)
    before = pd.Series(ids[ids >= 0]).value_counts(); after = pd.Series(new[new >= 0]).value_counts()
    seg["absorbed"] = seg.seg_id.map(lambda s: s in before.index and s not in after.index)
    absorbed_into = {}
    for s in seg.loc[seg.absorbed, "seg_id"]:
        absorbed_into[s] = int(pd.Series(new[ids == s]).mode().iloc[0])
    seg["absorbed_into"] = seg.seg_id.map(absorbed_into)
    seg.to_csv(CACHE / f"{plot}_flip_h{hcut}.csv", index=False, float_format="%.3f")
    pd.DataFrame(tiles).to_csv(CACHE / f"{plot}_flip_h{hcut}_tiles.csv", index=False)
    np.savez(out_npz, lab=new, dist_own=dist_own, dist_any=dist_any, reassignable=re_ok & inreg)
    print(f"{plot} h{hcut}: {len(tiles)} tiles, {time.time() - t_all:.0f} s, gained {int(gained.sum()):,}, "
          f"absorbed {int(seg.absorbed.sum())}", flush=True)



def cyl_labels(ids, z, tol):
    """c_cyl: re-assigned points beyond tol of their owner's cylinder surface -> -1."""
    lab = z["lab"].copy()
    re = z["reassignable"] & (lab != ids) & (lab >= 0)
    far = re & ~(z["dist_own"].astype(np.float32) <= tol)
    lab[far] = -1
    return lab

