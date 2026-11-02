"""Bucket files of the plot's points (tmp), for the tree files of review."""
import json
import shutil
import time
from collections import OrderedDict

import numpy as np
import pandas as pd
from scipy import ndimage

from rayroots.plot.roots import plotdata as cm

DT = np.dtype([("x", "<f4"), ("y", "<f4"), ("z", "<f4"), ("id", "<i4"), ("fi", "<i8")])


def store(plot, d, bucket, margin):
    meta_p = d / "store.json"
    if meta_p.exists():
        return json.load(open(meta_p))
    t0 = time.time()
    bdir = d / "buckets"
    if bdir.exists():
        shutil.rmtree(bdir)
    bdir.mkdir(parents=True)
    hull = cm.hull_xy(plot)
    lo, hi = hull.min(0) - margin, hull.max(0) + margin
    K = 1 << 24
    bmin = np.full((K, 3), np.inf, np.float32); bmax = np.full((K, 3), -np.inf, np.float32)
    cmin = np.full(3, np.inf); cmax = np.full(3, -np.inf)
    for s, c, ids in cm.chunks(cm.seg_path(plot), ["x", "y", "z", "alpha", "red", "green", "blue"]):
        k = np.flatnonzero(c["alpha"] > 0)
        xyz = np.c_[c["x"][k], c["y"][k], c["z"][k]].astype(np.float32)
        cmin = np.minimum(cmin, xyz.min(0)); cmax = np.maximum(cmax, xyz.max(0))
        i = ids[k]
        sg = i >= 0
        for a in range(3):
            np.minimum.at(bmin[:, a], i[sg], xyz[sg, a]); np.maximum.at(bmax[:, a], i[sg], xyz[sg, a])
        m = (xyz[:, :2] >= lo).all(1) & (xyz[:, :2] < hi).all(1)
        rec = np.empty(int(m.sum()), DT)
        rec["x"], rec["y"], rec["z"] = xyz[m, 0], xyz[m, 1], xyz[m, 2]
        rec["id"], rec["fi"] = i[m], (k[m] + s)
        bij = np.floor((xyz[m, :2] - lo) / bucket).astype(np.int64)
        key = bij[:, 0] * 10000 + bij[:, 1]
        o = np.argsort(key, kind="stable")
        uk, first = np.unique(key[o], return_index=True)
        for kk, a, b in zip(uk, first, np.r_[first[1:], len(o)]):
            with open(bdir / f"b_{kk // 10000}_{kk % 10000}.bin", "ab") as f:
                rec[o[a:b]].tofile(f)
        print(f"  store {plot}: {s + len(c):,}", flush=True)
    u = np.flatnonzero(np.isfinite(bmin[:, 0]))
    pd.DataFrame({"seg_id": u, "xmin": bmin[u, 0], "ymin": bmin[u, 1], "zmin": bmin[u, 2], "xmax": bmax[u, 0],
                  "ymax": bmax[u, 1], "zmax": bmax[u, 2]}).to_csv(d / "seg_bbox.csv", index=False,
                                                                   float_format="%.3f")
    meta = {"plot": plot, "lo": lo.tolist(), "hi": hi.tolist(), "cloud_min": cmin.tolist(),
            "cloud_max": cmax.tolist(), "seconds": round(time.time() - t0)}
    json.dump(meta, open(meta_p, "w"))
    return meta


class Buckets:
    """Points of the stored buckets in a box; an LRU cache of bucket arrays (budget in bytes)."""

    def __init__(self, d, meta, bucket, budget=1.2e9):
        self.dir, self.meta, self.bucket = d / "buckets", meta, bucket
        self.lo = np.array(self.meta["lo"])
        self.cache, self.budget = OrderedDict(), budget

    def _get(self, i, j):
        k = (i, j)
        if k in self.cache:
            self.cache.move_to_end(k)
            return self.cache[k]
        p = self.dir / f"b_{i}_{j}.bin"
        a = np.fromfile(p, DT) if p.exists() else np.zeros(0, DT)
        self.cache[k] = a
        while sum(v.nbytes for v in self.cache.values()) > self.budget and len(self.cache) > 1:
            self.cache.popitem(last=False)
        return a

    def box(self, lo, hi):
        i0, j0 = np.floor((np.asarray(lo) - self.lo) / self.bucket).astype(int)
        i1, j1 = np.floor((np.asarray(hi) - self.lo) / self.bucket).astype(int)
        parts = [self._get(i, j) for i in range(max(i0, 0), i1 + 1) for j in range(max(j0, 0), j1 + 1)]
        a = np.concatenate(parts) if parts else np.zeros(0, DT)
        m = (a["x"] >= lo[0]) & (a["x"] <= hi[0]) & (a["y"] >= lo[1]) & (a["y"] <= hi[1])
        return a[m]


def boundary_counts(d, meta, cell):
    """Points per segment in the cells of the cloud's outer boundary: occupancy of the stored points on a cell
    grid, holes filled, boundary = filled cells next to empty ones; cells at the edge of the stored area (margin
    around the boundary, not the cloud's edge) excluded. Two streaming passes over the bucket files."""
    lo, hi = np.array(meta["lo"]), np.array(meta["hi"])
    shape = np.ceil((hi - lo) / cell).astype(int) + 1
    files = sorted((d / "buckets").glob("b_*.bin"))
    occ = np.zeros(shape, bool)
    for fp in files:
        a = np.fromfile(fp, DT)
        ij = np.floor((np.c_[a["x"], a["y"]] - lo) / cell).astype(int)
        occ[ij[:, 0], ij[:, 1]] = True
    fill = ndimage.binary_fill_holes(ndimage.binary_closing(occ, iterations=2))
    bnd = fill & ~ndimage.binary_erosion(fill, border_value=1)
    bnd[:3, :] = bnd[-3:, :] = bnd[:, :3] = bnd[:, -3:] = False
    bnd = ndimage.binary_dilation(bnd)
    cnt = {}
    for fp in files:
        a = np.fromfile(fp, DT)
        ij = np.floor((np.c_[a["x"], a["y"]] - lo) / cell).astype(int)
        i = a["id"][bnd[ij[:, 0], ij[:, 1]] & (a["id"] >= 0)]
        for k, n in zip(*np.unique(i, return_counts=True)):
            cnt[int(k)] = cnt.get(int(k), 0) + int(n)
    return cnt
