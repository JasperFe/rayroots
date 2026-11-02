# Frozen method (v1 + v2). Do not edit here: a change is a new version.
"""Segments table, region points, stem tracks and the assignment region."""
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from rayroots.plot.roots import plotdata as cm

BUFFER = 2.0          # m around the scan-position hull: trees of the plot (as SU7_psp1)
R_MAX = 7.0           # m: largest root reach a fold can pick (root_window_check.py grid)
WINDOW = BUFFER + R_MAX   # m: segments up to here compete (their stems seed), whether in the plot or not
REGION = WINDOW + 3.0     # m: points kept around the hull (stem columns of the outermost competitors)
Z_MIN = 0.15          # m above terrain: ground removed (all methods)
HTOP = 8.6            # m above terrain: top of the kept points (S0 + slices)
S0 = 8.0              # m: highest start of a stem track (above any plausible root insertion)
H_START = 0.5         # m: lowest slice (skips pneumatophores)
SLICE, STEP, RASTER, MIN_CELLS = 0.05, 0.02, 0.03, 4
REACH_H = 1.0         # m: a root slice needs a non-stem tree blob that connects, slice by slice through non-stem tree blobs
                      # (1 cell dilation, one missing slice allowed), down to below this height: stilt roots reach
                      # the ground, crown clutter around the stem ends in the air
REACH_GAP = 1         # slices a chain may skip (occlusion)
ROOT_RUN = 5          # slices (0.1 m): a candidate root's chain must continue in each of these below it
CONV_DZ = 0.40        # m: convergence measured over this height 
CONV_MIN = 0.05       # m: the chain's largest blob CONV_DZ lower lies at least this much further from the stem
                      # 
REACH_DIL = 1         # cells of dilation between slices of a chain
GAP_STEM = 1          # cells of dilation for the stem track (a stem slice breaks into fragments on 3 cm raster)
STEM_MAX_W = 1.2      # m: a stem blob is narrower than this (root_workclouds.py STEM_MAX_W)
COL_LO = 1.5          # m: lowest layer of the stem-column search (above most stilt roots' lower parts)
RAD_GROW = 1.5        # a stem's radius stays below this x its median radius in the 1 m above
STEM_LOST = 0.5       # m: track stops after this many metres without an own blob near the centre
STEM_STEP = 0.06      # m: largest centre move per slice (0.02 m): a stem leans < ~70 deg
MASK_PAD = 0.05       # m added to another stem's radius when it is masked
MIN_HEIGHT = 2.0      # m: segments lower than this are not trees (no track; they can be absorbed)
OUT = None            # set by lib.use_tmp (the step's tmp folder)
HEIGHTS = np.round(np.arange(H_START, S0 + 1e-9, STEP), 4)
ASSIGN_REGION = BUFFER + R_MAX + 1.0   # m around the hull: points the assignment needs


def segments(plot):
    st = cm.segment_stats(plot)
    g = cm.ground(plot)
    st["base_hag"] = st.z_min - g(st.base_x.values, st.base_y.values)
    st["height"] = st.z_max - st.z_min
    st["dist_to_hull"] = cm.hull_distance(st[["base_x", "base_y"]].values, cm.hull_xy(plot))
    known = cm.known_ids(plot)
    st["known"] = st.seg_id.isin(known)
    st["in_plot"] = st.dist_to_hull <= BUFFER
    st["window"] = ((st.dist_to_hull <= WINDOW) | st.known) & (st.height >= MIN_HEIGHT)
    return st


def region(plot):
    path = cm.seg_path(plot)
    cache = OUT / f"{path.stem}_region.npz"
    if cache.exists():
        return dict(np.load(cache))
    g = cm.ground(plot)
    hull = cm.hull_xy(plot)
    st = segments(plot)
    known_xy = st.loc[st.known, ["base_x", "base_y"]].values
    lo, hi = hull.min(0) - REGION, hull.max(0) + REGION
    out = []
    for s, c, ids in cm.chunks(path, ["x", "y", "z", "red", "green", "blue", "alpha"]):
        x, y = c["x"].astype(np.float64), c["y"].astype(np.float64)
        m = (c["alpha"] > 0) & (x >= lo[0] - 15) & (x <= hi[0] + 15) & (y >= lo[1] - 15) & (y <= hi[1] + 15)
        hag = np.full(len(c), np.nan, np.float32)
        hag[m] = c["z"][m] - g(x[m], y[m])
        m &= (hag > Z_MIN) & (hag < HTOP)
        w = np.flatnonzero(m)
        keep = cm.hull_distance(np.c_[x[w], y[w]], hull) <= REGION
        if len(known_xy):        # known trees outside the window (SU7_psp3_199): their surroundings too
            keep |= cKDTree(known_xy).query(np.c_[x[w], y[w]], distance_upper_bound=R_MAX + 3)[0] <= R_MAX + 3
        w = w[keep]
        out.append(((w + s).astype(np.int32), np.c_[x[w], y[w], c["z"][w]].astype(np.float32), hag[w],
                    ids[w].astype(np.int32)))
        print(f"  region {plot}: {s + len(c):,}", end="\r", flush=True)
    print()
    idx, xyz, hag, ids = (np.concatenate(a) for a in zip(*out))
    r = {"idx": idx, "xyz": xyz, "hag": hag, "ids": ids}
    np.savez(cache, **r)
    return r


def by_segment(r):
    """Region point indices sorted by segment id, and {id: (start, stop)} into that order."""
    o = np.argsort(r["ids"], kind="stable")
    u, first, cnt = np.unique(r["ids"][o], return_index=True, return_counts=True)
    return o, dict(zip(u.tolist(), zip(first.tolist(), (first + cnt).tolist())))


class CellIndex:
    """Region points bucketed on a CELL m grid: fast 'points in a box'."""
    CELL = 2.0

    def __init__(self, xy):
        ij = np.floor(np.asarray(xy, np.float32) / self.CELL).astype(np.int32)    # plots span < +-10 km
        key = (ij[:, 0] + 5000) * 10000 + (ij[:, 1] + 5000)
        self.order = np.argsort(key, kind="stable").astype(np.int32)
        self.sorted = key[self.order]

    def box(self, lo, hi):
        i0, j0 = np.floor(np.asarray(lo) / self.CELL).astype(np.int64)
        i1, j1 = np.floor(np.asarray(hi) / self.CELL).astype(np.int64)
        out = []
        for i in range(i0, i1 + 1):
            a = (i + 5000) * 10000 + (j0 + 5000)
            b = (i + 5000) * 10000 + (j1 + 5000)
            out.append(self.order[np.searchsorted(self.sorted, a):np.searchsorted(self.sorted, b, "right")])
        return np.concatenate(out).astype(np.int64) if out else np.zeros(0, np.int64)


def blobs(q, origin, shape, gap):
    """Labelled raster of points q (n, 2) on a fixed frame; blobs smaller than MIN_CELLS original cells -> 0."""
    ij = np.floor((q - origin) / RASTER).astype(int)
    ok = (ij >= 0).all(1) & (ij[:, 0] < shape[0]) & (ij[:, 1] < shape[1])
    grid = np.zeros(shape, bool)
    grid[ij[ok, 0], ij[ok, 1]] = True
    dil = ndimage.binary_dilation(grid, np.ones((3, 3)), iterations=gap) if gap else grid
    lab, n = ndimage.label(dil, np.ones((3, 3)))
    if n:
        size = np.bincount(lab[grid], minlength=n + 1)
        small = np.flatnonzero(size < MIN_CELLS)
        lab[np.isin(lab, small)] = 0
        lab[~grid] = 0                        # blobs are original cells; dilation only joined them
    return lab


def stem_column(xy, hag):
    """Centre of the segment's stem: the 0.3 m column (3 x 3 cells of 0.1 m) of its own points between COL_LO and
    S0 that covers the most 0.5 m layers, most points on a tie (the trunk is the tallest dense column; a crown
    branch or a thin pole is not)."""
    m = (hag >= COL_LO) & (hag <= S0)
    if m.sum() < 50:
        return None
    q, h = xy[m], hag[m]
    lo = q.min(0)
    ij = np.floor((q - lo) / 0.1).astype(int)
    shape = tuple(ij.max(0) + 1)
    if shape[0] * shape[1] > 4_000_000:
        return None
    layer = np.floor((h - COL_LO) / 0.5).astype(int)
    nl = layer.max() + 1
    cov = np.zeros(shape, int); npts = np.zeros(shape)
    for k in range(nl):
        g = np.zeros(shape); np.add.at(g, (ij[layer == k, 0], ij[layer == k, 1]), 1)
        g3 = ndimage.uniform_filter(g, 3, mode="constant") * 9
        cov += g3 >= 5; npts += g3
    best = np.argwhere(cov == cov.max())
    b = best[np.argmax(npts[tuple(best.T)])]
    c = lo + (b + 0.5) * 0.1
    near = np.linalg.norm(q - c, axis=1) <= 0.3
    return q[near].mean(0) if near.any() else c


def track(xy, hag):
    """Stem centre (x, y), radius per HEIGHTS slice from the segment's own points; NaN where not tracked.
    Starts at the stem column (stem_column) at its highest slice <= S0 and follows it down: per slice the blob
    cells near the previous centre; centre moves <= STEM_STEP per slice; radius <= RAD_GROW x the median radius of
    the slices above (a stem does not widen much; roots merging with it would)."""
    cen = np.full((len(HEIGHTS), 2), np.nan); rad = np.full(len(HEIGHTS), np.nan)
    c_start = stem_column(xy, hag)
    if c_start is None:
        return cen, rad
    o = np.argsort(hag); hs, q = hag[o], xy[o]
    near0 = np.linalg.norm(q - c_start, axis=1) <= 0.3
    top = hs[near0 & (hs <= S0)]
    if not len(top):
        return cen, rad
    k_top = int(np.searchsorted(HEIGHTS, top.max()) - 1)
    origin = xy.min(0) - 2 * RASTER
    shape = tuple(np.floor((xy.max(0) - origin) / RASTER).astype(int) + 3)
    if shape[0] * shape[1] > 6_000_000:
        return cen, rad
    prev, lost = c_start, 0.0
    for k in range(max(k_top, 0), -1, -1):
        h = HEIGHTS[k]
        i0, i1 = np.searchsorted(hs, [h - SLICE / 2, h + SLICE / 2])
        lab = blobs(q[i0:i1], origin, shape, GAP_STEM) if i1 > i0 else None
        above = rad[k + 1:k + 1 + 50]
        r_med = np.nanmedian(above) if np.isfinite(above).any() else 0.15
        r0 = rad[k + 1] if k + 1 < len(rad) and np.isfinite(rad[k + 1]) else r_med
        p, cells = None, None
        if lab is not None and lab.any():
            c = np.argwhere(lab > 0)
            pc = (c + 0.5) * RASTER + origin
            d = np.linalg.norm(pc - prev, axis=1)
            near = d <= r0 * 1.5 + STEM_STEP
            if near.any():
                b = np.unique(lab[tuple(c[near].T)])
                cells = pc[np.isin(lab[tuple(c.T)], b) & (d <= RAD_GROW * r_med + 0.1)]
                if len(cells):
                    p = cells.mean(0)
                    mv = np.linalg.norm(p - prev)
                    if mv > STEM_STEP:
                        p = prev + (p - prev) * STEM_STEP / mv
        if p is None:
            lost += STEP
            if lost > STEM_LOST:
                break
            cen[k] = prev; rad[k] = r0
            continue
        lost = 0.0
        d = np.linalg.norm(cells - p, axis=1)
        cen[k] = p
        rad[k] = float(np.clip(np.percentile(d, 90), 0.03, min(STEM_MAX_W / 2, RAD_GROW * r_med)))
        prev = p
    return cen, rad


def stems(plot, r=None):
    p = OUT / f"{plot}_stems.npz"
    if p.exists():
        z = np.load(p)
        return {int(k[1:]): (z[k][:, :2], z[k][:, 2]) for k in z.files}
    r = r if r is not None else region(plot)
    st = segments(plot)
    o, span = by_segment(r)
    out = {}
    for sid in st.loc[st.window, "seg_id"]:
        a, b = span.get(int(sid), (0, 0))
        q = o[a:b]
        cen, rad = track(r["xyz"][q, :2].astype(np.float64), r["hag"][q].astype(np.float64))
        if np.isfinite(rad).any():
            out[int(sid)] = (cen, rad)
    np.savez(p, **{f"s{k}": np.c_[c, rr] for k, (c, rr) in out.items()})
    return out



def assign_region(plot):
    """Region points within ASSIGN_REGION of the hull (cached), as float32."""
    p = cm.CACHE / f"{plot}_region_assign.npz"
    if p.exists():
        return dict(np.load(p))
    r = region(plot)
    keep = cm.hull_distance(r["xyz"][:, :2].astype(np.float64), cm.hull_xy(plot)) <= ASSIGN_REGION
    r = {k: v[keep] for k, v in r.items()}
    np.savez(p, **r)
    return r

