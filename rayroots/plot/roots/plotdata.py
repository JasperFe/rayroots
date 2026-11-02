# Frozen method (v1 + v2). Do not edit here: a change is a new version.
"""The one cloud of this run: segmented file, boundary hull, ground, segment ids."""
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

from rayroots.plot.roots import ground as rw

CACHE = None          # set by cli.use_tmp (the run's tmp folder)
CHUNK = 10_000_000
SEG = BOUNDARY = MESH = None   # set by setup


def setup(segmented, boundary=None, mesh=None):
    """segmented: the RCT file; boundary: csv path or None (whole cloud); mesh: path = ground from the RCT mesh."""
    global SEG, BOUNDARY, MESH
    SEG, BOUNDARY, MESH = segmented, boundary, mesh


def seg_path(plot):
    return SEG


def known_ids(plot):
    """Empty in this tool (no known trees). Kept so that cuts.py stays the frozen code."""
    return set()


def read_boundary(path):
    """(n, 2) float64 from a csv: columns x and y (any case), else the first two columns; comma or whitespace."""
    d = pd.read_csv(path, sep=r"[,\s]+", engine="python", header=None, dtype=str, comment="#")
    head = [str(v).strip().lower() for v in d.iloc[0]]
    if "x" in head and "y" in head:
        return d.iloc[1:, [head.index("x"), head.index("y")]].astype(np.float64).to_numpy()
    try:
        return d.iloc[:, :2].astype(np.float64).to_numpy()
    except ValueError:
        return d.iloc[1:, :2].astype(np.float64).to_numpy()   # another header


def hull_xy(plot):
    """Convex hull of the boundary points; without a boundary the bounding rectangle of the returned points."""
    if BOUNDARY is not None:
        p = read_boundary(BOUNDARY)
        return p[ConvexHull(p).vertices]
    cache = CACHE / "cloud_rectangle.npy"
    if cache.exists():
        return np.load(cache)
    lo, hi = np.full(2, np.inf), np.full(2, -np.inf)
    mm = rw.ply_memmap(SEG)
    for s in range(0, len(mm), CHUNK):
        c = np.asarray(mm[s:s + CHUNK][["x", "y", "alpha"]])
        c = c[c["alpha"] > 0]
        if len(c):
            lo = np.minimum(lo, [c["x"].min(), c["y"].min()])
            hi = np.maximum(hi, [c["x"].max(), c["y"].max()])
    rect = np.array([lo, [hi[0], lo[1]], hi, [lo[0], hi[1]]], np.float64)
    np.save(cache, rect)
    return rect


def ground(plot):
    return rw.Ground(SEG) if MESH is None else rw.MeshGround(MESH)


def memmap(path):
    """(memmap of all vertices, header bytes)."""
    with open(path, "rb") as f:
        head = b""
        while not head.endswith(b"end_header\n"):
            head += f.readline()
    return rw.ply_memmap(path), head


def chunks(path, fields=None):
    """All points (unbounded rays included) chunk by chunk: (start index, records, segment ids)."""
    mm, _ = memmap(path)
    for s in range(0, len(mm), CHUNK):
        c = np.asarray(mm[s:s + CHUNK] if fields is None else mm[s:s + CHUNK][fields])
        yield s, c, rw.decode_ids(c["red"], c["green"], c["blue"])


def encode_ids(ids):
    """Inverse of rw.decode_ids: segment id -> (r, g, b) as raycloudtools convertIntToColour; -1 -> black."""
    x = np.asarray(ids, np.int64) + 1
    rgb = [np.zeros(len(x), np.int64) for _ in range(3)]
    for i in range(24):
        rgb[i % 3] |= ((x >> i) & 1) << (7 - i // 3)
    return [c.astype(np.uint8) for c in rgb]


def hull_distance(p, hull):
    """Signed distance of points p (n, 2) to a convex polygon (negative inside)."""
    a, b = hull, np.roll(hull, -1, axis=0)
    ab = b - a
    out = np.empty(len(p))
    for s in range(0, len(p), 200_000):
        q = p[s:s + 200_000]
        t = np.clip(np.einsum("nkd,kd->nk", q[:, None] - a, ab) / (ab ** 2).sum(1), 0, 1)
        d = np.linalg.norm(q[:, None] - (a + t[..., None] * ab), axis=2).min(1)
        cross = ab[:, 0] * (q[:, None, 1] - a[:, 1]) - ab[:, 1] * (q[:, None, 0] - a[:, 0])
        inside = np.all(cross >= 0, axis=1) | np.all(cross <= 0, axis=1)
        out[s:s + 200_000] = np.where(inside, -d, d)
    return out


def segment_stats(plot, base_dz=0.3):
    """Per segment: n, z_min, z_max, base centroid (points within base_dz of z_min), returned points only. Cached."""
    path = seg_path(plot)
    cache = CACHE / f"{path.stem}_segstats_dz{base_dz}.csv"
    if cache.exists():
        return pd.read_csv(cache)
    K = 1 << 24
    n = np.zeros(K, np.int64); zmin = np.full(K, np.inf); zmax = np.full(K, -np.inf)
    for _, c, ids in chunks(path, ["z", "red", "green", "blue", "alpha"]):
        k = (ids >= 0) & (c["alpha"] > 0)
        i, z = ids[k], c["z"][k].astype(np.float64)
        n += np.bincount(i, minlength=K)
        np.minimum.at(zmin, i, z); np.maximum.at(zmax, i, z)
    sx = np.zeros(K); sy = np.zeros(K); nb = np.zeros(K, np.int64)
    for _, c, ids in chunks(path, ["x", "y", "z", "red", "green", "blue", "alpha"]):
        k = (ids >= 0) & (c["alpha"] > 0)
        i = ids[k]
        b = c["z"][k] <= zmin[i] + base_dz
        i = i[b]
        sx += np.bincount(i, c["x"][k][b].astype(np.float64), minlength=K)
        sy += np.bincount(i, c["y"][k][b].astype(np.float64), minlength=K)
        nb += np.bincount(i, minlength=K)
    u = np.flatnonzero(n)
    st = pd.DataFrame({"seg_id": u, "base_x": sx[u] / nb[u], "base_y": sy[u] / nb[u], "z_min": zmin[u],
                       "z_max": zmax[u], "n_points": n[u]})
    st.to_csv(cache, index=False)
    return st
