# Frozen method (v1 + v2). Do not edit here: a change is a new version.
"""PLY reading, RCT colour -> segment id, and the ground under the cloud: grid (frozen) or RCT terrain mesh."""
import numpy as np
from scipy import ndimage
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

GROUND_CELL, GROUND_K = 0.5, 20
CHUNK = 5_000_000
CACHE = None          # set by cli.use_tmp (the run's tmp folder)


def ply_memmap(path):
    types = {"float": "f4", "double": "f8", "uchar": "u1", "int": "i4", "uint": "u4",
             "short": "i2", "ushort": "u2", "char": "i1", "float32": "f4", "float64": "f8"}
    with open(path, "rb") as f:
        head = []
        while (s := f.readline().decode("ascii", "replace").strip()) != "end_header":
            head.append(s)
        off = f.tell()
    n = int(next(s for s in head if s.startswith("element vertex")).split()[2])
    props, in_vertex = [], False
    for s in head:
        if s.startswith("element"):
            in_vertex = s.split()[1] == "vertex"
        elif s.startswith("property") and in_vertex:
            props.append((s.split()[-1], "<" + types[s.split()[1]]))
    return np.memmap(path, dtype=np.dtype(props), mode="r", offset=off, shape=(n,))


def decode_ids(r, g, b):
    """raycloudtools convertColourToInt: bit-interleaved RGB -> index; black -> -1."""
    x = np.zeros(len(r), np.int64)
    chans = (np.asarray(r, np.int64), np.asarray(g, np.int64), np.asarray(b, np.int64))
    for i in range(24):
        x |= ((chans[i % 3] >> (7 - i // 3)) & 1) << i
    return x - 1


def chunks(path):
    """Returned points (alpha > 0) of a ray cloud, chunk by chunk, with their segment ids."""
    mm = ply_memmap(path)
    for s in range(0, len(mm), CHUNK):
        c = np.asarray(mm[s:s + CHUNK])
        c = c[c["alpha"] > 0]
        print(f"  {path.name}: {min(s + CHUNK, len(mm)):,} / {len(mm):,}", end="\r")
        yield c, decode_ids(c["red"], c["green"], c["blue"])
    print()


class Ground:
    """Ground height on a GROUND_CELL grid from the plot's own lowest points."""

    def __init__(self, path):
        p = CACHE / f"ground_{path.stem}.npz"
        if p.exists():
            z = np.load(p)
            self.x0, self.y0, self.g = float(z["x0"]), float(z["y0"]), z["g"]
            return
        OFF = 4000
        keep_c, keep_z = np.zeros(0, np.int64), np.zeros(0, np.float32)
        for c, _ in chunks(path):
            ix = np.floor(c["x"] / GROUND_CELL).astype(np.int64) + OFF
            iy = np.floor(c["y"] / GROUND_CELL).astype(np.int64) + OFF
            cc = np.concatenate([keep_c, ix * 2 * OFF + iy]); zz = np.concatenate([keep_z, c["z"]])
            o = np.lexsort((zz, cc)); cc, zz = cc[o], zz[o]
            first = np.r_[0, np.flatnonzero(np.diff(cc)) + 1]
            rank = np.arange(len(cc)) - np.repeat(first, np.diff(np.r_[first, len(cc)]))
            keep_c, keep_z = cc[rank < GROUND_K], zz[rank < GROUND_K]
        u, first, cnt = np.unique(keep_c, return_index=True, return_counts=True)
        gz = keep_z[first + cnt - 1]                       # K-th lowest (or highest kept if fewer)
        ix, iy = u // (2 * OFF) - OFF, u % (2 * OFF) - OFF
        self.x0, self.y0 = ix.min() * GROUND_CELL, iy.min() * GROUND_CELL
        g = np.full((ix.max() - ix.min() + 1, iy.max() - iy.min() + 1), np.nan, np.float32)
        g[ix - ix.min(), iy - iy.min()] = gz
        # cells far above their neighbourhood (no ground seen: stems, root thickets) -> neighbourhood median
        med = ndimage.generic_filter(g, np.nanmedian, size=5, mode="nearest")
        bad = np.isnan(g) | (g > med + 0.3)
        g[bad] = med[bad]
        if np.isnan(g).any():                              # fill gaps from the nearest cell
            idx = ndimage.distance_transform_edt(np.isnan(g), return_distances=False, return_indices=True)
            g = g[tuple(idx)]
        self.g = ndimage.median_filter(g, size=3, mode="nearest")
        np.savez(p, x0=self.x0, y0=self.y0, g=self.g)

    def __call__(self, x, y):
        """Bilinear ground height at x, y."""
        fi = (np.asarray(x, np.float64) - self.x0) / GROUND_CELL - 0.5
        fj = (np.asarray(y, np.float64) - self.y0) / GROUND_CELL - 0.5
        return ndimage.map_coordinates(self.g, [fi, fj], order=1, mode="nearest")


class MeshGround(Ground):
    """New, not frozen: the RCT terrain mesh sampled at the centres of a GROUND_CELL grid; lookup as Ground.
    Linear between the mesh vertices, nearest vertex outside them."""

    def __init__(self, mesh_path):
        p = CACHE / f"ground_mesh_{mesh_path.stem}.npz"
        if p.exists():
            z = np.load(p)
            self.x0, self.y0, self.g = float(z["x0"]), float(z["y0"]), z["g"]
            return
        v = ply_memmap(mesh_path)
        xy = np.c_[v["x"], v["y"]].astype(np.float64)
        z = np.asarray(v["z"], np.float64)
        lo = np.floor(xy.min(0) / GROUND_CELL) * GROUND_CELL
        hi = np.ceil(xy.max(0) / GROUND_CELL) * GROUND_CELL
        self.x0, self.y0 = float(lo[0]), float(lo[1])
        xs = np.arange(lo[0], hi[0], GROUND_CELL) + GROUND_CELL / 2
        ys = np.arange(lo[1], hi[1], GROUND_CELL) + GROUND_CELL / 2
        X, Y = np.meshgrid(xs, ys, indexing="ij")
        g = LinearNDInterpolator(xy, z)(X, Y)
        out = np.isnan(g)
        g[out] = NearestNDInterpolator(xy, z)(X[out], Y[out])
        self.g = g.astype(np.float32)
        np.savez(p, x0=self.x0, y0=self.y0, g=self.g)
