"""rayroots-single split: find the cut and split a finished tree into stems and roots (Rhizophora) or forks (Avicennia).
The core comes first: horizontal slices -> blobs per slice -> below the last prop root a slice has many blobs, above it
one per stem -> last root = bottom of the lowest clean run -> cut = last root + 0.6 m.
"""
# %% Imports
import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import ConvexHull, cKDTree
from scipy.spatial.distance import pdist

from rayroots.plot.finish import read_ply
from rayroots.single import figure, settings

# %% Settings of the noise filter and the ground
# 0.03 / 0.07 suits dense roots only; on whole trees it cuts off the sparser crown
noise_voxel = 0.05       # voxel size [m]
noise_link = 0.12        # max gap bridged between voxels [m]
z_ref_pct = 1      # percentile of z used as ground reference (ground "percentile")
flag_removed = 0.05  # flag a tree that loses more than this fraction of points
# ground "stem_base": where the stem meets the ground = the lowest point of the stem column, followed down slice by
# slice, so stray points and a sparse low tail do not decide it
col_raster = 0.03     # m, raster of the slices ; a slice's blobs continue the column if their cells
                      # touch the column's cells one slice up (one cell of slack)
col_slice = 0.05      # m, slice thickness and step going down
col_start = 0.5       # m above the 1st percentile of z, where the column is picked up (largest blob): the 1st
                      # percentile lies on the lower stem of every Avicennia tree (its lower stem holds few points)
sparse_h = 0.5        # m; flag "sparse base" when the column's lowest sparse_h m holds fewer than sparse_min points:
sparse_min = 1500     # its height is less certain. Set on the first 12 Avicennia trees: well-scanned bases hold
                      # 2800-7000 points there, SU8_PSP4_101 / _86 / SU7_psp3_123 ~900-1100 (a density relative to
                      # the stem higher up does not catch 101: its whole lower 1.5 m is sparse)

# %% Settings of the cut and the split
# cut detection
slice_thick = 0.05   # m
slice_step = 0.02    # m
raster = 0.03        # m
min_cells = 4        # blobs smaller than this are noise
persist = 0.40       # m of constant blob count needed above the last root
margin = 0.60        # m, cut = last root + margin (kept inside the clean run)
max_lean = 45        # degrees from vertical a stem may lean (flatter blobs are roots)
max_converge = 0.05  # m, stems may not move closer together over the run (roots converge on the stem)
scan_to = 10.0       # m above ground, highest the search looks
z0_pct = 0.1         # percentile of z used as the algorithm's ground (h_* are relative to this)

# root cleaning
voxel = 0.03         # m
link = 0.07          # m, max gap bridged

# stem separation (coarser: tops are big)
voxel_top = 0.05     # m
link_top = 0.10      # m

# a tree is flagged for checking when its clean stem run is shorter than this (m)
check_run = 0.60
# a "stem" with less than this fraction of the top's points is not a stem: merged into the nearest real
# stem if wood-connected to one (voxel_top/link_top graph), else it is a root tip poking above the cut -> roots
# (split "roots" only)
min_stem_frac = 0.01

# DBH split (split "dbh"), field rule: a fork below 1.3 m makes separate stems; above the DBH point a fork is a
# branch. Blobs are tracked bottom-up from base_h; slice_thick, slice_step, persist, max_lean, max_converge and the
# stem separation as above
min_dbh = 0.05                    # m; detect every stem the field inventory counts (smallest Avicennia stems ~5 cm)
raster_dbh = min_dbh / 3          # m; a min_dbh stem spans at least 3 cells
min_cells_dbh = int(np.pi * min_dbh / raster_dbh / 2)   # = 4: half the circumference of a min_dbh stem (one-sided scan)
dbh_h = 1.3                       # m above z_ref, field DBH height: forks up to here make stems
base_h = slice_thick / 2          # m above z_ref (the ground) of the first slice: tracking starts at the ground
fork_near = 0.15                  # m; a fork within this distance of dbh_h is flagged "fork near DBH"
late_start = 0.3                  # m; flag "late start" when the stem is picked up higher than this above the ground
                                  # (sparse base: a fork below the start height cannot be seen)


# %% Core: slices, blobs, the cut
def slice_profile(xyz, z0):
    """Per slice: height, number of blobs, blob centres (xy)."""
    order = np.argsort(xyz[:, 2])
    zs = xyz[order, 2]
    xy0 = xyz[:, :2].min(0)
    heights = np.arange(slice_thick / 2, scan_to, slice_step)
    counts, centres = [], []
    for h in heights:
        i0, i1 = np.searchsorted(zs, [z0 + h - slice_thick / 2, z0 + h + slice_thick / 2])
        q = xyz[order[i0:i1], :2]
        if len(q) == 0:
            counts.append(0); centres.append(np.empty((0, 2))); continue
        ij = np.floor((q - xy0) / raster).astype(int)
        grid = np.zeros(ij.max(0) + 1, bool)
        grid[ij[:, 0], ij[:, 1]] = True
        lab, n = ndimage.label(grid, np.ones((3, 3)))
        size = np.bincount(lab.ravel())[1:]
        big = np.where(size >= min_cells)[0] + 1
        c = np.array(ndimage.center_of_mass(grid, lab, big)).reshape(-1, 2)
        counts.append(len(big)); centres.append(c * raster + xy0 + raster / 2)
    return heights, np.array(counts), centres

def find_cut(heights, counts, centres):
    """Lowest run of `persist` m where the blob count stays N and every blob is near-vertical.
    Returns h_last_root, h_cut, n_stems, run length, stem centres at the cut."""
    need = int(round(persist / slice_step))
    for i in range(len(heights) - need):
        w = counts[i:i + need]
        n = w[0]
        if n == 0 or (w != n).any():
            continue
        # lean test: every blob at the bottom of the run has a blob close above it at the top
        bottom = centres[i]
        d, m = cKDTree(centres[i + need - 1]).query(bottom)
        if d.max() > persist * np.tan(np.radians(max_lean)):
            continue
        # convergence test: going up, roots move towards the stem, stems do not move towards each other
        top = centres[i + need - 1][m]
        gap_bottom = np.linalg.norm(bottom[:, None] - bottom[None], axis=2)
        gap_top = np.linalg.norm(top[:, None] - top[None], axis=2)
        if (gap_bottom - gap_top).max() > max_converge:
            continue
        j = i
        while j + 1 < len(counts) and counts[j + 1] == n:
            j += 1
        h_cut = min(heights[i] + margin, heights[j] - slice_thick)
        k = np.argmin(abs(heights - h_cut))
        return heights[i], h_cut, n, heights[j] - heights[i], centres[k]
    return np.nan, np.nan, 0, 0.0, None


# %% Reading and writing
ply_types = {"float": "<f4", "double": "<f8", "uchar": "u1", "char": "i1", "ushort": "<u2",
             "short": "<i2", "int": "<i4", "uint": "<u4"}
ply_names = {np.dtype(v).str: k for k, v in ply_types.items()}

def read_cloud(path):
    """Binary little-endian .ply as written by standardise."""
    raw = open(path, "rb").read()
    end = raw.index(b"end_header\n") + len(b"end_header\n")
    header = raw[:end].decode().splitlines()
    props = [(l.split()[2], ply_types[l.split()[1]]) for l in header if l.startswith("property")]
    n = int([l for l in header if l.startswith("element vertex")][0].split()[2])
    return np.frombuffer(raw, dtype=props, count=n, offset=end)

def write_ply(path, arr):
    lines = ["ply", "format binary_little_endian 1.0", f"element vertex {len(arr)}"]
    lines += [f"property {ply_names[arr.dtype[k].str]} {k}" for k in arr.dtype.names]
    lines += ["end_header"]
    with open(path, "wb") as f:
        f.write(("\n".join(lines) + "\n").encode())
        f.write(np.ascontiguousarray(arr).tobytes())


# %% Noise filter and ground
def largest_component(xyz, voxel, link):
    """Voxelise, link voxel centres closer than `link`, keep the largest component.
    Returns point mask, number of components and the point count of the second largest."""
    ijk = np.floor((xyz - xyz.min(0)) / voxel).astype(np.int64)
    _, point_vox = np.unique(ijk, axis=0, return_inverse=True)
    point_vox = point_vox.ravel()
    n_vox = point_vox.max() + 1

    centres = np.zeros((n_vox, 3))
    np.add.at(centres, point_vox, xyz)
    centres /= np.bincount(point_vox)[:, None]

    pairs = cKDTree(centres).query_pairs(link, output_type="ndarray")
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n_vox, n_vox))
    n_comp, comp = connected_components(graph, directed=False)

    comp_points = np.bincount(comp[point_vox], minlength=n_comp)  # points per component
    order = np.argsort(comp_points)[::-1]
    keep = comp[point_vox] == order[0]
    second = int(comp_points[order[1]]) if n_comp > 1 else 0
    return keep, n_comp, second

def slice_blobs(q, xy0):
    """Occupied col_raster cells of a slice's points, grouped into 8-connected blobs (list of cell sets)."""
    if not len(q):
        return []
    ij = np.unique(np.floor((q[:, :2] - xy0) / col_raster).astype(int), axis=0)
    off = ij.min(0)
    g = np.zeros(ij.max(0) - off + 1, bool)
    g[tuple((ij - off).T)] = True
    lab = ndimage.label(g, np.ones((3, 3)))[0][tuple((ij - off).T)]
    return [set(map(tuple, ij[lab == k])) for k in np.unique(lab)]


def stem_base(xyz):
    """Ground where the stem meets it: pick up the stem column (largest blob) col_start above the 1st percentile
    of z, follow it down in col_slice slices (blobs touching the column one slice up), and take the lowest point
    of the column. Returns ground z and the number of column points in its lowest sparse_h m."""
    xy0 = xyz[:, :2].min(0)
    order = np.argsort(xyz[:, 2]); zs = xyz[order, 2]
    sl = lambda z: xyz[order[np.searchsorted(zs, z):np.searchsorted(zs, z + col_slice)]]
    z = np.percentile(xyz[:, 2], 1) + col_start
    col = max(slice_blobs(sl(z), xy0), key=len)
    counts, z_low = [], z
    while True:
        z -= col_slice
        q = sl(z)
        grown = {(i + di, j + dj) for i, j in col for di in (-1, 0, 1) for dj in (-1, 0, 1)}
        nxt = [b for b in slice_blobs(q, xy0) if b & grown]
        if not nxt:
            break
        col = set().union(*nxt)
        cells = {c: None for c in col}
        ij = np.floor((q[:, :2] - xy0) / col_raster).astype(int)
        inside = np.array([tuple(c) in cells for c in ij])
        counts.append((z, int(inside.sum())))
        z_low = q[inside, 2].min()
    zc = np.array([c[0] for c in counts]); nc = np.array([c[1] for c in counts])
    return float(z_low), int(nc[zc < z_low + sparse_h].sum())

# %% Stems and roots
# %% Function: keep roots connected to the stem(s)
def connected_to_stem(xyz, z_top):
    ijk = np.floor((xyz - xyz.min(0)) / voxel).astype(int)
    _, pv = np.unique(ijk, axis=0, return_inverse=True)
    pv = pv.ravel()
    cen = np.zeros((pv.max() + 1, 3))
    np.add.at(cen, pv, xyz)
    cen /= np.bincount(pv)[:, None]
    pairs = cKDTree(cen).query_pairs(link, output_type="ndarray")
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(cen), len(cen)))
    _, piece = connected_components(g, directed=False)
    seed = cen[:, 2] > z_top - slice_thick
    return np.isin(piece, piece[seed])[pv]

# %% Function: assign every top point to a stem (shortest path through the wood)
def voxel_graph(xyz):
    """Wood graph for the stem assignment: voxel_top voxels, edges between voxel centres closer than link_top
    (weight = distance). Returns voxel per point, voxel centres, graph (csr)."""
    ijk = np.floor((xyz - xyz.min(0)) / voxel_top).astype(int)
    _, pv = np.unique(ijk, axis=0, return_inverse=True)
    pv = pv.ravel()
    cen = np.zeros((pv.max() + 1, 3))
    np.add.at(cen, pv, xyz)
    cen /= np.bincount(pv)[:, None]
    pairs = cKDTree(cen).query_pairs(link_top, output_type="ndarray")
    w = np.linalg.norm(cen[pairs[:, 0]] - cen[pairs[:, 1]], axis=1)
    g = coo_matrix((w, (pairs[:, 0], pairs[:, 1])), shape=(len(cen), len(cen))).tocsr()
    return pv, cen, g


def split_stems(xyz, z_cut, stem_xy):
    """Stem label per point: shortest path through the wood (voxel_graph) from seed voxels just above the cut,
    each seed given to the nearest stem centre. Voxels no seed reaches get the nearest labelled voxel's stem."""
    pv, cen, g = voxel_graph(xyz)

    # seeds: voxels just above the cut, each given to the nearest stem centre
    seeds = np.where(cen[:, 2] < z_cut + 0.10)[0]
    seed_stem = cKDTree(stem_xy).query(cen[seeds, :2])[1]
    _, _, src = dijkstra(g, directed=False, indices=seeds, min_only=True, return_predecessors=True)

    lab = np.full(len(cen), -1)
    ok = src >= 0
    lab[ok] = seed_stem[np.searchsorted(seeds, src[ok])]
    # voxels not connected to any stem (loose leaves, gaps): nearest labelled voxel
    if (~ok).any():
        lab[~ok] = lab[ok][cKDTree(cen[ok]).query(cen[~ok])[1]]
    return lab[pv]

# %% Function: tiny "stems": merge into a connected real stem, else send to the roots
def resolve_small_stems(xyz, stem, min_frac):
    """Stems with < min_frac of the points are fragments. A fragment wood-connected to a real stem
    (same component of the voxel_top/link_top graph) goes to the real stem with the closest point;
    an unconnected one gets label -1 (-> roots file).
    Returns labels (real stems 0..n-1, -1 = to roots), number merged, number sent to the roots."""
    n = stem.max() + 1
    frac = np.bincount(stem, minlength=n) / len(stem)
    real, small = np.where(frac >= min_frac)[0], np.where(frac < min_frac)[0]
    if len(small) == 0:
        return stem, 0, 0
    ijk = np.floor((xyz - xyz.min(0)) / voxel_top).astype(int)
    _, pv = np.unique(ijk, axis=0, return_inverse=True)
    pv = pv.ravel()
    cen = np.zeros((pv.max() + 1, 3))
    np.add.at(cen, pv, xyz)
    cen /= np.bincount(pv)[:, None]
    pairs = cKDTree(cen).query_pairs(link_top, output_type="ndarray")
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(cen), len(cen)))
    piece = connected_components(g, directed=False)[1][pv]  # component per point
    real_pieces = set(piece[np.isin(stem, real)])

    trees = {r: cKDTree(xyz[stem == r]) for r in real}
    new = stem.copy()
    n_merged = n_roots = 0
    for s in small:
        m = stem == s
        if real_pieces & set(piece[m]):
            d = [trees[r].query(xyz[m])[0].min() for r in real]
            new[m] = real[np.argmin(d)]
            n_merged += 1
        else:
            new[m] = -1
            n_roots += 1
    out = np.full(len(stem), -1)
    ok = new >= 0
    out[ok] = np.searchsorted(real, new[ok])
    return out, n_merged, n_roots

# %% Functions: DBH split
def blob_diameter(cells):
    """Largest extent of a blob (cell indices) in m: max distance between cell centres + one cell."""
    c = cells.astype(float)
    if len(c) > 3:
        try:
            c = c[ConvexHull(c).vertices]
        except Exception:      # collinear cells
            pass
    return ((pdist(c).max() if len(c) > 1 else 0) + 1) * raster_dbh


def dbh_slices(xyz, z_ref):
    """Blobs per horizontal slice from base_h to dbh_h + fork_near + persist above z_ref (forks are judged up to
    dbh_h + fork_near, each part needs persist above its fork), on one raster_dbh grid.
    Returns heights, per slice a list of blobs {centre xy, diam m, cells (set of grid indices)}, grid origin."""
    order = np.argsort(xyz[:, 2])
    zs = xyz[order, 2]
    xy0 = xyz[:, :2].min(0)
    heights = np.arange(base_h, dbh_h + fork_near + persist + 1e-9, slice_step)
    out = []
    for h in heights:
        i0, i1 = np.searchsorted(zs, [z_ref + h - slice_thick / 2, z_ref + h + slice_thick / 2])
        ij = np.unique(np.floor((xyz[order[i0:i1], :2] - xy0) / raster_dbh).astype(int), axis=0)
        blobs = []
        if len(ij):
            off = ij.min(0)
            g = np.zeros(ij.max(0) - off + 1, bool)
            g[tuple((ij - off).T)] = True
            lab = ndimage.label(g, np.ones((3, 3)))[0][tuple((ij - off).T)]
            for k in np.unique(lab):
                c = ij[lab == k]
                if len(c) >= min_cells_dbh:
                    blobs.append({"centre": (c.mean(0) + 0.5) * raster_dbh + xy0, "diam": blob_diameter(c),
                                  "cells": set(map(tuple, c))})
        out.append(blobs)
    return heights, out, xy0


def grow(cells):
    """Cells plus their 8 neighbours (one cell of slack between slices)."""
    return {(i + di, j + dj) for i, j in cells for di in (-1, 0, 1) for dj in (-1, 0, 1)}


def cells_centre(cells, xy0):
    return (np.array(list(cells)).mean(0) + 0.5) * raster_dbh + xy0


def follow_up(slices, k, cells, n):
    """Follow cells up from slice k for n slices (union of the overlapping blobs each slice).
    Returns the cells per slice (k .. k+n), or None if the chain ends."""
    chain = [cells]
    for kk in range(k + 1, min(k + n + 1, len(slices))):
        nxt = [b["cells"] for b in slices[kk] if b["cells"] & grow(chain[-1])]
        if not nxt:
            return None
        chain.append(set().union(*nxt))
    return chain if len(chain) == n + 1 else None


def valid_parts(slices, k, parts, xy0):
    """The parts of a dividing blob at slice k that are stems: each persists upward for `persist` m as its own
    blob (its chain never touches a sibling's), is >= min_dbh across (median over the chain) and leans <= max_lean;
    together they do not converge (as the Rhizophora stems). Returns the indices of the valid parts."""
    need = int(round(persist / slice_step))
    chains = {i: follow_up(slices, k, p["cells"], need) for i, p in enumerate(parts)}
    ok = []
    for i, ch in chains.items():
        if ch is None:
            continue
        others = [c for j, c in chains.items() if j != i and c is not None]
        if any(ch[t] & grow(o[t]) for o in others for t in range(len(ch))):
            continue   # rejoins a sibling: two sides of one stem, not a fork
        if np.median([blob_diameter(np.array(list(c))) for c in ch]) < min_dbh:
            continue
        if np.linalg.norm(cells_centre(ch[-1], xy0) - cells_centre(ch[0], xy0)) > persist * np.tan(np.radians(max_lean)):
            continue
        ok.append(i)
    if len(ok) >= 2:
        bottom = np.array([cells_centre(chains[i][0], xy0) for i in ok])
        top = np.array([cells_centre(chains[i][-1], xy0) for i in ok])
        gap_b = np.linalg.norm(bottom[:, None] - bottom[None], axis=2)
        gap_t = np.linalg.norm(top[:, None] - top[None], axis=2)
        if (gap_b - gap_t).max() > max_converge:
            return []
    return ok


def track_forks(heights, slices, xy0, max_stems=None):
    """Track blobs bottom-up from base_h. A track whose blob divides into >= 2 valid parts (valid_parts) at or
    below dbh_h forks: each part becomes its own track and is followed on (a second fork splits it again); above
    dbh_h (up to dbh_h + fork_near) a fork is a branch, only recorded. max_stems (n_stems override): forks are taken
    bottom-up only while the number of stems stays <= max_stems.
    Returns tracks [{cells per slice index, start (slice of its fork), parent, alive at dbh_h}], fork events."""
    k_top = int(np.argmin(abs(heights - dbh_h)))
    k_end = int(np.argmin(abs(heights - (dbh_h + fork_near))))
    # start: the lowest slice with blobs >= min_dbh that continue up to dbh_h (a sparse base has none at the ground);
    # blobs whose chains meet are one stem (clutter beside the stem rides up into it), stems already apart there stay
    # separate
    for k0 in range(k_top):
        start = [(b, follow_up(slices, k0, b["cells"], k_top - k0)) for b in slices[k0] if b["diam"] >= min_dbh]
        start = [(b, ch) for b, ch in start if ch is not None]
        if start:
            break
    group = list(range(len(start)))
    root = lambda i: i if group[i] == i else root(group[i])
    for i in range(len(start)):
        for j in range(i + 1, len(start)):
            if any(a & grow(c) for a, c in zip(start[i][1], start[j][1])):
                group[root(j)] = root(i)
    tracks = []
    for r in sorted({root(i) for i in range(len(start))}):
        cells = set().union(*(start[i][0]["cells"] for i in range(len(start)) if root(i) == r))
        tracks.append({"cells": {k0: cells}, "start": k0, "parent": None, "alive": True})
    forks = []
    for k in range(k0, k_end):
        n_alive = sum(t["alive"] for t in tracks)
        for t in [t for t in tracks if t["alive"] and k in t["cells"]]:
            parts = [b for b in slices[k + 1] if b["cells"] & grow(t["cells"][k])]
            if not parts:
                t["cells"][k + 1] = t["cells"][k]          # gap in the scan: carry on
                continue
            ok = valid_parts(slices, k + 1, parts, xy0) if len(parts) >= 2 else []
            branch = k + 1 > k_top
            if len(ok) < 2:
                t["divided"] = False
            elif t.get("divided"):
                ok = []                                  # same division as one slice lower, already recorded
            else:
                t["divided"] = True
                forks.append({"h": heights[k + 1], "n_parts": len(ok), "branch": branch,
                              "taken": not branch and (max_stems is None or n_alive + len(ok) - 1 <= max_stems)})
            if len(ok) >= 2 and forks[-1]["taken"]:
                t["alive"] = False
                t["fork"] = forks[-1]
                for i in ok:
                    tracks.append({"cells": {k + 1: parts[i]["cells"]}, "start": k + 1, "parent": t, "alive": True})
                n_alive += len(ok) - 1
            else:
                t["cells"][k + 1] = set().union(*(p["cells"] for p in parts))
        if k + 1 <= k_top:
            tracks = rejoin(tracks, k + 1)
    return tracks, forks


def rejoin(tracks, k):
    """Alive tracks whose cells touch at slice k are one stem: a fork whose parts rejoin is not a fork (the parent
    track lives on, or the two parts merge if the fork had more); other tracks that touch merge into the one that
    started lower (flagged via fork note "rejoined")."""
    alive = [t for t in tracks if t["alive"] and k in t["cells"]]
    for a in alive:
        for b in alive:
            if a is b or not (a["alive"] and b["alive"]) or not (a["cells"][k] & grow(b["cells"][k])):
                continue
            par = a["parent"]
            siblings = [t for t in tracks if t["parent"] is par and t["alive"]] if par is not None else []
            if par is not None and b["parent"] is par and len(siblings) == 2:
                for kk in set(a["cells"]) | set(b["cells"]):     # the fork was not a fork: the parent lives on
                    par["cells"][kk] = a["cells"].get(kk, set()) | b["cells"].get(kk, set())
                par["alive"] = True
                par["fork"]["taken"], par["fork"]["note"] = False, "rejoined"
                tracks = [t for t in tracks if t is not a and t is not b]
                return rejoin(tracks, k)
            keep, gone = (a, b) if a["start"] <= b["start"] else (b, a)
            for kk, c in gone["cells"].items():
                keep["cells"][kk] = keep["cells"].get(kk, set()) | c
            if gone["parent"] is not None:
                gone["parent"]["fork"]["note"] = "part rejoined"
                gone["parent"]["fork"]["n_parts"] -= 1
            tracks = [t for t in tracks if t is not gone]
            return rejoin(tracks, k)
    return tracks


def assign_dbh(xyz, z_ref, heights, tracks, xy0):
    """Stem label per point (-1 = base). Seeds: the points in each track's cells at each slice from its start up to
    dbh_h; stems = tracks alive at the end, forked (parent) tracks are base. The rest by the shortest path through
    the wood (voxel_graph); a stem's points below its own fork height and points no seed reaches go to the base."""
    stems = [t for t in tracks if t["alive"]]
    k_top = int(np.argmin(abs(heights - dbh_h)))
    lab_of = {id(t): (stems.index(t) if t in stems else -1) for t in tracks}
    pt_lab = np.full(len(xyz), -2)                          # -2 = not a seed
    k_pt = np.rint((xyz[:, 2] - z_ref - heights[0]) / slice_step).astype(int)
    ij = np.floor((xyz[:, :2] - xy0) / raster_dbh).astype(np.int64)
    key = ij[:, 0] * 1_000_000 + ij[:, 1]
    for t in tracks:
        for k, cells in t["cells"].items():
            if k > k_top:
                continue
            m = k_pt == k
            ck = np.array([i * 1_000_000 + j for i, j in cells])
            sel = np.where(m)[0][np.isin(key[m], ck)]
            pt_lab[sel] = lab_of[id(t)]
    pv, cen, g = voxel_graph(xyz)
    vox_lab = np.full(len(cen), -2)
    seeded = pt_lab > -2
    vox_lab[pv[seeded]] = pt_lab[seeded]
    seeds = np.where(vox_lab > -2)[0]
    _, _, src = dijkstra(g, directed=False, indices=seeds, min_only=True, return_predecessors=True)
    lab = np.full(len(cen), -1)
    ok = src >= 0
    lab[ok] = vox_lab[src[ok]]
    out = lab[pv]
    split_h = np.array([heights[t["start"]] for t in stems])
    low = (out >= 0) & (xyz[:, 2] < z_ref + split_h[np.maximum(out, 0)])
    out[low] = -1                                           # below its own fork: base
    return out, split_h


# %% Functions: one tree per split rule
def split_roots(tree_id, pc, xyz, h_override):
    """split "roots": stems above the cut, connected roots below. Returns outputs {name: points}, log fields,
    profile."""
    z0 = np.percentile(xyz[:, 2], z0_pct)
    heights, counts, centres = slice_profile(xyz, z0)
    h_root, h_cut, n_stems, run, stem_xy = find_cut(heights, counts, centres)
    flags = []
    if h_override is not None:
        h_cut = h_override
        k = np.argmin(abs(heights - h_cut))
        stem_xy, n_stems = centres[k], counts[k]
        flags.append("override")
    if n_stems == 0:
        raise RuntimeError("no clean stem window found, set h_cut_above_z0 in overrides.csv")
    z_cut = z0 + h_cut
    top = xyz[:, 2] >= z_cut

    # roots: keep what is connected to the stem(s)
    keep = connected_to_stem(xyz[~top], z_cut)

    # stems (always _stem1, _stem2, ...); unconnected fragments (label -1) join the roots
    if n_stems == 1:
        stem = np.zeros(top.sum(), int)
    else:
        stem = split_stems(xyz[top], z_cut, stem_xy)
    n_detected = n_stems
    stem, n_merged, n_frag = resolve_small_stems(xyz[top], stem, min_stem_frac)
    n_stems = n_detected - n_merged - n_frag
    if n_merged:
        flags.append("small stem merged")
    if n_frag:
        flags.append("fragment to roots")
    if run < check_run:
        flags.append("short run")
    out = {"roots": np.concatenate([pc[~top][keep], pc[top][stem < 0]])}
    out |= {f"stem{s + 1}": pc[top][stem == s] for s in range(n_stems)}
    row = {"z0": round(float(z0), 3), "z_cut": round(float(z_cut), 3),
           "h_cut_above_z0": round(float(h_cut), 3), "h_last_root": round(float(h_root), 3),
           "clean_run": round(float(run), 2), "n_stems_detected": int(n_detected), "n_stems": int(n_stems),
           "n_top": int((stem >= 0).sum()), "n_roots": int(keep.sum() + (stem < 0).sum()),
           "n_fragment_to_roots": int((stem < 0).sum()),
           "n_roots_removed": int((~keep).sum()), "roots_kept": round(float(keep.mean()), 4),
           "points_per_stem": "/".join(str((stem == s).sum()) for s in range(n_stems)), "flags": ";".join(flags)}
    profile = pd.DataFrame({"h_above_z0": heights, "n_blobs": counts})
    return out, row, profile


def split_dbh(tree_id, pc, xyz, z_ref, n_override):
    """split "dbh": forks below 1.3 m make stems (track_forks), points by assign_dbh, below the stems = base.
    One stem: no split, the stem is everything from z_ref up. n_override (manual_overrides.csv, n_stems): forks
    are taken bottom-up up to that many stems. Returns outputs {name: points}, log fields, profile."""
    heights, slices, xy0 = dbh_slices(xyz, z_ref)
    tracks, forks = track_forks(heights, slices, xy0, n_override)
    stems = [t for t in tracks if t["alive"]]
    n = len(stems)
    flags = []
    if n == 0:
        raise RuntimeError("no stem blob between the ground and 1.3 m that continues up to 1.3 m")
    if n_override is not None:
        flags.append("override n_stems")
        if n != n_override:
            raise RuntimeError(f"n_stems override {n_override}, but the fork tests give {n} stems (see the QC figure)")
    h_start = heights[min(t["start"] for t in tracks if t["parent"] is None)]
    if h_start > late_start:
        flags.append("late start")
    if any(abs(f["h"] - dbh_h) <= fork_near and f.get("note") != "rejoined" for f in forks):
        flags.append("fork near DBH")
    k_dbh = int(np.argmin(abs(heights - dbh_h)))
    d_dbh = [blob_diameter(np.array(list(t["cells"][k_dbh]))) if k_dbh in t["cells"] else np.nan for t in stems]
    if n == 1:
        z_cuts = [z_ref]
        keep = xyz[:, 2] >= z_ref
        out = {"stem1": pc[keep]}
        stem = np.zeros(keep.sum(), int)
        n_base, n_below = 0, int((~keep).sum())          # below z_ref: ground points, not modelled
    else:
        stem, split_h = assign_dbh(xyz, z_ref, heights, tracks, xy0)
        z_cuts = list(z_ref + split_h)
        above = xyz[:, 2] >= z_ref                      # below the ground: not modelled (as for one stem)
        out = {"base": pc[(stem < 0) & above]}
        out |= {f"stem{s + 1}": pc[stem == s] for s in range(n)}
        n_base, n_below = len(out["base"]), int((~above).sum())
    row = {"z_cut": round(float(min(z_cuts)), 3), "z_cut_per_stem": "/".join(f"{z:.3f}" for z in z_cuts),
           "h_split_per_stem": "/".join(f"{z - z_ref:.2f}" for z in z_cuts),
           "fork_heights": "/".join(f"{f['h']:.2f}({f['n_parts']}" + ("" if f["taken"] and "note" not in f else
                                    f", {f.get('note') or ('branch' if f['branch'] else 'not taken')}") + ")" for f in forks),
           "h_track_start": round(float(h_start), 2), "n_stems_detected": n, "n_stems": n,
           "D_extent_dbh_cm": "/".join(f"{100 * d:.1f}" for d in d_dbh),
           "n_top": int((stem >= 0).sum()), "n_base": n_base, "n_below_zref_dropped": n_below,
           "points_per_stem": "/".join(str((stem == s).sum()) for s in range(n)), "flags": ";".join(flags)}
    n_tracks = [sum(k in t["cells"] for t in tracks) for k in range(len(heights))]
    profile = pd.DataFrame({"h_above_zref": heights, "n_blobs": [len(sl) for sl in slices],
                            "n_stem_blobs": [sum(b["diam"] >= min_dbh for b in sl) for sl in slices],
                            "n_tracks": n_tracks})
    return out, row, profile

# %% Run
LOG = "split_log.csv"


def write_xyz(path, xyz):
    """Binary .ply, xyz as float64."""
    header = f"ply\nformat binary_little_endian 1.0\nelement vertex {len(xyz)}\n" \
             "property double x\nproperty double y\nproperty double z\nend_header\n"
    with open(path, "wb") as f:
        f.write(header.encode())
        f.write(np.ascontiguousarray(xyz, dtype="<f8").tobytes())


def standardise(path, ground):
    """Largest connected component, ground height z_ref. Returns cloud (structured, double), log fields."""
    a = read_ply(path)
    xyz = np.c_[a["x"], a["y"], a["z"]].astype(float)
    keep, n_comp, second = largest_component(xyz, noise_voxel, noise_link)
    kept = xyz[keep]
    removed = 1 - keep.mean()
    flags = ["removed>5%"] if removed > flag_removed else []
    if ground == "percentile":
        z_ref = round(float(np.percentile(kept[:, 2], z_ref_pct)), 3)
    else:
        z_ref, n_low = stem_base(kept)
        z_ref = round(z_ref, 3)
        if n_low < sparse_min:
            flags.append("sparse base")
    pc = np.empty(len(kept), [("x", "<f8"), ("y", "<f8"), ("z", "<f8")])
    pc["x"], pc["y"], pc["z"] = kept.T
    return pc, kept, {"n_in": len(xyz), "n_out": len(kept), "frac_removed": round(float(removed), 4), "z_ref": z_ref}, flags


def split_tree(tree, path, genus, rule, s, over):
    """One tree: standardise, split, write the parts and the QC figure. Returns a log row."""
    for sub in ("standard", "stems", "roots", "base", "profile"):
        (s / sub).mkdir(parents=True, exist_ok=True)
    row = {"tree": tree, "genus": genus, "split": rule["split"]}
    pc, xyz, fields, fl01 = standardise(path, rule["ground"])
    z_ref = fields["z_ref"]
    write_xyz(s / "standard" / f"{tree}.ply", xyz)
    if rule["split"] == "roots":
        out, f2, profile = split_roots(tree, pc, xyz, over.get("h_cut_above_z0"))
    else:
        n = over.get("n_stems")
        out, f2, profile = split_dbh(tree, pc, xyz, z_ref, None if n is None else int(n))
    for old in [*(s / "stems").glob(f"{tree}_stem*.ply"), s / "base" / f"{tree}_base.ply", s / "roots" / f"{tree}_roots.ply"]:
        old.unlink(missing_ok=True)
    for name, pts in out.items():
        folder = "stems" if name.startswith("stem") else name
        write_ply(s / folder / f"{tree}_{name}.ply", pts)
    profile.to_csv(s / "profile" / f"{tree}_profile.csv", index=False)
    flags = ";".join(x for x in fl01 + [f2.pop("flags")] if x)
    row |= {"status": "ok"} | fields | f2 | {"h_cut_above_zref": round(f2["z_cut"] - z_ref, 3), "flags": flags}
    figure.qc(s, row)
    return row


def run(project=".", tree=None, settings_path=None, redo=False, src=None, genus=None, out=None):
    """All finished trees of tree_register.csv, or the --input files; a tree that is done is skipped. --tree T redoes T,
    --redo all."""
    cfg = settings.load(settings_path)
    if src is None:
        part, todo = settings.part(project), settings.trees(project, cfg, tree)
    else:
        part, todo = settings.direct(cfg, src, genus, out, tree)
    s = part / "work"
    s.mkdir(parents=True, exist_ok=True)
    ov = settings.overrides(part)
    log_p = s / LOG
    log = pd.read_csv(log_p) if log_p.exists() else pd.DataFrame(columns=["tree", "status"])
    done = set(log.loc[log.status == "ok", "tree"])
    rows = []
    for t in todo:
        if t.tree in done and tree is None and not redo:
            continue
        print(t.tree, "...", flush=True)
        try:
            row = split_tree(t.tree, t.path, t.genus, cfg["genus"][t.genus], s, ov.get(t.tree, {}))
        except Exception as e:
            row = {"tree": t.tree, "genus": t.genus, "status": "error", "flags": repr(e)}
        rows.append(row)
        print(f"  {row['status']}: cut {row.get('h_cut_above_zref')} m above the ground, {row.get('n_stems')} stem(s) {row.get('flags', '')}", flush=True)
    if rows:
        new = pd.DataFrame(rows)
        log = pd.concat([log[~log.tree.isin(new.tree)], new], ignore_index=True).sort_values("tree")
        log.to_csv(log_p, index=False)
    print(f"{len(rows)} trees split, {len(todo) - len(rows)} already done; {log_p}")
    return 0
