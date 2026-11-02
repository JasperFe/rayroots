# Frozen. Do not edit here: a change is a new version.
"""Shared RayCloudTools helpers for the stems: run rayimport / rayextract trees, read tree files, flat terrain,
and stem diameters from the points. In a script: `import rct_tools as rt`."""
import shutil
import subprocess
from pathlib import Path
import numpy as np
from scipy import ndimage
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

# ensemble (steps 03 and 04 share it, so stem and roots of a tree get the same cut shift and seed per run).
# RCT is deterministic and its parameters spread volume by < 2 %, so only the input varies:
ENS_FRAC = 0.9                     # fraction of points kept per subsample (which points happened to be scanned)
ENS_SEEDS = [0, 1, 2, 3, 4]        # fixed seeds: the ensemble is reproducible
ENS_SHIFTS = [-0.10, 0.0, 0.10]    # m, cut height shifts; the cut is uncertain by about this much
ENS_RUNS = [(shift, seed) for shift in ENS_SHIFTS for seed in ENS_SEEDS]


def ens_name(shift, seed):
    """Run name shared by stem and roots, e.g. cut+0.10_seed3."""
    return f"cut{shift:+.2f}_seed{seed}"


# reference diameters from the points (defaults used by the analysis tests and 03b)
SLAB = 0.10          # m, slab thickness along the stem axis for the circle fit
AXIS_HALF = 0.3      # m, stem axis fitted through tracked slice centres within this distance
FORK_FRAC = 0.3      # a second blob near the stem larger than this fraction of the stem = junction
JUNCTION_GAP = 0.15  # m, no reference this close to a junction (stem and branch merged there)
LOST_FRAC = 0.5      # after a junction, a followed blob below this fraction of the stem area is a branch: stop
MAX_RMSE = 0.015     # m, circle fits worse than this give no reference
MAX_LEAN_REF = 45    # degrees; an axis flatter than this follows a branch, not the stem (as max_lean in 02)


# ---------------------------------------------------------------- io
def read_xyz(path):
    """xyz (float64, n x 3) of a binary little-endian .ply."""
    raw = open(path, "rb").read()
    end = raw.index(b"end_header\n") + len(b"end_header\n")
    header = raw[:end].decode().splitlines()
    types = {"float": "<f4", "double": "<f8", "uchar": "u1", "int": "<i4", "uint": "<u4", "ushort": "<u2"}
    props, in_vertex = [], False
    for l in header:
        if l.startswith("element"):
            in_vertex = l.split()[1] == "vertex"
        elif l.startswith("property") and in_vertex:
            props.append((l.split()[2], types[l.split()[1]]))
    n = int([l for l in header if l.startswith("element vertex")][0].split()[2])
    a = np.frombuffer(raw, dtype=props, count=n, offset=end)
    return np.c_[a["x"], a["y"], a["z"]].astype(float)


def write_xyz(path, xyz):
    header = (f"ply\nformat binary_little_endian 1.0\nelement vertex {len(xyz)}\n"
              "property double x\nproperty double y\nproperty double z\nend_header\n")
    with open(path, "wb") as f:
        f.write(header.encode())
        f.write(np.ascontiguousarray(xyz, dtype="<f8").tobytes())


GROUND_OFFSET = 0.001   # m; the flat ground goes this far above the cloud's lowest point (see write_flat_ground)


def write_flat_ground(path, xyz, step=0.1, margin=0.5, z_floor=None):
    """Flat triangle grid over the cloud's xy extent (+ margin), GROUND_OFFSET above its lowest point (or above
    z_floor, e.g. z_ref for a cloud with ground points below it): the cut plane of a cut-off top, or the ground. Tree bases snap to mesh vertices (raysegment.cpp), so the spacing bounds how
    far a base can sit from the stem. Why the offset: RCT's height field (raymesh.cpp toHeightField) only finds
    triangles inside the cloud's box, whose floor is the lowest point; a plane at or below it misses every cell
    centre and RCT takes a fallback that reads past a triangle's corners (undefined memory; warning "mesh does
    not intersect any pixel centres")."""
    xy_lo, xy_hi = xyz[:, :2].min(0) - margin, xyz[:, :2].max(0) + margin
    z = max(xyz[:, 2].min(), z_floor if z_floor is not None else -np.inf) + GROUND_OFFSET
    xs, ys = np.arange(xy_lo[0], xy_hi[0] + step, step), np.arange(xy_lo[1], xy_hi[1] + step, step)
    X, Y = np.meshgrid(xs, ys)
    v = np.c_[X.ravel(), Y.ravel(), np.full(X.size, z)].astype("<f4")
    nx, ny = len(xs), len(ys)
    j, i = np.meshgrid(np.arange(ny - 1), np.arange(nx - 1), indexing="ij")
    a = (j * nx + i).ravel(); b = a + 1; c = a + nx + 1; d = a + nx
    faces = np.c_[np.full(2 * len(a), 3), np.r_[a, a], np.r_[b, c], np.r_[c, d]].astype("<i4")
    header = (f"ply\nformat binary_little_endian 1.0\nelement vertex {len(v)}\nproperty float x\nproperty float y\n"
              f"property float z\nelement face {len(faces)}\nproperty list int int vertex_indices\nend_header\n")
    with open(path, "wb") as f:
        f.write(header.encode()); f.write(v.tobytes()); f.write(faces.tobytes())


def read_trees(path):
    """rayextract trees.txt: one line per tree, segments 'x,y,z,radius,parent_id,section_id'.
    Returns a list of dicts (xyz, r, par, volume = sum of cylinder volumes), largest first."""
    out = []
    for line in open(path):
        if line.startswith("#") or line.startswith("x,"):
            continue
        a = np.array([[float(v) for v in s.strip(",").split(",")] for s in line.split()])
        if a.ndim != 2 or not len(a):
            continue
        xyz, r, par = a[:, :3], a[:, 3], a[:, 4].astype(int)
        vol = sum(np.pi * r[i] ** 2 * np.linalg.norm(xyz[i] - xyz[par[i]]) for i in range(len(a)) if par[i] >= 0)
        out.append({"xyz": xyz, "r": r, "par": par, "volume": vol})
    return sorted(out, key=lambda t: -t["volume"])


# ---------------------------------------------------------------- RCT runs
def rct(cmd, cwd=None, timeout=3600):
    """Run an RCT tool; returns '' on success, else the last line of its output."""
    try:
        r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=timeout, cwd=cwd)
    except subprocess.TimeoutExpired:
        return f"timeout {timeout} s"
    if r.returncode != 0:
        out = (r.stderr.strip() or r.stdout.strip()).splitlines()
        return f"exit {r.returncode}: {out[-1] if out else ''}"
    return ""


def rayimport(ply, raycloud, ray_vector="0,0,-0.1"):
    """rayimport writes <name>_raycloud.ply next to its input; move it to `raycloud`. Returns error or ''."""
    made = Path(ply).with_name(Path(ply).stem + "_raycloud.ply")
    made.unlink(missing_ok=True)
    err = rct(["rayimport", ply, "ray", ray_vector, "--max_intensity", "0"])
    if err or not made.exists():
        return err or "rayimport: no output"
    shutil.move(made, raycloud)
    return ""


def rayextract_trees(raycloud, mesh, settings, out_prefix, keep=("trees.txt",)):
    """rayextract trees on `raycloud` with ground `mesh`; outputs (written next to the ray cloud) are
    renamed to <out_prefix>_<suffix> for the suffixes in `keep` and the rest deleted.
    Returns error or ''."""
    raycloud = Path(raycloud)
    stem = raycloud.stem
    made = {s: raycloud.with_name(f"{stem}_{s}") for s in ("trees.txt", "trees_mesh.ply", "segmented.ply")}
    for p in made.values():
        p.unlink(missing_ok=True)
    args = [a for k, v in settings.items() for a in (f"--{k}", str(v))]
    err = rct(["rayextract", "trees", raycloud, mesh, *args])
    ok = made["trees.txt"].exists()
    for s, p in made.items():
        if p.exists():
            if s in keep and ok:
                p.rename(f"{out_prefix}_{s}")
            else:
                p.unlink()
    return "" if ok else (err or "rayextract trees: no trees file")


# ---------------------------------------------------------------- diameters
def rct_diameter(tree, z):
    """Follow the trunk (thickest child) up from the base; diameter interpolated at height z."""
    kids = {}
    for i, p in enumerate(tree["par"]):
        if p >= 0:
            kids.setdefault(p, []).append(i)
    i = 0
    while i in kids:
        j = max(kids[i], key=lambda k: tree["r"][k])
        z0, z1 = tree["xyz"][i, 2], tree["xyz"][j, 2]
        if z1 >= z:
            f = 0 if z1 == z0 else np.clip((z - z0) / (z1 - z0), 0, 1)
            return 2 * (tree["r"][i] + f * (tree["r"][j] - tree["r"][i]))
        i = j
    return np.nan


def blobs(q, cell=0.02, min_cells=4):
    """2D connected blobs of a slice: list of (centre, area m2, member mask)."""
    lo = q.min(0)
    ij = np.floor((q - lo) / cell).astype(int)
    g = np.zeros(ij.max(0) + 1, bool); g[ij[:, 0], ij[:, 1]] = True
    lab, n = ndimage.label(g, np.ones((3, 3)))
    size = np.bincount(lab.ravel())
    pl = lab[ij[:, 0], ij[:, 1]]
    return [(q[pl == k].mean(0), size[k] * cell ** 2, pl == k) for k in range(1, n + 1) if size[k] >= min_cells]


def fit_circle(p):
    """Geometric least-squares circle (robust soft_l1) from an algebraic start. Returns centre, r, rmse."""
    A = np.c_[2 * p, np.ones(len(p))]
    c0 = np.linalg.lstsq(A, (p ** 2).sum(1), rcond=None)[0]
    r0 = np.sqrt(c0[2] + c0[0] ** 2 + c0[1] ** 2)
    res = least_squares(lambda c: np.hypot(*(p - c[:2]).T) - c[2], [c0[0], c0[1], r0], loss="soft_l1", f_scale=0.01)
    return res.x[:2], abs(res.x[2]), float(np.sqrt(np.mean(res.fun ** 2)))


def track_stem(pts, z_cut, top=4.6, step=0.05, start_xy=None):
    """Follow the main stem up from the cut: per slice the largest blob within 0.4 m of the previous
    centre (first slice: of `start_xy` if given, else any blob). A second nearby blob larger than FORK_FRAC of the
    stem marks a junction; tracking continues on the larger blob unless it has less than LOST_FRAC of the stem's
    area (a fork: the stem is lost).
    Returns {height above cut: centre xyz}, junction heights, height where the stem was lost (or None)."""
    centres, junctions, prev, areas = {}, [], None if start_xy is None else np.asarray(start_xy, float), []
    for h in np.arange(0.05, top, step):
        q = pts[np.abs(pts[:, 2] - z_cut - h) < 0.025]
        bl = blobs(q[:, :2]) if len(q) else []
        near = bl if prev is None else [b for b in bl if np.linalg.norm(b[0] - prev) < 0.4]
        if not near:
            continue
        near.sort(key=lambda b: -b[1])
        if len(near) > 1 and prev is not None and near[1][1] > FORK_FRAC * near[0][1]:
            junctions.append(round(h, 2))
        stem_area = np.median(areas[-10:]) if areas else near[0][1]
        if junctions and near[0][1] < LOST_FRAC * stem_area:
            return centres, junctions, round(h, 2)
        if not junctions or junctions[-1] < round(h, 2) - JUNCTION_GAP:
            areas.append(near[0][1])   # stem area, not taken where stem and branch merge
        centres[round(h, 2)] = np.r_[near[0][0], z_cut + h]
        prev = near[0][0]
    return centres, junctions, None


def point_diameter(pts, z_cut, h, centres):
    """Circle-fit diameter at h above the cut, in the plane perpendicular to the stem axis fitted through
    the tracked centres within +-AXIS_HALF m (a horizontal slice through a leaning stem is an ellipse).
    Returns diameter, fit rmse, axis lean (degrees), number of points."""
    cs = np.array([c for hh, c in centres.items() if abs(hh - h) <= AXIS_HALF])
    if len(cs) < 4:
        return np.nan, np.nan, np.nan, 0
    m = cs.mean(0)
    axis = np.linalg.svd(cs - m)[2][0]
    axis = axis if axis[2] > 0 else -axis
    c = m + axis * (z_cut + h - m[2]) / axis[2]          # axis point at height h
    near = pts[np.linalg.norm(np.cross(pts - c, axis), axis=1) < 0.4]
    s = near[np.abs((near - c) @ axis) < SLAB / 2]
    u = np.cross(axis, [1, 0, 0]); u /= np.linalg.norm(u); v = np.cross(axis, u)
    p2 = np.c_[(s - c) @ u, (s - c) @ v]
    bl = blobs(p2) if len(p2) else []
    if not bl:
        return np.nan, np.nan, np.nan, 0
    p2 = p2[min(bl, key=lambda b: np.linalg.norm(b[0]))[2]]   # blob on the axis
    _, r, rmse = fit_circle(p2)
    return 2 * r, rmse, np.degrees(np.arccos(axis[2])), len(p2)


def reference_diameters(pts, z_cut, heights):
    """Stem diameter from the points at each height above the cut, or NaN with a note why not.
    Returns (list of row dicts, stem centre xy at the cut)."""
    centres, junctions, lost = track_stem(pts, z_cut)
    rows = []
    for h in heights:
        row = {"h": h, "D_ref_cm": np.nan, "note": ""}
        close = [j for j in junctions if abs(j - h) <= JUNCTION_GAP]
        if lost is not None and h + AXIS_HALF > lost:
            row["note"] = f"stem lost at {lost} m (fork)"
        elif close:
            row["note"] = f"junction at {close[0]} m"
        else:
            d, rmse, lean, n = point_diameter(pts, z_cut, h, centres)
            note = ("no fit" if not np.isfinite(d) else "fit rmse too high" if rmse > MAX_RMSE
                    else "axis leans > 45 deg (branch)" if lean > MAX_LEAN_REF else "")
            row |= {"D_ref_cm": round(d * 100, 1) if not note else np.nan, "D_fit_cm": round(d * 100, 1),
                    "fit_rmse_mm": round(rmse * 1000, 1), "lean_deg": round(lean, 1), "n_pts": n, "note": note}
        rows.append(row)
    return rows, centres[min(centres)][:2]


def clean_diameter(pts, z0, h_target, h_min, start_xy=None, window=0.30, max_rel_rmse=None, step=0.05):
    """Field rule for DBH: the point diameter at h_target above z0, or, where a junction lies within JUNCTION_GAP
    (or the stem is lost / not fitted there), at the clean height nearest h_target within +-window (below first on
    a tie: below a swelling, else above), never below h_min and never outside the window.
    Fit note when rmse > max_rel_rmse x the fitted radius (None: rmse > MAX_RMSE).
    Returns a dict (h, D_fit_cm, fit_rmse_mm, lean_deg, n_pts, xy: stem centre at h, note: fit quality) or None."""
    centres, junctions, lost = track_stem(pts, z0, top=h_target + window + AXIS_HALF + step, start_xy=start_xy)
    if not centres:
        return None
    n = int(round(window / step))
    for k in sorted(range(-n, n + 1), key=lambda k: (abs(k), k > 0)):
        h = round(h_target + k * step, 2)
        if h < h_min - 1e-9:
            continue
        if (lost is not None and h + AXIS_HALF > lost) or any(abs(j - h) <= JUNCTION_GAP for j in junctions):
            continue
        d, rmse, lean, n_pts = point_diameter(pts, z0, h, centres)
        if not np.isfinite(d):
            continue
        bad_fit = rmse > (MAX_RMSE if max_rel_rmse is None else max_rel_rmse * d / 2)
        note = "fit rmse too high" if bad_fit else "axis leans > 45 deg (branch)" if lean > MAX_LEAN_REF else ""
        xy = centres[min(centres, key=lambda c: abs(c - h))][:2]
        return {"h": float(h), "D_fit_cm": round(d * 100, 1), "fit_rmse_mm": round(rmse * 1000, 1),
                "lean_deg": round(lean, 1), "n_pts": n_pts, "xy": xy, "note": note}
    return None


def diameter_at(trees, radii, z, xy, max_dist=0.3):
    """Diameter at height z of the segment (any tree) crossing z nearest to the stem centre xy (within max_dist m),
    with the given radii per tree (RCT's: [t['r'] for t in trees]). NaN if none."""
    best, out = max_dist, np.nan
    for t, r in zip(trees, radii):
        i = np.where(t["par"] >= 0)[0]
        p = t["par"][i]
        z0, z1 = t["xyz"][p, 2], t["xyz"][i, 2]
        m = (np.minimum(z0, z1) <= z) & (np.maximum(z0, z1) >= z) & (z0 != z1)
        if not m.any():
            continue
        i, p = i[m], p[m]
        f = ((z - t["xyz"][p, 2]) / (t["xyz"][i, 2] - t["xyz"][p, 2]))[:, None]
        c = t["xyz"][p, :2] + f * (t["xyz"][i, :2] - t["xyz"][p, :2])
        d = np.linalg.norm(c - xy, axis=1)
        k = d.argmin()
        if d[k] < best:
            best, out = d[k], 2 * (r[p[k]] + f[k, 0] * (r[i[k]] - r[p[k]]))
    return float(out)


FOOT_H = 0.3   # m above the plane over which the stem radius at the plane is taken (rct_tools.stem_at_cut)


def stem_at_cut(pts, z_cut):
    """Stem centre (xy, first slice track_stem follows) and radius (m) at the plane (the cut, or the ground for an
    Avicennia single stem). Radius = half the largest horizontal extent of the followed blob in the slices up to
    FOOT_H above the plane: the plane lies at the lowest point, so on a sloping or sparse stem foot the first slices
    only cut its downhill edge (SU7_psp3_123: 0.03 m at 0.05 m, 0.12 m at 0.3 m); a one-sided scan still spans the
    stem's width, and a circle fit fails on a flared, cluttered foot."""
    centres, _, _ = track_stem(pts, z_cut, top=FOOT_H + 0.05)
    if not centres:
        return np.full(2, np.nan), np.nan
    r = 0.0
    for h, c in centres.items():
        q = pts[np.abs(pts[:, 2] - z_cut - h) < 0.025][:, :2]
        blob = min(blobs(q), key=lambda b: np.linalg.norm(b[0] - c[:2]))
        p = q[blob[2]] - q[blob[2]].mean(0)
        axis = np.linalg.svd(p, full_matrices=False)[2][0]   # longest horizontal direction of the blob
        r = max(r, float(np.ptp(p @ axis) / 2))
    return centres[min(centres)][:2], r


def base_check(tree, pts, z_cut, margin):
    """Distance of an RCT tree's base to the stem centre at the plane, the stem radius there, and whether the base
    lies within radius + margin. On a ground plane the base snaps to the edge of the stem foot (the lowest point),
    on a horizontal cut (Rhizophora) near the centre: radius + margin holds for both."""
    c, r = stem_at_cut(pts, z_cut)
    d = float(np.linalg.norm(tree["xyz"][0, :2] - c))
    return d, r, bool(d <= r + margin)


# ---------------------------------------------------------------- input variation
def shift_cut(stem_xyz, tree_xyz, stems_xyz, idx, z_cut, shift, band_link=0.3, sib_cuts=None):
    """The stem cloud for a cut moved by `shift` m. Up: drop the band above the old cut. Down: add the band
    below the cut from the whole tree (`tree_xyz`); each band point goes to the nearest stem of the
    tree (`stems_xyz`, this stem is number `idx`) within `band_link` m of its points just above its own cut
    (`sib_cuts`, one per stem; default: all at z_cut), else it is root or base and left out."""
    if shift >= 0:
        return stem_xyz[stem_xyz[:, 2] >= z_cut + shift]
    from scipy.spatial import cKDTree
    band = tree_xyz[(tree_xyz[:, 2] >= z_cut + shift) & (tree_xyz[:, 2] < z_cut)]
    if not len(band):
        return stem_xyz
    cuts = [z_cut] * len(stems_xyz) if sib_cuts is None else sib_cuts
    dist = np.array([_dist_to(s[s[:, 2] < zc + 0.2], band) for s, zc in zip(stems_xyz, cuts)])
    mine = (dist.argmin(0) == idx) & (dist.min(0) <= band_link)
    return np.vstack([stem_xyz, band[mine]])


def _dist_to(pts, q):
    """Distance from each q to the nearest of pts (inf if pts is empty)."""
    from scipy.spatial import cKDTree
    return cKDTree(pts).query(q)[0] if len(pts) else np.full(len(q), np.inf)


def shift_cut_roots(roots_xyz, stems_xyz, z_cut, shift, band_link=0.3):
    """The roots (or base) cloud for a cut moved by `shift` m, the counterpart of shift_cut for the stems.
    `z_cut`: one cut, or one per stem (Avicennia stems start at their own fork).
    Up: the stems' points between the old and the new cut join the roots. Down: root points between the new
    and the old cut that lie within `band_link` m of a stem's points just above its cut go to that stem
    (the rest of the band stays root)."""
    cuts = np.broadcast_to(np.asarray(z_cut, float), (len(stems_xyz),))
    if shift > 0:
        band = [s[s[:, 2] < zc + shift] for s, zc in zip(stems_xyz, cuts)]
        return np.vstack([roots_xyz] + band)
    if shift == 0:
        return roots_xyz
    to_stem = np.zeros(len(roots_xyz), bool)
    for s, zc in zip(stems_xyz, cuts):
        in_band = (roots_xyz[:, 2] >= zc + shift) & (roots_xyz[:, 2] < zc) if np.ndim(z_cut) else roots_xyz[:, 2] >= zc + shift
        if in_band.any():
            to_stem[np.where(in_band)[0][_dist_to(s[s[:, 2] < zc + 0.2], roots_xyz[in_band]) <= band_link]] = True
    return roots_xyz[~to_stem]


# ---------------------------------------------------------------- measured radii (post-processing of an RCT tree file)
# RCT's skeleton (which segment goes where) is kept; each segment's radius is measured from the points around it.
# Only segments with RCT radius >= MEAS_MIN_R_RCT are measured. Fallback chain per segment: "fit" (robust circle perpendicular to the segment) -> "median" (median distance of
# its points to the segment axis; RCT centres its axes on the points, so this also works on one-sided scans)
# -> "rct" (RCT's radius). Leaves are not separated from wood: around twigs they inflate a point-based radius,
# so a measured radius above LEAF_CAP_FACTOR x RCT + LEAF_CAP_ABS counts as leaf contamination -> next level.
# The cap ties the measured radius to RCT's; it is only used where there are leaves (stems with crowns), not on
# the roots (on 199's roots it made the measured volume follow RCT's split-dependent radii: CV 16 % vs 20 %).
MEAS_MIN_R_RCT = 0.02      # m; segments with a thinner RCT radius keep it: there point noise and leaves (stems) or
                           # clutter (roots) dominate the points around the segment
MEAS_MIN_FIT = 10          # points needed for a circle fit
MEAS_MEDIAN_BIAS = 0.993   # median distance / clean-fit radius on segments with both (SU7_psp3_364_136, 598
                           # segments, MAD 1.7 %; analysis/median_fallback_check.py); the median level divides by it
MEAS_MIN_MEDIAN = 5        # points needed for the median distance
MEAS_MAX_REL_RMSE = 0.2    # circle fit rmse / radius above this: no clean section
MEAS_MIN_COVERAGE = 180    # degrees of the circle that must carry points
MEAS_MAX_DIST = 0.4        # m; points further from the segment axis are not used (no root or stem is that thick)
SHELL_FACTOR = 2.0         # a segment's own surface is the nearest shell of points around its axis: points beyond
                           # SHELL_FACTOR x its 20th-percentile distance belong to neighbouring roots, clutter or
                           # leaves (without this, 199's measured root area in the check slice was ~3x the points')
LEAF_CAP_FACTOR = 2.0
LEAF_CAP_ABS = 0.01        # m
SAMPLE_STEP = 0.02         # m between samples along the segments for the point-to-segment assignment


def segments(trees):
    """All segments (parent node -> node) of all trees: start, end, RCT radius (at the node), tree, node."""
    a, b, r, t_id, node = [], [], [], [], []
    for k, t in enumerate(trees):
        i = np.where(t["par"] >= 0)[0]
        a.append(t["xyz"][t["par"][i]]); b.append(t["xyz"][i]); r.append(t["r"][i])
        t_id.append(np.full(len(i), k)); node.append(i)
    return (np.vstack(a), np.vstack(b), np.concatenate(r), np.concatenate(t_id), np.concatenate(node))


def fit_circle_trimmed(p, iters=2):
    """Algebraic circle fit with iterative trimming of outliers (> 2.5 robust sd). Returns centre, r, rmse, mask."""
    keep = np.ones(len(p), bool)
    for it in range(iters + 1):
        q = p[keep]
        A = np.c_[2 * q, np.ones(len(q))]
        c = np.linalg.lstsq(A, (q ** 2).sum(1), rcond=None)[0]
        r = np.sqrt(max(c[2] + c[0] ** 2 + c[1] ** 2, 0))
        res = np.abs(np.hypot(*(p - c[:2]).T) - r)
        if it == iters:
            break
        mad = 1.4826 * np.median(res[keep])
        new = res <= max(2.5 * mad, 0.003)
        if new.sum() < 5 or (new == keep).all():
            break
        keep = new
    return c[:2], r, float(np.sqrt(np.mean(res[keep] ** 2))), keep


def measured_radii(trees, xyz, leaf_cap=True):
    """Measured radius per segment of all trees (see the fallback chain above); leaf_cap=False for leafless
    input (roots).
    Returns (per-segment DataFrame, list of per-tree radius arrays indexed like tree['r'])."""
    import pandas as pd
    A, B, r_rct, t_id, node = segments(trees)
    d = B - A
    L = np.linalg.norm(d, axis=1)
    axis = d / np.maximum(L, 1e-12)[:, None]
    # assign each point to the nearest segment (via samples along the segments)
    n_s = np.maximum(1, np.ceil(L / SAMPLE_STEP).astype(int))
    seg_of_sample = np.repeat(np.arange(len(A)), n_s)
    frac = (np.arange(n_s.sum()) - np.repeat(np.cumsum(n_s) - n_s, n_s) + 0.5) / np.repeat(n_s, n_s)
    samples = A[seg_of_sample] + frac[:, None] * d[seg_of_sample]
    seg = seg_of_sample[cKDTree(samples).query(xyz)[1]]
    rel = xyz - A[seg]
    along = (rel * axis[seg]).sum(1)
    radial = rel - along[:, None] * axis[seg]
    dist = np.linalg.norm(radial, axis=1)
    use = dist <= MEAS_MAX_DIST
    order = np.argsort(seg[use], kind="stable")
    seg_u, rad_u, dist_u = seg[use][order], radial[use][order], dist[use][order]
    starts = np.searchsorted(seg_u, np.arange(len(A)))
    ends = np.searchsorted(seg_u, np.arange(len(A)), side="right")

    r_meas = r_rct.copy()
    level = np.full(len(A), "rct", dtype=object)
    n_pts = ends - starts
    cap = LEAF_CAP_FACTOR * r_rct + LEAF_CAP_ABS if leaf_cap else np.full(len(A), MEAS_MAX_DIST)
    r_fit, r_q = np.full(len(A), np.nan), np.full((len(A), 3), np.nan)   # clean fit; shell distance q10, q25, median
    for s in np.where((n_pts >= MEAS_MIN_MEDIAN) & (r_rct >= MEAS_MIN_R_RCT))[0]:
        rs, ds = rad_u[starts[s]:ends[s]], dist_u[starts[s]:ends[s]]
        shell = ds <= SHELL_FACTOR * np.quantile(ds, 0.2) + 0.002
        rs, ds = rs[shell], ds[shell]
        if len(ds) < MEAS_MIN_MEDIAN:
            continue
        r_q[s] = np.quantile(ds, [0.1, 0.25, 0.5])
        if len(ds) >= MEAS_MIN_FIT:
            u = np.cross(axis[s], [1, 0, 0] if abs(axis[s][0]) < 0.9 else [0, 1, 0]); u /= np.linalg.norm(u)
            v = np.cross(axis[s], u)
            p2 = np.c_[rs @ u, rs @ v]
            c, r, rmse, keep = fit_circle_trimmed(p2)
            q = p2[keep] - c
            ang = np.sort(np.arctan2(q[:, 1], q[:, 0]))
            coverage = 360 - np.degrees(np.diff(np.r_[ang, ang[0] + 2 * np.pi]).max()) if len(ang) > 1 else 0
            if (r > 0 and rmse <= MEAS_MAX_REL_RMSE * r and coverage >= MEAS_MIN_COVERAGE
                    and np.linalg.norm(c) <= 0.5 * r + 0.01 and r <= cap[s]):
                r_fit[s] = r
                r_meas[s], level[s] = r, "fit"
                continue
        med = float(np.median(ds)) / MEAS_MEDIAN_BIAS
        if med <= cap[s]:
            r_meas[s], level[s] = med, "median"
    table = pd.DataFrame({"tree": t_id, "node": node, "length_m": L, "r_rct_m": r_rct, "r_meas_m": r_meas,
                          "level": level, "n_points": n_pts, "r_fit_m": r_fit, "r_q10_m": r_q[:, 0],
                          "r_q25_m": r_q[:, 1], "r_median_m": r_q[:, 2], "z_a": A[:, 2], "z_b": B[:, 2],
                          "ax": axis[:, 0], "ay": axis[:, 1], "az": axis[:, 2]})
    return table, radii_from_table(trees, table)


def radii_from_table(trees, table):
    """Per-tree radius arrays (indexed like tree['r']) with the measured radii of a segment table (measured_radii,
    or a saved <stem>_segments.csv of the same trees file)."""
    t_id, node, r_meas = table["tree"].to_numpy(), table["node"].to_numpy(), table["r_meas_m"].to_numpy()
    radii = []
    for k, t in enumerate(trees):
        rr = t["r"].copy()
        m = t_id == k
        rr[node[m]] = r_meas[m]
        base_kids = node[m][t["par"][node[m]] < 0]          # node 0 (the base) takes its thickest child's radius
        if len(base_kids):
            rr[t["par"] < 0] = rr[base_kids].max()
        radii.append(rr)
    return radii


def measured_volume(table):
    return float((np.pi * table["r_meas_m"] ** 2 * table["length_m"]).sum())


def level_shares(table):
    """Share of segments and of measured volume per fallback level."""
    vol = np.pi * table["r_meas_m"] ** 2 * table["length_m"]
    out = {}
    for lv in ("fit", "median", "rct"):
        m = table["level"] == lv
        out[f"seg_share_{lv}"] = round(float(m.mean()), 3)
        out[f"vol_share_{lv}"] = round(float(vol[m].sum() / vol.sum()), 3) if vol.sum() else np.nan
    return out


def diameter_with(tree, radii, z):
    """rct_diameter with other radii (e.g. measured) on the same skeleton."""
    return rct_diameter({"xyz": tree["xyz"], "r": radii, "par": tree["par"]}, z)


def area_profile(table, z_lo, z_hi, band=0.1, n_dir=36):
    """Per height band: total cross-section (pi r^2 of every segment crossing the band's mid-height, perpendicular
    to the segment) and frontal area (side-view silhouette of the segment parts inside the band, averaged over
    n_dir horizontal view directions), with measured and with RCT radii. Heights in the frame of the table."""
    import pandas as pd
    za, zb = table["z_a"].to_numpy(), table["z_b"].to_numpy()
    L = table["length_m"].to_numpy()
    ax = table[["ax", "ay", "az"]].to_numpy()
    phi = np.linspace(0, np.pi, n_dir, endpoint=False)
    nview = np.c_[np.cos(phi), np.sin(phi), np.zeros(n_dir)]
    side = np.sqrt(np.clip(1 - (ax @ nview.T) ** 2, 0, 1)).mean(1)   # |sin(axis, view)| averaged over views
    lo_z, hi_z = np.minimum(za, zb), np.maximum(za, zb)
    rows = []
    for z0 in np.arange(z_lo, z_hi, band):
        zm = z0 + band / 2
        cross = (lo_z <= zm) & (hi_z > zm)
        # length of each segment inside the band (segments are straight; share of their z extent)
        inside = np.clip(np.minimum(hi_z, z0 + band) - np.maximum(lo_z, z0), 0, None)
        dz = np.maximum(hi_z - lo_z, 1e-9)
        l_in = np.where(hi_z - lo_z > 1e-9, L * inside / dz, L * ((lo_z >= z0) & (lo_z < z0 + band)))
        row = {"z_lo": z0, "z_hi": z0 + band, "n_crossing": int(cross.sum())}
        for name, col in (("meas", "r_meas_m"), ("rct", "r_rct_m")):
            r = table[col].to_numpy()
            row[f"cross_section_{name}_m2"] = float((np.pi * r[cross] ** 2).sum())
            row[f"frontal_{name}_m2"] = float((2 * r * l_in * side).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def summarise(prefix, central, values, nd=4):
    """Central value, mean, sd, 2.5-97.5 % range and CV (%) of an ensemble, as log columns <prefix>_..."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    out = {f"{prefix}_central": round(float(central), nd)}
    if len(v):
        mean, sd = float(v.mean()), float(v.std(ddof=1)) if len(v) > 1 else np.nan
        out |= {f"{prefix}_mean": round(mean, nd), f"{prefix}_sd": round(sd, nd + 1),
                f"{prefix}_p2.5": round(float(np.quantile(v, 0.025)), nd),
                f"{prefix}_p97.5": round(float(np.quantile(v, 0.975)), nd),
                f"{prefix}_cv_pct": round(100 * sd / mean, 1) if mean else np.nan}
    return out


RADIUS_CLASSES = {"thick": (0.05, np.inf), "mid": (0.01, 0.05), "twig": (0.0, 0.01)}   # by RCT radius, m


def volume_by_class(table):
    """RCT and measured volume per RCT radius class (where a measured volume differs: trunk, branches, twigs)."""
    out = {}
    for name, (lo, hi) in RADIUS_CLASSES.items():
        m = (table["r_rct_m"] >= lo) & (table["r_rct_m"] < hi)
        L = table.loc[m, "length_m"]
        out[f"vol_rct_{name}_m3"] = round(float((np.pi * table.loc[m, "r_rct_m"] ** 2 * L).sum()), 4)
        out[f"vol_meas_{name}_m3"] = round(float((np.pi * table.loc[m, "r_meas_m"] ** 2 * L).sum()), 4)
    return out


def cylinder_polygons(tree, axis, radii=None):
    """Side outline of every segment (parent -> node) projected on the plane (axis, z): a quadrilateral whose
    width is the segment diameter, perpendicular to the segment in that plane. For figures (03b, 04b)."""
    r = tree["r"] if radii is None else radii
    polys = []
    for i, p in enumerate(tree["par"]):
        if p < 0:
            continue
        a = np.array([tree["xyz"][p, axis], tree["xyz"][p, 2]])
        b = np.array([tree["xyz"][i, axis], tree["xyz"][i, 2]])
        d = b - a
        n = np.array([-d[1], d[0]]) / (np.linalg.norm(d) or 1)
        polys.append([a + n * r[p], b + n * r[i], b - n * r[i], a - n * r[p]])
    return polys
