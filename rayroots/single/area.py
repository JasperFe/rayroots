"""vol_area (area v1): root volume from the cross-section area of the points, for oval and fused roots.
Roots only, on the central run. vol_rct and vol_meas are not touched.

Per root segment with r_meas >= MIN_R: three thin slabs perpendicular to the RCT axis (25 / 50 / 75 % of its length).
In each slab all root points (not only the segment's own), rastered in the slab plane; the blob that contains the
axis (or the nearest one) is kept and an ellipse is fitted to its points (fits the middle of the noisy ring, works on
one-sided arcs). Slab area = the ellipse when it covers >= MIN_COVERAGE degrees and is plausible, else the measured
circle. Segment area = median of the slabs (>= MIN_SLABS ellipses, else the circle). Segments that cross the same
blob (stacked roots, twigs inside an oval) share its area: larger circle first, no more than the blob in total.
Thinner segments keep pi r_meas^2. vol_area = sum of area x length."""
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

VERSION = "area_v1"
MIN_R = 0.03            # m, r_meas from which the slab areas are used
SLAB = 0.02             # m, slab thickness
SLAB_AT = (0.25, 0.5, 0.75)
CELL = 0.01             # m, raster in the slab plane
MAX_DIST = 0.4          # m, points further from the axis are not used (as MEAS_MAX_DIST)
NEAR = 0.02             # m, a blob that does not contain the axis is used when this close to the circle r_meas
MIN_POINTS = 12
MIN_COVERAGE = 180      # degrees of the ellipse that carry points
MAX_RATIO = 4.0         # axis ratio
MAX_GROWTH = 4.0        # ellipse area / pi r_meas^2
MIN_SLABS = 2
TRIM = 2.5              # x rmse, points dropped before the refit


# %% Ellipse fit
def fit_conic(p):
    """Direct least-squares ellipse (Fitzgibbon; numerically stable form of Halir & Flusser). Returns the conic
    coefficients (a, b, c, d, e, f) of a x^2 + b xy + c y^2 + d x + e y + f = 0, or None."""
    x, y = p[:, 0], p[:, 1]
    D1 = np.c_[x * x, x * y, y * y]
    D2 = np.c_[x, y, np.ones_like(x)]
    S1, S2, S3 = D1.T @ D1, D1.T @ D2, D2.T @ D2
    try:
        T = -np.linalg.solve(S3, S2.T)
    except np.linalg.LinAlgError:
        return None
    M = S1 + S2 @ T
    M = np.array([M[2] / 2, -M[1], M[0] / 2])
    w, v = np.linalg.eig(M)
    v = np.real(v)
    cond = 4 * v[0] * v[2] - v[1] ** 2
    k = np.flatnonzero(cond > 0)
    if not len(k):
        return None
    a1 = v[:, k[0]]
    return np.r_[a1, T @ a1]


def conic_geometry(q):
    """Centre, semi-axes (major, minor) and angle of the major axis of an ellipse conic; None if not an ellipse."""
    a, b, c, d, e, f = q
    den = b * b - 4 * a * c
    if den >= 0:
        return None
    x0 = (2 * c * d - b * e) / den
    y0 = (2 * a * e - b * d) / den
    M = np.array([[a, b / 2], [b / 2, c]])
    f0 = a * x0 ** 2 + b * x0 * y0 + c * y0 ** 2 + d * x0 + e * y0 + f
    w, v = np.linalg.eigh(M)
    if f0 == 0 or np.any(-f0 / w <= 0):
        return None
    ax = np.sqrt(-f0 / w)                       # semi-axis along each eigenvector
    i = int(np.argmax(ax))
    return np.array([x0, y0]), float(ax[i]), float(ax[1 - i]), float(np.arctan2(v[1, i], v[0, i]))


def fit_ellipse(p, iters=2):
    """Ellipse fit with trimming of far points. Returns dict(centre, a, b, angle, rmse, keep, coverage_deg) or None."""
    keep = np.ones(len(p), bool)
    g = None
    for _ in range(iters + 1):
        if keep.sum() < MIN_POINTS:
            return None
        q = fit_conic(p[keep])
        g = conic_geometry(q) if q is not None else None
        if g is None:
            return None
        c, A, B, ang = g
        u = np.c_[np.cos(ang), np.sin(ang)][0]
        r = (p - c) @ np.c_[u, [-u[1], u[0]]]      # in the ellipse frame
        rho = np.hypot(r[:, 0] / A, r[:, 1] / B)   # 1 on the ellipse
        dist = np.abs(rho - 1) * np.sqrt(A * B)    # approximate distance to the ellipse
        rmse = float(np.sqrt(np.mean(dist[keep] ** 2)))
        new = dist <= max(TRIM * rmse, 0.003)
        if np.array_equal(new, keep):
            break
        keep = new
    c, A, B, ang = refine(p[keep], c, A, B, ang)
    u = np.c_[np.cos(ang), np.sin(ang)][0]
    r = (p - c) @ np.c_[u, [-u[1], u[0]]]
    rmse = float(np.sqrt(np.mean(sampson(r[keep], A, B) ** 2)))
    t = np.sort(np.arctan2(r[keep, 1] / B, r[keep, 0] / A))
    cov = 360 - np.degrees(np.diff(np.r_[t, t[0] + 2 * np.pi]).max()) if len(t) > 1 else 0.0
    return {"centre": c, "a": A, "b": B, "angle": ang, "rmse": rmse, "keep": keep, "coverage_deg": float(cov)}


def sampson(r, A, B):
    """Distance of points (in the ellipse frame) to the ellipse, first order (Sampson)."""
    F = (r[:, 0] / A) ** 2 + (r[:, 1] / B) ** 2 - 1
    g = 2 * np.hypot(r[:, 0] / A ** 2, r[:, 1] / B ** 2)
    return F / np.maximum(g, 1e-12)


def refine(p, c, A, B, ang):
    """Geometric fit from the algebraic one: the algebraic fit shrinks the ellipse on partial arcs with noise."""
    from scipy.optimize import least_squares

    def res(x):
        u = np.r_[np.cos(x[4]), np.sin(x[4])]
        r = (p - x[:2]) @ np.c_[u, [-u[1], u[0]]]
        return sampson(r, np.exp(x[2]), np.exp(x[3]))

    x = least_squares(res, np.r_[c, np.log(A), np.log(B), ang], method="lm").x
    A, B, ang = np.exp(x[2]), np.exp(x[3]), x[4]
    if B > A:
        A, B, ang = B, A, ang + np.pi / 2
    return x[:2], float(A), float(B), float(ang)


# %% One slab
def blob_of_axis(p2, r_circle):
    """Raster the slab points (2D, axis at 0,0); the blob containing the axis, or the nearest within r_circle + NEAR.
    Returns (bool mask of the points in the blob, filled blob raster, raster origin) or None."""
    ij = np.floor(p2 / CELL).astype(int)
    o = ij.min(0) - 2
    ij -= o
    img = np.zeros(ij.max(0) + 3, bool)
    img[ij[:, 0], ij[:, 1]] = True
    img = ndimage.binary_closing(img, np.ones((3, 3)), border_value=0) | img      # bridge 1 cm gaps
    lab, n = ndimage.label(img, np.ones((3, 3)))
    if n == 0:
        return None
    a = np.floor(-o.astype(float)).astype(int)
    pick = 0
    for k in range(1, n + 1):
        filled = ndimage.binary_fill_holes(lab == k)
        if 0 <= a[0] < img.shape[0] and 0 <= a[1] < img.shape[1] and filled[a[0], a[1]]:
            pick = k
            break
    if not pick:
        cells = np.argwhere(lab > 0)
        dist = np.hypot(*((cells + 0.5) * CELL + o * CELL).T)
        j = int(np.argmin(dist))
        if dist[j] > r_circle + NEAR:
            return None
        pick = lab[tuple(cells[j])]
    filled = ndimage.binary_fill_holes(lab == pick)
    return lab[ij[:, 0], ij[:, 1]] == pick, filled, o


def slab_area(p2, r_circle):
    """One slab: (area, level, axis ratio, ellipse dict or None, filled raster, origin)."""
    circle = np.pi * r_circle ** 2
    if len(p2) < MIN_POINTS:
        return circle, "circle", np.nan, None, None, None
    b = blob_of_axis(p2, r_circle)
    if b is None:
        return circle, "circle", np.nan, None, None, None
    m, filled, o = b
    e = fit_ellipse(p2[m]) if m.sum() >= MIN_POINTS else None
    if e is None:
        return circle, "circle", np.nan, None, filled, o
    area = np.pi * e["a"] * e["b"]
    ci = np.floor(e["centre"] / CELL).astype(int) - o
    inside = 0 <= ci[0] < filled.shape[0] and 0 <= ci[1] < filled.shape[1] and filled[ci[0], ci[1]]
    ok = (e["coverage_deg"] >= MIN_COVERAGE and e["a"] / e["b"] <= MAX_RATIO and area <= MAX_GROWTH * circle
          and e["a"] <= MAX_DIST and (inside or np.linalg.norm(e["centre"]) <= r_circle))
    if not ok:
        return circle, "circle", e["a"] / e["b"], e, filled, o
    return area, "ellipse", e["a"] / e["b"], e, filled, o


def frame(axis):
    u = np.cross(axis, [1, 0, 0] if abs(axis[0]) < 0.9 else [0, 1, 0])
    u /= np.linalg.norm(u)
    return u, np.cross(axis, u)


# %% All root segments
def root_areas(trees, xyz, table, segs):
    """Adds area_own_m2, area_m2, area_level, axis_ratio, n_slabs_ok to the segments table (same row order as
    rct_tools.segments). segs = rct_tools.segments(trees). Returns (table, slabs) with slabs: one dict per slab
    (for the shared areas and the QC figure)."""
    A, B, _, t_id, node = segs
    d = B - A
    L = np.linalg.norm(d, axis=1)
    axis = d / np.maximum(L, 1e-12)[:, None]
    r_meas = table["r_meas_m"].to_numpy()
    base_level = np.where(table["level"].to_numpy() == "rct", "rct", "circle")
    own = np.pi * r_meas ** 2
    level = base_level.astype(object)
    ratio = np.full(len(A), np.nan)
    n_ok = np.zeros(len(A), int)
    kd = cKDTree(xyz)
    slabs = []
    for s in np.flatnonzero(r_meas >= MIN_R):
        u, v = frame(axis[s])
        areas, ratios = [], []
        for f in SLAB_AT:
            M = A[s] + f * d[s]
            reach = min(MAX_DIST, MAX_GROWTH * r_meas[s] + NEAR)
            q = xyz[kd.query_ball_point(M, np.hypot(reach, SLAB / 2))] - M
            t = q @ axis[s]
            q = q[np.abs(t) < SLAB / 2]
            p2 = np.c_[q @ u, q @ v]
            p2 = p2[np.hypot(*p2.T) <= reach]
            area, lv, rt_, e, filled, o = slab_area(p2, r_meas[s])
            slabs.append({"seg": s, "M": M, "axis": axis[s], "u": u, "v": v, "area": area, "level": lv, "ellipse": e,
                          "filled": filled, "origin": o, "points": p2})
            if lv == "ellipse":
                areas.append(area)
                ratios.append(rt_)
        n_ok[s] = len(areas)
        if len(areas) >= MIN_SLABS:
            own[s] = float(np.median(areas))
            level[s] = "ellipse"
            ratio[s] = float(np.median(ratios))
    final = share(own, np.pi * r_meas ** 2, slabs, A, B, t_id, node, trees)
    out = table.copy()
    out["area_own_m2"], out["area_m2"], out["area_level"] = own, final, level
    out["axis_ratio"], out["n_slabs_ok"] = np.round(ratio, 3), n_ok
    return out, slabs


def share(own, circle, slabs, A, B, t_id, node, trees):
    """Segments crossing the same ellipse blob share its area, larger circle first. A segment's neighbours on the
    same RCT branch (parent, children) are not counted as crossing its slabs."""
    final = own.copy()
    key = {(int(t), int(n)): i for i, (t, n) in enumerate(zip(t_id, node))}
    par = np.array([key.get((int(t), int(trees[t]["par"][n])), -1) for t, n in zip(t_id, node)])
    for sl in slabs:
        if sl["level"] != "ellipse" or own[sl["seg"]] <= 0:
            continue
        s, M, ax = sl["seg"], sl["M"], sl["axis"]
        ta, tb = (A - M) @ ax, (B - M) @ ax
        j = np.flatnonzero((ta * tb <= 0) & (ta != tb))
        j = j[(j != s) & (par[j] != s) & (j != par[s])]
        if not len(j):
            continue
        X = A[j] + (B[j] - A[j]) * (ta[j] / (ta[j] - tb[j]))[:, None] - M
        ci = np.floor(np.c_[X @ sl["u"], X @ sl["v"]] / CELL).astype(int) - sl["origin"]
        f = sl["filled"]
        ins = (ci[:, 0] >= 0) & (ci[:, 1] >= 0) & (ci[:, 0] < f.shape[0]) & (ci[:, 1] < f.shape[1])
        ins[ins] = f[ci[ins, 0], ci[ins, 1]]
        members = np.r_[s, j[ins]]
        if len(members) == 1:
            continue
        left = sl["area"]
        for m in members[np.argsort(-circle[members], kind="stable")]:
            give = min(own[m], left)
            final[m] = min(final[m], give)
            left -= give
    return final
