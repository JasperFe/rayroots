"""Checks of vol_area (rayroots/single/area.py): synthetic rings and an oval root, the SU6_psp4_17 arch sections,
the area profile of the points against the cylinders, root totals, and a QC figure of the segments that change most.
Usage: area_check.py <rayroots-single --out folder> <out dir> [--arch-tree SU6_psp4_17]
Run in the conda env rayroots, after rayroots-single model."""
# %% Imports
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Circle, Ellipse  # noqa: E402
from scipy.spatial import ConvexHull, cKDTree  # noqa: E402

from rayroots.single import area  # noqa: E402
from rayroots.single import rct_tools as rt  # noqa: E402

C_RCT, C_MEAS, C_AREA, GREY = "#c8642a", "#2a6db0", "#1a9850", "#8c8c8c"
ARCH = (1, [34, 48, 65])            # tree index and nodes of the arch sections a-c (ref_scripts/arch.py)
BAND, NDIR, RASTER = 0.10, 36, 0.01
rng = np.random.default_rng(0)


def table(df):
    """Markdown table (no tabulate in the env)."""
    rows = [list(map(str, df.columns)), ["---"] * len(df.columns)] + [[str(v) for v in r] for r in df.itertuples(index=False)]
    return "\n".join("| " + " | ".join(r) + " |" for r in rows)


# %% Synthetic
def ring(a, b, cov, n=400, noise=0.005):
    th = rng.uniform(0, np.radians(cov), n) + rng.uniform(0, 2 * np.pi)
    p = np.c_[a * np.cos(th), b * np.sin(th)]
    ang = rng.uniform(0, np.pi)
    R = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    return p @ R.T + rng.normal(0, noise, p.shape)


def synthetic_rings(reps=20):
    rows = []
    for r_eq in (0.10, 0.05, 0.04):
        for ratio in (1, 1.5, 2):
            for cov in (360, 180):
                b = r_eq / np.sqrt(ratio)
                a = b * ratio
                res = [area.slab_area(ring(a, b, cov), b)[:2] for _ in range(reps)]
                A = np.array([x[0] for x in res])
                rows.append({"r_eq_cm": r_eq * 100, "ratio": ratio, "coverage_deg": cov, "true_cm2": round(np.pi * a * b * 1e4, 1),
                             "median_cm2": round(np.median(A) * 1e4, 1), "err_pct": round(100 * (np.median(A) / (np.pi * a * b) - 1), 1),
                             "p10_pct": round(100 * (np.quantile(A, 0.1) / (np.pi * a * b) - 1), 1),
                             "p90_pct": round(100 * (np.quantile(A, 0.9) / (np.pi * a * b) - 1), 1),
                             "n_ellipse": sum(x[1] == "ellipse" for x in res), "n": reps})
    return pd.DataFrame(rows)


def synthetic_root(ratio, cov, a_eq=0.06, length=1.0, seg_len=0.05, density=40000):
    """A straight oval root along x, scanned on `cov` degrees, skeleton as RCT would give it (circle = minor axis).
    Returns true volume, vol_meas-like circle volume, vol_area."""
    b = a_eq / np.sqrt(ratio)
    a = b * ratio
    n = int(density * length * np.pi * (a + b))
    x = rng.uniform(0, length, n)
    th = rng.uniform(0, np.radians(cov), n)
    p = np.c_[x, a * np.cos(th), b * np.sin(th)] + rng.normal(0, 0.005, (n, 3))
    nodes = np.arange(0, length + 1e-9, seg_len)
    xyz = np.c_[nodes, np.zeros((len(nodes), 2))]
    tree = {"xyz": xyz, "r": np.full(len(nodes), b), "par": np.r_[-1, np.arange(len(nodes) - 1)]}
    segs = rt.segments([tree])
    tab = pd.DataFrame({"length_m": np.linalg.norm(segs[1] - segs[0], axis=1), "r_meas_m": np.full(len(segs[0]), b),
                        "level": "fit"})
    out, _ = area.root_areas([tree], p, tab, segs)
    return np.pi * a * b * length, float((np.pi * b ** 2 * tab.length_m).sum()), float((out.area_m2 * out.length_m).sum())


# %% A tree's roots from the rayroots-single output
def load(single, tree):
    w = single / "work"
    log = pd.read_csv(w / "split_log.csv").set_index("tree").loc[tree]
    z_cut, z_ref = float(log["z_cut"]), float(log["z_ref"])
    pts = rt.read_xyz(w / "roots" / f"{tree}_roots.ply")
    pts = np.c_[pts[:, :2], 2 * z_cut - pts[:, 2]]                    # the flipped frame RCT worked in
    d = w / "rct" / tree
    trees = rt.read_trees(d / f"{tree}_roots_trees.txt")
    tab = pd.read_csv(d / f"{tree}_roots_segments.csv")
    segs = rt.segments(trees)
    assert len(tab) == len(segs[0]) and (tab.node.to_numpy() == segs[4]).all()
    return {"pts": pts, "trees": trees, "tab": tab, "segs": segs, "z_cut": z_cut, "z_ref": z_ref}


def crossing(segs, M, d, radius):
    A, B = segs[0], segs[1]
    ta, tb = (A - M) @ d, (B - M) @ d
    j = np.flatnonzero((ta * tb <= 0) & (ta != tb))
    X = A[j] + (B[j] - A[j]) * (ta[j] / (ta[j] - tb[j]))[:, None] - M
    k = np.linalg.norm(X, axis=1) <= radius
    return j[k], X[k]


def arch_sections(R):
    """As ref_scripts/arch.py: slab |t| < 4 cm, clumps removed, convex hull; cylinders crossing the plane within 0.2 m."""
    t_id, node = R["segs"][3], R["segs"][4]
    kd = cKDTree(R["pts"])
    rows, shown = [], []
    for i, n in enumerate(ARCH[1]):
        s = int(np.flatnonzero((t_id == ARCH[0]) & (node == n))[0])
        A, B = R["segs"][0][s], R["segs"][1][s]
        d = (B - A) / np.linalg.norm(B - A)
        M = (A + B) / 2
        u, w = area.frame(d)
        q = R["pts"][kd.query_ball_point(M, 0.3)] - M
        t = q @ d
        rad = q - np.outer(t, d)
        m = (np.abs(t) < 0.04) & (np.linalg.norm(rad, axis=1) < 0.2)
        xy = np.c_[rad[m] @ u, rad[m] @ w]
        kk = cKDTree(xy)
        xy = xy[np.array([len(kk.query_ball_point(p, 0.02)) >= 5 for p in xy])]
        hull = ConvexHull(xy)
        j, X = crossing(R["segs"], M, d, 0.2)
        tb = R["tab"].iloc[j]
        rows.append({"section": "abc"[i], "node": n, "outline_cm2": round(hull.volume * 1e4),
                     "rct_cm2": round((np.pi * tb.r_rct_m ** 2).sum() * 1e4), "meas_cm2": round((np.pi * tb.r_meas_m ** 2).sum() * 1e4),
                     "area_cm2": round(tb.area_m2.sum() * 1e4), "n_cylinders": len(j),
                     "level": R["tab"].area_level.iat[s], "axis_ratio": R["tab"].axis_ratio.iat[s]})
        shown.append((s, xy, hull, np.c_[X @ u, X @ w], tb))
    return pd.DataFrame(rows), shown


# %% Area profile: side-view outline per 10 cm band, points against cylinders
def surface_samples(A, B, r, step=0.005):
    """Points on the mantle of each cylinder, every `step` along and around."""
    out = []
    for a, b, rr in zip(A, B, r):
        d = b - a
        L = np.linalg.norm(d)
        if L == 0 or rr <= 0:
            continue
        ax = d / L
        u, v = area.frame(ax)
        na, nr = max(2, int(np.ceil(L / step)) + 1), max(8, int(np.ceil(2 * np.pi * rr / step)))
        s = np.linspace(0, L, na)[:, None, None]
        th = np.linspace(0, 2 * np.pi, nr, endpoint=False)[None, :, None]
        out.append((a + s * ax + rr * (np.cos(th) * u + np.sin(th) * v)).reshape(-1, 3))
    return np.vstack(out)


def outline(p, h, bands):
    """Side-view silhouette area per band (occupied RASTER cells), mean over NDIR horizontal view directions."""
    b = np.floor(h / BAND).astype(int)
    keep = (b >= 0) & (b < len(bands))
    p, h, b = p[keep], h[keep], b[keep]
    iz = np.floor(h / RASTER).astype(np.int64)
    out = np.zeros(len(bands))
    for phi in np.linspace(0, np.pi, NDIR, endpoint=False):
        s = np.floor((p[:, 0] * np.cos(phi) + p[:, 1] * np.sin(phi)) / RASTER).astype(np.int64)
        key = np.unique((iz - iz.min()) * (1 << 32) + (s - s.min()))
        zb = ((key >> 32) + iz.min()) // round(BAND / RASTER)
        out += np.bincount(np.clip(zb, 0, len(bands) - 1), minlength=len(bands)) * RASTER ** 2
    return out / NDIR


def area_profile(R):
    A, B = R["segs"][0].copy(), R["segs"][1].copy()
    for X in (A, B):                                              # back to the real frame
        X[:, 2] = 2 * R["z_cut"] - X[:, 2]
    pts = R["pts"].copy()
    pts[:, 2] = 2 * R["z_cut"] - pts[:, 2]
    top = R["z_cut"] - R["z_ref"]
    bands = np.arange(0, top, BAND)
    prof = pd.DataFrame({"h_lo_m": bands.round(2), "points_m2": outline(pts, pts[:, 2] - R["z_ref"], bands)})
    tab = R["tab"]
    for name, r in (("rct", tab.r_rct_m), ("meas", tab.r_meas_m), ("area", np.sqrt(tab.area_m2 / np.pi))):
        q = surface_samples(A, B, r.to_numpy())
        prof[f"{name}_m2"] = outline(q, q[:, 2] - R["z_ref"], bands)
    return prof


# %% QC figure: the 3 segments whose area changes most
def qc_figure(R, tree, out):
    tab = R["tab"]
    gain = (tab.area_m2 - np.pi * tab.r_meas_m ** 2) * tab.length_m
    top = np.argsort(-np.abs(gain.to_numpy()))[:3]
    _, slabs = area.root_areas(R["trees"], R["pts"], tab[["r_meas_m", "level", "length_m"]], R["segs"])
    fig, ax = plt.subplots(1, 3, figsize=(12, 4.4))
    for a, s in zip(ax, top):
        sl = [x for x in slabs if x["seg"] == s][1]                   # the middle slab
        p2 = sl["points"]
        a.scatter(p2[:, 0], p2[:, 1], s=2, c=GREY, lw=0)
        r = tab.iloc[s]
        a.add_patch(Circle((0, 0), r.r_rct_m, fc=C_RCT, ec="none", alpha=0.3))
        a.add_patch(Circle((0, 0), r.r_meas_m, fill=False, ec=C_MEAS, lw=1.8))
        e = sl["ellipse"]
        if e is not None:
            a.add_patch(Ellipse(e["centre"], 2 * e["a"], 2 * e["b"], angle=np.degrees(e["angle"]), fill=False,
                                ec=C_AREA, lw=1.8, ls="-" if sl["level"] == "ellipse" else ":"))
        if sl["filled"] is not None:
            f, o = sl["filled"], sl["origin"]
            xs = (np.arange(f.shape[0]) + o[0] + 0.5) * area.CELL
            ys = (np.arange(f.shape[1]) + o[1] + 0.5) * area.CELL
            a.contour(xs, ys, f.T.astype(float), levels=[0.5], colors="k", linewidths=0.6, linestyles="--")
        lim = max(0.12, 1.3 * (e["a"] if e is not None else r.r_meas_m))
        a.set_xlim(-lim, lim)
        a.set_ylim(-lim, lim)
        a.set_aspect("equal")
        a.set_xticks([])
        a.set_yticks([])
        a.set_title(f"tree {r.tree} node {r.node}: {r.area_level}, ratio {r.axis_ratio:.2f}", fontsize=9)
        a.set_xlabel(f"RCT {np.pi * r.r_rct_m ** 2 * 1e4:.0f}  meas {np.pi * r.r_meas_m ** 2 * 1e4:.0f}  "
                     f"area {r.area_m2 * 1e4:.0f} cm²  (own {r.area_own_m2 * 1e4:.0f})", fontsize=9)
    fig.suptitle(f"{tree}: segments whose root volume changes most. Orange RCT, blue measured circle, green ellipse "
                 f"(dotted: rejected), dashed blob", fontsize=9)
    fig.savefig(out, dpi=120, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def profile_figure(prof, tree, out):
    fig, ax = plt.subplots(figsize=(4.5, 5))
    h = prof.h_lo_m + BAND / 2
    ax.plot(prof.points_m2, h, "k", lw=1.6, label="points")
    for name, c in (("rct", C_RCT), ("meas", C_MEAS), ("area", C_AREA)):
        ax.plot(prof[f"{name}_m2"], h, c=c, lw=1.2, label=f"vol_{name} cylinders")
    ax.set_xlabel("side-view outline per 10 cm band (m²)")
    ax.set_ylabel("height above ground (m)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title(tree, fontsize=10, loc="left")
    fig.savefig(out, dpi=120, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# %% Run
if __name__ == "__main__":
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("single", type=Path)
    a.add_argument("out", type=Path)
    a.add_argument("--arch-tree", default="SU6_psp4_17")
    x = a.parse_args()
    x.out.mkdir(parents=True, exist_ok=True)
    md = ["# vol_area checks", ""]
    syn = synthetic_rings()
    md += ["## Synthetic rings (5 mm noise, 20 rings each; ellipse area of one slab)", "", table(syn), ""]
    rows = []
    for ratio in (1, 1.5, 2):
        for cov in (360, 180):
            v_true, v_circ, v_area = synthetic_root(ratio, cov)
            rows.append({"ratio": ratio, "coverage_deg": cov, "true_m3": round(v_true, 5), "circle_minor_m3": round(v_circ, 5),
                         "vol_area_m3": round(v_area, 5), "err_pct": round(100 * (v_area / v_true - 1), 1)})
    md += ["## Synthetic oval root, 1 m, whole method (skeleton circle = minor axis)", "",
           table(pd.DataFrame(rows)), ""]
    vol = pd.read_csv(x.single / "tree_volumes.csv")
    vol = vol[vol["row"] == "tree"]
    tot, prof_dev = [], []
    for tree in vol.tree:
        R = load(x.single, tree)
        prof = area_profile(R)
        prof.to_csv(x.out / f"{tree}_area_profile.csv", index=False)
        profile_figure(prof, tree, x.out / f"{tree}_area_profile.png")
        qc_figure(R, tree, x.out / f"{tree}_area_qc.png")
        p = prof.points_m2.sum()
        prof_dev.append({"tree": tree, **{f"{m}_vs_points_pct": round(100 * (prof[f"{m}_m2"].sum() / p - 1), 1)
                                          for m in ("rct", "meas", "area")},
                         "bands_area_closer_than_meas": f"{int(((prof.area_m2 - prof.points_m2).abs() < (prof.meas_m2 - prof.points_m2).abs()).sum())}/{len(prof)}",
                         "bands_area_above_points": int((prof.area_m2 > prof.points_m2).sum())})
        lv = R["tab"][R["tab"].r_meas_m >= area.MIN_R].area_level.value_counts().to_dict()
        r = vol[vol.tree == tree].iloc[0]
        tot.append({"tree": tree, "vol_roots_rct_m3": r.vol_roots_rct_m3, "vol_roots_meas_m3": r.vol_roots_meas_m3,
                    "vol_roots_area_m3": r.vol_roots_area_m3, "area/meas": round(r.vol_roots_area_m3 / r.vol_roots_meas_m3, 3),
                    "segments >= 3 cm": int((R["tab"].r_meas_m >= area.MIN_R).sum()), "levels": lv,
                    "vol_rct_m3": r.vol_rct_m3, "vol_meas_m3": r.vol_meas_m3, "vol_area_m3": r.vol_area_m3})
        if tree == x.arch_tree:
            arch, _ = arch_sections(R)
            md += [f"## Arch sections a-c, {tree} (as ref_scripts/arch.py)", "", table(arch), ""]
    md += ["## Root totals", "", table(pd.DataFrame(tot)), "",
           "## Area profile: side-view outline summed over the 10 cm bands, against the points", "",
           table(pd.DataFrame(prof_dev)), ""]
    (x.out / "area_check.md").write_text("\n".join(md))
    print("\n".join(md))
