"""One QC figure per tree (B_trees/qc/<tree>.png): (a) side view with slices, last root and cut; (b) blobs per height;
(c) top views of a slice below the last root and at the cut; (d) stem and flipped roots with RCT cylinders, once modelled."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

from rayroots.single import rct_tools as rt  # noqa: E402

CUT, ROOT = "#c0143c", "#2a78d6"
LIGHT, MID, WOOD = "#c9c9c9", "#8c8c8c", "#5b3a1e"
N_SHOW = 150_000
SLICE = 0.05
TOP_H = 3.0                 # m shown above the cut in (a) and (d)
NV = 12                     # slices drawn in (a)
PANEL_H = 3.6               # inches
TOP_W = 2.0                 # inches, the two slices
MAX_W = 4.0                 # inches, widest side-view panel


def sub(a, n=N_SHOW):
    return a[np.random.default_rng(0).permutation(len(a))[:n]]


def cylinders(ax, trees, shift_x=0.0, fz=lambda z: z, top=np.inf, color=WOOD):
    q = []
    for t in trees:
        for i in np.flatnonzero(t["par"] >= 0):
            a, b = t["xyz"][t["par"][i]], t["xyz"][i]
            pa, pb = np.array([a[0] + shift_x, fz(a[2])]), np.array([b[0] + shift_x, fz(b[2])])
            n = np.hypot(*(pb - pa))
            if n == 0 or min(pa[1], pb[1]) > top:
                continue
            w = np.array([-(pb - pa)[1], (pb - pa)[0]]) / n * t["r"][i]
            q.append([pa - w, pa + w, pb + w, pb - w])
    ax.add_collection(PolyCollection(q, facecolors=color, edgecolors="none"))


def top_view(ax, xyz, z, title, c, half):
    """One slice; c, half: centre and half width of the window, the same for both slices."""
    k = np.abs(xyz[:, 2] - z) < SLICE / 2
    ax.scatter(xyz[k, 0], xyz[k, 1], s=0.6, c=MID, lw=0, rasterized=True)
    ax.set_xlim(c[0] - half, c[0] + half)
    ax.set_ylim(c[1] - half, c[1] + half)
    ax.set_aspect("equal", anchor="N")
    ax.axis("off")
    ax.set_title(title, fontsize=9)


def qc(s, r):
    """r: the tree's split_log row. Draws (d) when the model exists."""
    tree, z_ref, z_cut = r["tree"], float(r["z_ref"]), float(r["z_cut"])
    xyz = rt.read_xyz(s / "standard" / f"{tree}.ply")
    prof = pd.read_csv(s / "profile" / f"{tree}_profile.csv")
    rooted = r["split"] == "roots"
    h_c = z_cut - z_ref
    if rooted:
        h0 = float(r["z0"]) - z_ref
        h, n = prof.h_above_z0 + h0, prof.n_blobs
        h_r = float(r["h_last_root"]) + h0
    else:
        h, n = prof.h_above_zref, prof.n_tracks
        h_r = 1.3
    x, hag = xyz[:, 0], xyz[:, 2] - z_ref
    top = h_c + TOP_H
    low = hag < h_c
    xa = (float(x[hag < top].min()), float(x[hag < top].max()))
    y_top = max(h_c, TOP_H) if rooted else TOP_H             # (d): height above the cut
    stem = (hag >= h_c) & (hag < h_c + y_top)
    gap_d = 0.6
    xs0, ws = float(x[stem].min()), float(np.ptp(x[stem]))
    xr0, wr = (float(x[low].min()), float(np.ptp(x[low]))) if rooted and low.any() else (0.0, 0.0)
    xd = (0.0, ws + (gap_d + wr if wr else 0.0))
    yd = (-0.15, y_top)
    sa = min(PANEL_H / top, MAX_W / (xa[1] - xa[0]))         # inches per metre, one scale per panel
    sd = min(PANEL_H / (yd[1] - yd[0]), MAX_W / (xd[1] - xd[0]))
    wa, ha, wd, hd = sa * (xa[1] - xa[0]), sa * top, sd * (xd[1] - xd[0]), sd * (yd[1] - yd[0])
    w = [wa, 1.0, TOP_W, TOP_W, wd]
    gap = [0.65, 0.2, 0.2, 0.3]
    W, H = sum(w) + sum(gap) + 0.2, PANEL_H + 0.75
    left = 0.1 + np.r_[0, np.cumsum(np.array(w[:-1]) + gap)]
    fig = plt.figure(figsize=(W, H))

    def panel(i, height):
        return fig.add_axes([left[i] / W, (H - 0.3 - height) / H, w[i] / W, height / H])

    a = panel(0, ha)
    k = hag < top
    q = sub(np.c_[x[k], hag[k]])
    a.scatter(q[:, 0], q[:, 1], s=0.2, c=MID, lw=0, rasterized=True)
    for z in np.linspace(max(h_r - 1.5, 0.1), h_c + 0.8, NV):
        a.axhline(z, c=LIGHT, lw=0.4, zorder=0)
    a.axhline(h_r, c=ROOT, lw=1)
    a.axhline(h_c, c=CUT, lw=1)
    a.set_xlim(*xa)
    a.set_ylim(0, top)
    a.axis("off")
    a.set_title("(a) slices, last root, cut" if rooted else "(a) slices, 1.3 m, cut", fontsize=9)
    b = panel(1, PANEL_H)
    b.plot(n, h, c="k", lw=0.7)
    b.axhline(h_r, c=ROOT, lw=1)
    b.axhline(h_c, c=CUT, lw=1)
    b.text(0.97, h_c + 0.06, "cut", color=CUT, ha="right", va="bottom", transform=b.get_yaxis_transform())
    b.text(0.97, h_r + 0.06, "last root" if rooted else "1.3 m", color=ROOT, ha="right", va="bottom", transform=b.get_yaxis_transform())
    b.set_ylim(0, top)
    b.set_xlim(0, max(float(n[h < top].max()), 3) + 1)
    b.spines[["top", "right"]].set_visible(False)
    b.set_xlabel("blobs per slice" if rooted else "stems tracked")
    b.set_ylabel("height (m)")
    b.set_title("(b)", fontsize=9)
    h_low = max(min(2.0, h_r - 0.5), 0.1)                    # a slice with many blobs
    z_low, z_top = z_ref + h_low, z_cut + 0.05
    k = (np.abs(xyz[:, 2] - z_low) < SLICE / 2) | (np.abs(xyz[:, 2] - z_top) < SLICE / 2)
    lo, hi = xyz[k, :2].min(0), xyz[k, :2].max(0)
    c, half = (lo + hi) / 2, float((hi - lo).max()) / 2 * 1.03
    top_view(panel(2, TOP_W), xyz, z_low, f"(c) slice at {h_low:.1f} m", c, half)
    top_view(panel(3, TOP_W), xyz, z_top, "slice at the cut", c, half)
    d = panel(4, hd)
    trees = sorted((s / "rct" / tree).glob("*_trees.txt"))
    stems = [p for p in trees if "_stem" in p.name]
    roots = [p for p in trees if p.name.endswith("_roots_trees.txt")]
    d.axis("off")
    if stems:                                                # from the cut line up, same scale
        q = sub(np.c_[x[stem] - xs0, hag[stem] - h_c], 60_000)
        d.scatter(q[:, 0], q[:, 1], s=0.2, c=LIGHT, lw=0, rasterized=True)
        for p in stems:
            cylinders(d, rt.read_trees(p), -xs0, lambda z: z - z_cut, y_top)
        d.text(ws / 2, yd[0], "stem", ha="center", va="top", fontsize=9)
        if roots:
            off = ws + gap_d - xr0
            q = sub(np.c_[x[low] + off, h_c - hag[low]], 60_000)
            d.scatter(q[:, 0], q[:, 1], s=0.2, c=LIGHT, lw=0, rasterized=True)
            cylinders(d, rt.read_trees(roots[0]), off, lambda z: z - z_cut, y_top)
            d.text(off + xr0 + wr / 2, yd[0], "roots, flipped", ha="center", va="top", fontsize=9)
        d.axhline(0, c=CUT, lw=1)
        d.set_xlim(*xd)
        d.set_ylim(*yd)
        d.set_title("(d) RCT cylinders from the cut", fontsize=9)
    fig.suptitle("")
    p = s.parent / "qc" / f"{tree}.png"
    p.parent.mkdir(exist_ok=True)
    tmp = p.with_name(p.stem + ".tmp.png")
    fig.savefig(tmp, dpi=110, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    tmp.replace(p)
