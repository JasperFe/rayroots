"""Figures of the README, from the test crop (docs/img/*.png, < 500 kB each).
Usage: make_figures.py <project> <plot> <stem> --hand <hand-cleaned tree.ply> --gallery <gallery screenshot.png> --qc <one.png> <two.png>
<project>: run_plot.sh --project; the crop must be run with --keep. Run in the conda env rayroots. Figure 3 runs RCT."""
# %% Imports
import argparse
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from PIL import Image  # noqa: E402

from rayroots.plot.roots import flip  # noqa: E402
from rayroots.plot.roots import ground as rw  # noqa: E402
from rayroots.plot.roots import plotdata as cm  # noqa: E402
from rayroots.plot.roots import rct_tools as rt  # noqa: E402
from rayroots.plot import page  # noqa: E402
from rayroots.plot.finish import read_ply  # noqa: E402

GREY = "#8c8c8c"
PALE = "#c9c9c9"
CUT = "#c0143c"
SRC = {k: np.array(page.COL[k]) / 255 for k in page.COL}
TREE = "crop_56"                   # hand-root tree SU7_psp3_364_136, one of its two stems
OUT = Path(__file__).parent / "img"
matplotlib.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.titleweight": "normal"})


# %% Helpers
def side(ax, x, z, c=GREY, s=0.15):
    """Orthographic side view."""
    if len(np.shape(c)) == 1 and not isinstance(c, str):
        c = [c]
    ax.scatter(x, z, s=s, c=c, lw=0, rasterized=True)
    ax.set_aspect("equal")
    ax.axis("off")


def save(fig, name):
    p = OUT / name
    fig.savefig(p, dpi=170, facecolor="white", bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    Image.open(p).convert("RGB").quantize(128, dither=Image.Dither.NONE).save(p, optimize=True)
    print(f"{p.name}: {p.stat().st_size / 1000:.0f} kB")


def titles(fig, ax, names):
    """One title height for panels of different shape."""
    fig.canvas.draw()
    top = max(a.get_position().y1 for a in ax)
    for a, name in zip(ax, names):
        p = a.get_position()
        fig.text((p.x0 + p.x1) / 2, top + 0.02, name, ha="center", va="bottom", fontsize=10)


def sample(n, m, seed=0):
    return np.random.default_rng(seed).permutation(n)[:m]


def colours(ids):
    """A fixed colour per id; darker hues only, on white."""
    pal = np.array(plt.get_cmap("tab20")(np.arange(0, 20, 2)))[:, :3]
    return pal[np.asarray(ids) % len(pal)]


def trees_csv(rv):
    return pd.read_csv(rv.parent.parent / "work" / rv.name / "trees.csv")


def tree_file(rv, tree):
    return rv / trees_csv(rv).set_index("tree").at[tree, "file"]


def read_tree(rv, tree):
    a = read_ply(tree_file(rv, tree))
    return {k: np.asarray(a[k], float) for k in ("x", "y", "z", "hag", "root_src")}


def cylinders(path):
    """(p0, p1, r) of every segment of a trees.txt, per tree."""
    seg = []
    for line in open(path):
        if line.startswith("#") or line.startswith("x,"):
            continue
        a = np.array([[float(v) for v in t.strip(",").split(",")] for t in line.split()])
        if a.ndim == 2 and len(a):
            xyz, r, par = a[:, :3], a[:, 3], a[:, 4].astype(int)
            seg += [(xyz[par[i]], xyz[i], r[i]) for i in np.flatnonzero(par >= 0)]
    return seg


def draw_cylinders(ax, seg, u=0, v=2, fz=lambda z: z, color="#5b3a1e"):
    """Each cylinder as its outline in the view plane."""
    q = []
    for p0, p1, r in seg:
        a, b = np.array([p0[u], fz(p0[v])]), np.array([p1[u], fz(p1[v])])
        d = b - a
        n = np.hypot(*d)
        if n == 0:
            continue
        w = np.array([-d[1], d[0]]) / n * r
        q.append([a - w, a + w, b + w, b - w])
    ax.add_collection(PolyCollection(q, facecolors=color, edgecolors="none"))


# %% Figure 1: overview strip
def fig1(project, plot, stem, hand, gallery):
    plots, rv = project / "A_plot" / "work" / plot, project / "A_plot" / "review" / plot
    t = read_tree(rv, TREE)
    y0 = np.median(t["y"][t["hag"] > 2])
    mm, _ = cm.memmap(plots / f"{stem}_segmented_roots.ply")
    a = np.asarray(mm[::12])
    a = a[np.abs(a["y"] - y0) < 3]
    a = a[sample(len(a), 500_000)]
    ids = rw.decode_ids(a["red"], a["green"], a["blue"])
    keep = a["alpha"] > 0
    a, ids = a[keep], ids[keep]
    x, z, src = a["x"], a["z"], a["root_src"]
    h = read_ply(hand)
    hx, hz = np.asarray(h["x"], float), np.asarray(h["z"], float)
    k = sample(len(hx), 200_000)
    tr = trees_csv(rv).sort_values("H_m", ascending=False).head(4)

    img = Image.open(gallery).convert("RGB")
    img = img.crop((0, 0, img.width, int(img.width * 0.72)))
    wid = [1, 1, 1, 1.15, 1.25, 0.42]
    fig, ax = plt.subplots(1, 6, figsize=(13, 3.2), gridspec_kw={"width_ratios": wid, "wspace": 0.05})
    side(ax[0], x, z)
    u = ids >= 0
    side(ax[1], x[~u], z[~u], PALE)
    side(ax[1], x[u & (src == 0)], z[u & (src == 0)], colours(ids[u & (src == 0)]))
    r = src > 0
    side(ax[2], x[~r], z[~r], PALE)
    for c in (1, 2, 3):
        side(ax[2], x[src == c], z[src == c], [SRC[c]], 0.5)
    off = 0.0
    for i, q in enumerate(tr.itertuples()):
        b = read_tree(rv, q.tree)
        s = sample(len(b["x"]), 60_000)
        bx = b["x"][s] - np.median(b["x"]) + off
        side(ax[3], bx, b["z"][s], colours([i])[0])
        off += max(np.ptp(b["x"][s]) * 0.55, 3)
    ax[4].imshow(img)
    ax[4].axis("off")
    side(ax[5], hx[k], hz[k])
    for i in range(3):
        ax[i].set_xlim(x.min(), x.max())
        ax[i].set_ylim(z.min() - 0.5, z.max() + 0.5)
    titles(fig, ax, ["Plot cloud", "RCT segments", "Root proposals", "Single trees", "Gallery, for review", "Cleaned tree"])
    save(fig, "1_overview.png")


# %% Figure 2: RCT segment against root proposals
def fig2(project, plot):
    t = read_tree(project / "A_plot" / "review" / plot, TREE)
    k = t["hag"] < 7
    x, y, z, s = t["x"][k], t["y"][k], t["z"][k], t["root_src"][k].astype(int)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw={"wspace": 0.03})
    r = s == 0
    side(ax[0], x[r], z[r], [SRC[0]], 0.5)
    ax[0].set_title("RCT segment only")
    side(ax[1], x[r], z[r], [SRC[0]], 0.5)
    for c in (1, 2, 3):
        side(ax[1], x[s == c], z[s == c], [SRC[c]], 0.8)
    ax[1].set_title("With root proposals")
    for a in ax:
        a.set_xlim(x.min(), x.max())
        a.set_ylim(z.min(), z.max())
    for c, name in ((3, "v1 and v2"), (1, "v1 only"), (2, "v2 only")):
        ax[1].scatter([], [], s=25, c=[SRC[c]], label=name)
    ax[1].legend(loc="upper right", frameon=False, handletextpad=0.1)
    save(fig, "2_roots.png")


# %% Figure 3: the flip, on the part of one tree below the cut
def fig3(project, plot, hcut=4.0):
    rv = project / "A_plot" / "review" / plot
    t = read_tree(rv, TREE)
    c = trees_csv(rv).set_index("tree").loc[TREE, ["x", "y"]].to_numpy(float)
    x, y, h = t["x"], t["y"], t["hag"]
    low = (h < hcut) & (h > 0.15)
    with tempfile.TemporaryDirectory(dir=project) as d:
        d = Path(d)
        zf = 2 * hcut - h[low]
        flip.write_ply_xyz(d / "w.ply", np.c_[x[low], y[low], zf])
        flip.write_patches(d / "g.ply", c[None], np.array([0.5]), zf.min() + rt.GROUND_OFFSET)
        err = rt.rayimport(d / "w.ply", d / "w_rc.ply") or rt.rayextract_trees(d / "w_rc.ply", d / "g.ply", flip.TREES, str(d / "out"))
        if err:
            raise SystemExit(f"RCT: {err}")
        seg = cylinders(d / "out_trees.txt")
    print(f"figure 3: {low.sum():,} points, {len(seg):,} cylinders")
    keep = h < 7.5
    fig, ax = plt.subplots(1, 4, figsize=(13, 4), gridspec_kw={"wspace": 0.04})
    lo = low & (h >= 0)
    up = keep & (h >= hcut)
    side(ax[0], x[up], h[up], PALE, 0.3)
    side(ax[0], x[lo], h[lo], GREY, 0.3)
    ax[0].axhline(hcut, color=CUT, lw=1)
    ax[0].text(x[keep].min(), hcut + 0.08, "cut, 4 m above ground", color=CUT, va="bottom")
    side(ax[1], x[low], zf, GREY, 0.3)
    ax[1].axhline(hcut, color=CUT, lw=1)
    side(ax[2], x[low], zf, PALE, 0.3)
    draw_cylinders(ax[2], seg)
    ax[2].axhline(hcut, color=CUT, lw=1)
    side(ax[3], x[up], h[up], PALE, 0.3)
    side(ax[3], x[low], h[low], PALE, 0.3)
    draw_cylinders(ax[3], seg, fz=lambda z: 2 * hcut - z)
    ax[3].axhline(hcut, color=CUT, lw=1)
    for a in ax:
        a.set_xlim(x[keep].min(), x[keep].max())
    for a in (ax[0], ax[3]):
        a.set_ylim(0, 7.5)
    for a in ax[1:3]:
        a.set_ylim(hcut - 0.5, 2 * hcut)
    titles(fig, ax, ["Side view", "Part below the cut, flipped", "RCT grows trees from the cut", "Mapped back: root points"])
    save(fig, "3_flip.png")


# %% Figure 4: gallery
def fig4(gallery):
    img = Image.open(gallery).convert("RGB")
    p = OUT / "4_gallery.png"
    img.resize((1100, round(img.height * 1100 / img.width)), Image.LANCZOS).quantize(128).save(p, optimize=True)
    print(f"{p.name}: {p.stat().st_size / 1000:.0f} kB")


# %% Figure 5: the cut, two rayroots-single QC figures
def fig5(qc):
    for src, name in zip(qc, ("5_cut.png", "5_cut_two_stems.png")):
        p = OUT / name
        Image.open(src).convert("RGB").quantize(128, dither=Image.Dither.NONE).save(p, optimize=True)
        print(f"{p.name}: {p.stat().st_size / 1000:.0f} kB")


# %% Run
if __name__ == "__main__":
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("project", type=Path)
    a.add_argument("plot")
    a.add_argument("stem")
    a.add_argument("--hand", type=Path, required=True)
    a.add_argument("--gallery", type=Path, required=True)
    a.add_argument("--qc", type=Path, nargs=2, required=True, help="QC figures: single-stem tree, two-stem tree")
    x = a.parse_args()
    OUT.mkdir(exist_ok=True)
    fig1(x.project, x.plot, x.stem, x.hand, x.gallery)
    fig2(x.project, x.plot)
    fig3(x.project, x.plot)
    fig4(x.gallery)
    fig5(x.qc)
