"""Grid vs mesh on one cloud: root rows that differ per src and per tree, and the two ground surfaces.
Usage: compare_ground.py <grid roots.npz> <mesh roots.npz> <segmented.ply> <mesh.ply> <report.md>"""
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from rayroots.plot.roots import ground as rw


def rows(p):
    z = np.load(p)
    return pd.DataFrame({"fi": z["fi"], "seg_id": z["seg_id"], "src": z["src"]})


g, m = rows(sys.argv[1]), rows(sys.argv[2])
j = g.merge(m, on=["fi", "seg_id"], how="outer", suffixes=("_grid", "_mesh"), indicator=True)
only_g, only_m = j[j._merge == "left_only"], j[j._merge == "right_only"]
both = j[j._merge == "both"]
changed = both[both.src_grid != both.src_mesh]
out = ["# Grid vs mesh on the crop", "",
       f"Rows (point, segment): grid {len(g):,}, mesh {len(m):,}; in both {len(both):,}, of which src changed "
       f"{len(changed):,}; only grid {len(only_g):,}; only mesh {len(only_m):,}.", ""]
out += ["| src | grid | mesh | only grid | only mesh |", "|---|---:|---:|---:|---:|"]
for c, n in ((1, "v1 only"), (2, "v2 only"), (3, "both")):
    out.append(f"| {n} | {(g.src == c).sum():,} | {(m.src == c).sum():,} | {(only_g.src_grid == c).sum():,} | "
               f"{(only_m.src_mesh == c).sum():,} |")
t = pd.DataFrame({"grid": g.groupby("seg_id").size(), "mesh": m.groupby("seg_id").size(),
                  "only_grid": only_g.groupby("seg_id").size(), "only_mesh": only_m.groupby("seg_id").size()})
t = t.fillna(0).astype(int).sort_values("grid", ascending=False)
t["diff_pct"] = (100 * (t.only_grid + t.only_mesh) / t[["grid", "mesh"]].max(axis=1)).round(1)
out += ["", "Per tree (rows):", "", "| segment | grid | mesh | only grid | only mesh | differ % |",
        "|---:|---:|---:|---:|---:|---:|"]
out += [f"| {s} | {r.grid:,} | {r.mesh:,} | {r.only_grid:,} | {r.only_mesh:,} | {r.diff_pct} |"
        for s, r in zip(t.index, t.itertuples())]
with tempfile.TemporaryDirectory() as d:
    rw.CACHE = Path(d)
    a, b = rw.Ground(Path(sys.argv[3])), rw.MeshGround(Path(sys.argv[4]))
    X, Y = np.meshgrid(a.x0 + (np.arange(a.g.shape[0]) + .5) * rw.GROUND_CELL,
                       a.y0 + (np.arange(a.g.shape[1]) + .5) * rw.GROUND_CELL, indexing="ij")
    dz = (b(X.ravel(), Y.ravel()) - a.g.ravel())
out += ["", f"Ground, mesh minus grid over the grid cells: median {np.median(dz):+.3f} m, "
        f"p5 {np.percentile(dz, 5):+.3f} m, p95 {np.percentile(dz, 95):+.3f} m."]
Path(sys.argv[5]).write_text("\n".join(out) + "\n")
print("\n".join(out[:9]) + "\n...\n" + out[-1])
