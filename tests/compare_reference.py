"""A hand-corrected reference tree against the tool's result on the crop. Points within 5 m of the stem and below
4 m of the tree's base, matched within 2 mm. Reference roots = points outside the tree's own RCT segments in the
reference file. Predicted = the crop segments holding the reference stem. Same idea as the study's untouched test,
not the same script; the crop's RCT segments differ from the full plot's.
Usage: compare_reference.py <roots.ply> <roots.npz> <reference tree.ply> <own ids, comma separated> [out.md]"""
import sys

import numpy as np
from scipy.spatial import cKDTree

from rayroots.plot.roots import ground as rw

ply, npz, ref, own = sys.argv[1], sys.argv[2], sys.argv[3], [int(v) for v in sys.argv[4].split(",")]
c = rw.ply_memmap(ply)
xyz = np.c_[c["x"], c["y"], c["z"]].astype(np.float64)
now = rw.decode_ids(c["red"], c["green"], c["blue"])
r = rw.ply_memmap(ref)
rx = np.c_[r["x"], r["y"], r["z"]].astype(np.float64)
rid = rw.decode_ids(r["red"], r["green"], r["blue"])
lo, hi = xyz[:, :2].min(0), xyz[:, :2].max(0)
box = ((rx[:, :2] >= lo) & (rx[:, :2] <= hi)).all(1)
rx, rid = rx[box], rid[box]
d, k = cKDTree(xyz).query(rx, distance_upper_bound=0.002)
ok = np.isfinite(d)
print(f"reference points in the crop {len(rx):,}, matched {ok.mean():.4f}")
rk, is_own = k[ok], np.isin(rid[ok], own)
z0 = rx[ok][:, 2].min()
stem = is_own & (xyz[rk, 2] > z0 + 2) & (xyz[rk, 2] < z0 + 8)
v = rw.decode_ids(c["red"], c["green"], c["blue"])[rk[stem]]
u, n = np.unique(v[v >= 0], return_counts=True)
segs = u[n >= 0.1 * n.sum()]
print("crop segments of the reference stem:", dict(zip(u.tolist(), n.tolist())), "->", segs.tolist())
cen = np.median(rx[ok][is_own & (xyz[rk, 2] > z0 + 3) & (xyz[rk, 2] < z0 + 4)][:, :2], axis=0)
z = np.load(npz)
in_s = np.isin(z["seg_id"], segs)
rct = np.zeros(len(xyz), bool)
rct[np.isin(now, segs) & (c["root_src"] == 0)] = True
sets = {"RCT": rct}
for name, codes in (("v1", (1, 3)), ("v2", (2, 3)), ("union", (1, 2, 3))):
    a = rct.copy()
    a[z["fi"][in_s & np.isin(z["src"], codes)]] = True
    sets[name] = a
dom = (np.linalg.norm(xyz[:, :2] - cen, axis=1) <= 5) & (xyz[:, 2] < z0 + 4)
truth = np.zeros(len(xyz), bool)
truth[rk] = True
lines = ["| method | precision | recall | F1 | TP | FP | FN |", "|---|---:|---:|---:|---:|---:|---:|"]
for name, a in sets.items():
    tp = int((a & truth & dom).sum()); fp = int((a & ~truth & dom).sum()); fn = int((~a & truth & dom).sum())
    pr, re = tp / max(tp + fp, 1), tp / max(tp + fn, 1)
    lines.append(f"| {name} | {pr:.3f} | {re:.3f} | {2 * pr * re / max(pr + re, 1e-9):.3f} | {tp:,} | {fp:,} | {fn:,} |")
print("\n".join(lines))
if len(sys.argv) > 5:
    open(sys.argv[5], "w").write("\n".join(lines) + "\n")
