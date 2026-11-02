"""Regression: rayroots on a full segmented file against the frozen roots.npz (fi, seg_id, src: array_equal).
Usage: regression.py <segmented.ply> <reference roots.npz> <boundary.csv> <out dir> [--known id id ...]
--known: test only, shows whether a difference comes from the known trees the tool leaves out."""
import argparse
import sys
from pathlib import Path

import numpy as np

from rayroots.plot.roots import run as cli
from rayroots.plot.roots import plotdata as cm

a = argparse.ArgumentParser()
a.add_argument("seg", type=Path); a.add_argument("ref", type=Path); a.add_argument("boundary", type=Path)
a.add_argument("out", type=Path); a.add_argument("--known", type=int, nargs="*", default=[])
x = a.parse_args()
if x.known:
    cm.known_ids = lambda plot: set(x.known)
rc = cli.run(x.seg, None, x.boundary, x.out)
if rc:
    sys.exit(rc)
stem = x.seg.stem.removesuffix("_segmented")
n, r = np.load(x.out / f"{stem}_roots.npz"), np.load(x.ref)
ok = True
for k in ("fi", "seg_id", "src"):
    same = np.array_equal(n[k], r[k])
    ok &= same
    print(f"{k}: array_equal {same}  (new {len(n[k]):,}, reference {len(r[k]):,})")
if not ok:
    a_, b_ = (set(zip(z["fi"].tolist(), z["seg_id"].tolist(), z["src"].tolist())) for z in (n, r))
    print(f"rows only in new {len(a_ - b_):,}, only in reference {len(b_ - a_):,}")
print("REGRESSION", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
