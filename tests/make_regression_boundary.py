"""Boundary of the frozen SU7_psp3 run: scan positions (RiSCAN csv) minus the ROUNDED rayimport offset, as the
roots.npz of 2026-10-21 used. Only for the regression; new runs use the exact LAS header offset.
Usage: make_regression_boundary.py <scan positions csv> <offset x> <offset y> <out csv>"""
import sys

import pandas as pd

src, ox, oy, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
sp = pd.read_csv(src, skipinitialspace=True)
sp.columns = [c.lstrip("#") for c in sp.columns]
sp = sp[(sp[["X", "Y", "Z"]] != 0).any(axis=1)]
pd.DataFrame({"x": sp.X - ox, "y": sp.Y - oy}).to_csv(out, index=False, float_format="%.17g")
print(len(sp), "positions ->", out)
