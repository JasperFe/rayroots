"""rayroots-single against logged reference results (read only).
Usage: compare_single.py <project> <reference tree_volumes.csv> <reference pipeline folder> <report.md>"""
import sys
from pathlib import Path

import pandas as pd

project, ref_v, ref_p, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3]), sys.argv[4]
log = pd.read_csv(project / "B_trees" / "work" / "split_log.csv").set_index("tree")
new = pd.read_csv(project / "B_trees" / "tree_volumes.csv")
ref = pd.read_csv(ref_v)
lines, bad = ["# rayroots-single against the reference logs", ""], 0


def check(tree, name, a, b, tol=0):
    global bad
    same = str(a) == str(b) or (pd.api.types.is_number(a) and pd.api.types.is_number(b) and abs(float(a) - float(b)) <= tol)
    bad += not same
    lines.append(f"| {tree} | {name} | {a} | {b} | {'same' if same else 'DIFF'} |")


lines += ["| tree | field | new | reference | |", "|---|---|---:|---:|---|"]
for tree, sp in (("SU7_psp3_364_136", "rhizophora"), ("SU7_psp3_119", "avicennia")):
    r2 = pd.read_csv(ref_p / sp / "02_split" / "log_02_split.csv").set_index("tree_id").loc[tree]
    for c in ("z_ref", "z_cut", "n_stems", "n_top", "points_per_stem"):
        check(tree, c, log.at[tree, c], r2[c], 1e-9)
    n, v = new[new.tree == tree], ref[ref.tree_id == tree]
    a, b = n[n["row"] == "tree"].iloc[0], v[v["row"] == "tree"].iloc[0]
    col = "_central" if "vol_rct_central_m3" in n and pd.notna(a.get("vol_rct_central_m3")) else ""
    for c in ("vol_rct", "vol_meas"):
        check(tree, f"{c} central", a[f"{c}{col}_m3"], b[f"{c}_central_m3"], 1e-9)
    if col:
        for c in ("vol_rct_m3", "vol_meas_m3", "vol_roots_rct_m3", "vol_roots_meas_m3"):
            check(tree, f"{c}, median of 15 runs", a[c], b[c], 1e-9)
lines += ["", "**PASS**" if not bad else f"**{bad} DIFFER**", ""]
Path(out).write_text("\n".join(lines))
print("\n".join(lines))
