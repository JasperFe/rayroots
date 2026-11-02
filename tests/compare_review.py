"""Regression of select + review against a reference review folder.
Usage: compare_review.py <A_plot/work/<plot>> <new trees_in_plot.csv> <reference review dir> <report.md> [--known id ...]
--known: segments the reference left out as known trees (test only). Without roots only tree set and point counts."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def ply_n(path):
    with open(path, "rb") as f:
        while (s := f.readline().decode("ascii", "replace").strip()) != "end_header":
            if s.startswith("element vertex"):
                return int(s.split()[2])


new, sel, ref, rep = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
known = {int(v) for v in sys.argv[sys.argv.index("--known") + 1:]} if "--known" in sys.argv else set()
n, r, s = pd.read_csv(new / "trees.csv"), pd.read_csv(ref / "trees.csv"), pd.read_csv(sel)
below, edge = pd.read_csv(new / "trees_below_3m.csv"), pd.read_csv(new / "trees_edge_truncated.csv")
out = [f"# Review regression: {new.name}", "", f"Reference: `{ref}`", ""]
rs = set(r.seg_id)
ns, bs, es, ks = set(n.seg_id), set(below.seg_id), set(edge.seg_id), known & set(s.seg_id)
out += [f"- new review {len(ns)}, below 3 m {len(bs)}, edge-truncated {len(es)}; reference {len(rs)}",
        f"- reference trees in the new review {len(rs & ns)}, in trees_below_3m {len(rs & bs)}, "
        f"edge-truncated {len(rs & es)}, missing everywhere {sorted(rs - ns - bs - es)}",
        f"- new review trees not in the reference: {sorted(ns - rs)}; of which known (left out there) "
        f"{sorted((ns - rs) & known)}; other {sorted(ns - rs - known)}",
        f"- known segments edge-truncated or below 3 m here: {sorted(ks - ns)}"]
same_set = (rs - bs - es) == (ns - ks)
j = n.merge(r, on="seg_id", suffixes=("", "_r"))
rooted = "n_v1_only" in r and r.n_v1_only.fillna(0).sum() > 0
cols = ["n_rct", "H_m", "D_cm"] + (["n_v1_only", "n_v2_only", "n_both"] if rooted else [])
rows = ["| field | trees compared | differ | max abs diff |", "|---|---:|---:|---:|"]
ok = same_set
for c in cols:
    a, b = j[c].to_numpy(float), j[f"{c}_r"].to_numpy(float)
    both = np.isfinite(a) & np.isfinite(b)
    tol = 0.0011 if c in ("H_m", "D_cm") else 0
    d = np.abs(a - b)[both]
    k = int((d > tol).sum())
    ok &= k == 0
    rows.append(f"| {c} | {int(both.sum())} | {k} | {d.max() if len(d) else 0:g} |")
    if (~both).any():
        x = j.loc[~both, ["tree", c, f"{c}_r"]]
        out.append(f"- {c}: value on one side only: {x.to_dict('records')}")
out += ["", *rows, ""]
files = {p.name.split("_", 1)[1][:-4]: p for p in (new.parent.parent / "review" / new.name / "trees").glob("*.ply")}
rf = {p.name.split("_", 1)[1][:-4]: p for p in (ref / "trees").glob("*.ply")}
nf = {t.tree: t.seg_id for t in n.itertuples()}
diff = [(t, ply_n(files[t]), ply_n(rf[t])) for t in sorted(set(files) & set(rf)) if ply_n(files[t]) != ply_n(rf[t])]
ok &= not diff
out.append(f"Tree files, vertex counts: {len(set(files) & set(rf))} compared, {len(diff)} differ {diff[:10]}")
if rooted:
    g = j[["tree", "genus", "priority", "genus_r", "priority_r"]].rename(columns={"genus": "genus (new)", "priority": "priority (new)", "genus_r": "genus (reference)", "priority_r": "priority (reference)"})
    out += ["", "## Genus and priority side by side (not compared)", "", "| tree | genus new | priority new | genus reference | priority reference |", "|---|---|---:|---|---:|"]
    out += [f"| {a} | {b} | {c} | {d} | {e} |" for a, b, c, d, e in g.sort_values("priority (reference)").itertuples(index=False)]
out += ["", f"**{'PASS' if ok else 'DIFFERENCES'}**"]
rep.write_text("\n".join(out) + "\n")
print("\n".join(l for l in out if "|" not in l or "field" in l or l.startswith("| n_") or l.startswith("| H") or l.startswith("| D")))
sys.exit(0 if ok else 1)
