"""The automatic cut against hand cuts. Only the cut, not the multi-stem split.
Full tree = roots un-flipped (z = 2 z_ref - z) + top; z_ref = the hand cut. Data read only, not in the repo.
Usage: validate_cut.py <folder with tops/ and roots/> <out dir>"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from rayroots.plot.finish import read_ply
from rayroots.single import split as sp

src, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
zdb = pd.read_csv(src / "roots" / "z_refdb_part2.csv")
zdb["tree"] = zdb.file.map(lambda f: Path(f).name.removesuffix("_roots_transformed.ply"))
zdb = zdb[[(src / "tops" / f"{t}_transformed.ply").exists() for t in zdb.tree]].sort_values("tree")


def xyz(path):
    a = read_ply(path)
    return np.c_[a["x"], a["y"], a["z"]].astype(float)


rows = []
for q in zdb.itertuples():
    r = xyz(src / "roots" / f"{q.tree}_roots_transformed.ply")
    r[:, 2] = 2 * q.z_ref - r[:, 2]
    t = np.vstack([r, xyz(src / "tops" / f"{q.tree}_transformed.ply")])
    t = t[sp.largest_component(t, sp.noise_voxel, sp.noise_link)[0]]
    z0 = np.percentile(t[:, 2], sp.z0_pct)
    h, c, cen = sp.slice_profile(t, z0)
    h_root, h_cut, n, run, _ = sp.find_cut(h, c, cen)
    row = {"tree": q.tree, "h_hand": round(q.z_ref - z0, 2), "n_points": len(t), "n_stems": n, "clean_run": round(run, 2)}
    if n:
        row |= {"h_auto": round(h_cut, 2), "diff_m": round(z0 + h_cut - q.z_ref, 2),
                "flag": "short run" if run < sp.check_run else ""}
    else:
        row |= {"flag": "no cut found"}
    rows.append(row)
    print(row, flush=True)
t = pd.DataFrame(rows)
t.to_csv(out / "cut_validation.csv", index=False)
d = t.diff_m.dropna()


def md(t):
    head = "| " + " | ".join(map(str, t.columns)) + " |\n|" + "---|" * len(t.columns) + "\n"
    return head + "".join("| " + " | ".join("" if pd.isna(v) else str(v) for v in r) + " |\n" for r in t.itertuples(index=False))


summary = (f"| trees | {len(t)} |\n|---|---:|\n| cut found | {len(d)} |\n| mean difference | {d.mean():.2f} m |\n"
           f"| mean absolute difference | {d.abs().mean():.2f} m |\n| within 0.5 m | {(d.abs() <= 0.5).sum()}/{len(d)} |\n"
           f"| flagged `short run` | {(t.flag == 'short run').sum()} |\n| no cut found | {(t.flag == 'no cut found').sum()} |\n")
table = md(t[["tree", "h_hand", "h_auto", "diff_m", "n_stems", "clean_run", "flag"]])
(out / "summary.md").write_text(summary + "\n" + table)
print(summary)
