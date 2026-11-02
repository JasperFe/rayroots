"""context: the points of the segmented plot file in a vertical cylinder around a tree's stem position, full height,
without the tree's own points; field seg_id kept (-1 = unsegmented). One pass over the file for all trees. Existing
files are not redone."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from rayroots.plot.roots import plotdata as cm
from rayroots.plot import page, settings

FIELDS = ["x", "y", "z", "alpha", "red", "green", "blue"]


def read_xyz(path):
    """x, y, z (float32) of a binary little-endian PLY written by review."""
    types = {"float": "<f4", "double": "<f8", "uchar": "u1", "int": "<i4"}
    with open(path, "rb") as f:
        props = []
        while (line := f.readline().decode("ascii").strip()) != "end_header":
            if line.startswith("element vertex"):
                n = int(line.split()[2])
            elif line.startswith("property"):
                props.append((line.split()[2], types[line.split()[1]]))
        a = np.fromfile(f, np.dtype(props), count=n)
    return np.c_[a["x"], a["y"], a["z"]].astype(np.float32)


def write_ply(path, xyz, seg_id):
    rec = np.empty(len(xyz), [("x", "<f4"), ("y", "<f4"), ("z", "<f4"), ("red", "u1"), ("green", "u1"),
                              ("blue", "u1"), ("seg_id", "<i4")])
    rec["x"], rec["y"], rec["z"] = xyz.T
    rec["red"], rec["green"], rec["blue"] = page.CTX_COL
    rec["seg_id"] = seg_id
    tmp = path.with_suffix(".part")
    with open(tmp, "wb") as f:
        f.write(("ply\nformat binary_little_endian 1.0\nelement vertex %d\nproperty float x\nproperty float y\n"
                 "property float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\n"
                 "property int seg_id\nend_header\n" % len(rec)).encode())
        rec.tofile(f)
    tmp.replace(path)                                  # complete files only


def run(plot, trees=None, project=".", radius=None, segmented=None, settings_path=None):
    radius = radius if radius is not None else settings.load(settings_path)["context"]["radius_m"]
    rv, work = settings.review_dir(project, plot), settings.work_dir(project, plot)
    if not (work / "trees.csv").exists():
        print(f"no review of {plot} in {rv}")
        return 2
    seg = Path(segmented) if segmented else Path(json.load(open(work / "run.json"))["segmented"])
    if not seg.is_file():
        print(f"not found: {seg} (moved? give --segmented)")
        return 2
    out = rv / "context"
    tr = pd.read_csv(work / "trees.csv").set_index("tree")
    dest = {t: out / f"{t}_context.ply" for t in tr.index}
    if not trees:
        ck = pd.read_csv(rv / "checklist.csv", dtype=str).fillna("")
        trees = ck.loc[ck.status.str.strip().str.lower() == "fix", "tree"].tolist()
    trees = [t if t in tr.index else f"{plot}_{t}" for t in trees]   # segment id -> tree name
    bad = [t for t in trees if t not in tr.index]
    if bad:
        print(f"{plot}: not in trees.csv: {', '.join(bad)}")
        return 2
    trees = [t for t in trees if not dest[t].exists()]
    if not trees:
        print(f"{plot}: no trees without context")
        return 0
    jobs = []
    for t in trees:
        own = read_xyz(rv / tr.at[t, "file"])
        c = np.array([tr.at[t, "x"], tr.at[t, "y"]])
        jobs.append({"tree": t, "c": c, "r": radius, "own": cKDTree(own), "xyz": [], "id": []})
        print(f"{t}: radius {radius:g} m", flush=True)
    for s, ch, ids in cm.chunks(seg, FIELDS):
        k = ch["alpha"] > 0                            # bounded points only (as the stored plot points)
        x, y = ch["x"].astype(np.float64), ch["y"].astype(np.float64)
        for j in jobs:
            m = k & (np.hypot(x - j["c"][0], y - j["c"][1]) <= j["r"])
            j["xyz"].append(np.c_[ch["x"][m], ch["y"][m], ch["z"][m]].astype(np.float32))
            j["id"].append(ids[m].astype(np.int32))
        print(f"  {plot}: {s + len(ch):,}", flush=True)
    for j in jobs:
        xyz, sid = np.concatenate(j["xyz"]), np.concatenate(j["id"])
        d = j["own"].query(xyz, distance_upper_bound=1e-6)[0]
        keep = ~np.isfinite(d)                         # not one of the tree's own points (exact coordinates)
        dest[j["tree"]].parent.mkdir(parents=True, exist_ok=True)
        write_ply(dest[j["tree"]], xyz[keep], sid[keep])
        print(f"{j['tree']}: context {keep.sum():,} points ({(~keep).sum():,} own excluded), "
              f"{dest[j['tree']].stat().st_size / 1e6:.0f} MB -> {dest[j['tree']]}", flush=True)
    return 0
