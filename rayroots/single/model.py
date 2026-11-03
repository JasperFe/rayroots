"""rayroots-single model: RCT per stem, base and root system, volumes and DBH -> B_trees/tree_volumes.csv.
Roots are flipped about the cut (z' = 2 z_cut - z), so RCT grows them like a tree. Default: the central run only, one
result per tree. --ensemble adds 15 runs (90 % subsamples x cut +-0.10 m); the volumes are then their median.
vol_rct: cylinders with RCT's radii. vol_meas: RCT's skeleton with radii measured from the points.
"""
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd

from rayroots.single import figure, settings
from rayroots.single import rct_tools as rt

STEM, ROOT, RAY_VECTOR = None, None, None    # the rct section of single.yaml, set by run()
FLAT_STEP, FLAT_MARGIN = 0.1, 0.5      # m, flat terrain: grid spacing, margin around the cloud
BASE_MARGIN = 0.15                     # m, flag "base off stem"
BAND_LINK = 0.3                        # m, cut-shift band: a point this close to a stem belongs to it
LEAF_CAP = {"stem": True, "base": True, "roots": False}
DBH_RULE = {"dbh": ("ground", 1.3), "roots": ("last root", 0.30)}   # field height
DBH_H_MIN = 0.10                       # m above the measurement start
DBH_WINDOW = 0.30                      # m, a clean height is searched only within +- this of the field height
DBH_MAX_REL_RMSE = 0.05                # fit flagged when rmse > this x fitted radius
MIN_RUNS = 10                          # of the 15 ensemble runs a value needs
FILE = "tree_volumes.csv"


def flip(xyz, z_plane):
    return np.c_[xyz[:, :2], 2 * z_plane - xyz[:, 2]]


def run_rct(xyz, d, name, trees_cfg, z_floor=None):
    """rayimport, flat terrain at the cloud's lowest point (or z_floor), rayextract trees -> d/<name>_trees.txt."""
    ply, rc, terrain = d / f"{name}.ply", d / f"{name}_raycloud.ply", d / f"{name}_flat_terrain.ply"
    rt.write_xyz(ply, xyz)
    err = rt.rayimport(ply, rc, RAY_VECTOR)
    ply.unlink()
    if err:
        return f"rayimport: {err}"
    rt.write_flat_ground(terrain, xyz, FLAT_STEP, FLAT_MARGIN, z_floor)
    err = rt.rayextract_trees(rc, terrain, trees_cfg, d / name, keep=("trees.txt",))
    rc.unlink()
    terrain.unlink()
    return err


def measure(trees_file, xyz, leaf):
    trees = rt.read_trees(trees_file)
    table, radii = rt.measured_radii(trees, xyz, leaf_cap=leaf)
    return {"n_trees": len(trees), "vol_rct": sum(t["volume"] for t in trees), "vol_meas": rt.measured_volume(table)}, trees, table


def model_part(part, name, xyz, d, shifted, flipped_at=None, z_floor=None):
    """Central run, and the ensemble when shifted is given (shift -> cloud). Returns central values, trees, runs."""
    trees_cfg = ROOT if part == "roots" else STEM
    frame = (lambda p, shift=0: flip(p, flipped_at + shift)) if flipped_at is not None else (lambda p, shift=0: p)
    err = run_rct(frame(xyz), d, name, trees_cfg, z_floor)
    if err:
        raise RuntimeError(f"{name}: {err}")
    central, trees, table = measure(d / f"{name}_trees.txt", frame(xyz), LEAF_CAP[part])
    table.to_csv(d / f"{name}_segments.csv", index=False)
    runs = []
    if shifted:
        e = d / "ensemble"
        e.mkdir(exist_ok=True)
        for shift, seed in rt.ENS_RUNS:
            run = rt.ens_name(shift, seed)
            cloud = shifted[shift]
            sub = cloud[np.random.default_rng(seed).random(len(cloud)) < rt.ENS_FRAC]
            err = run_rct(frame(sub, shift), e, f"{name}_{run}", trees_cfg, z_floor)
            if err:
                print(f"  {name} {run}: {err}", flush=True)
                continue
            m = measure(e / f"{name}_{run}_trees.txt", frame(sub, shift), LEAF_CAP[part])[0]
            runs.append({"run": run, "cut_shift_m": shift, "seed": seed, "n_points": len(sub), "rct": m["vol_rct"], "meas": m["vol_meas"]})
    return central, trees, pd.DataFrame(runs)


def stats(prefix, unit, central, values, nd):
    """<prefix>_<unit>: the central value, or the ensemble median when there are enough runs."""
    v = pd.Series(values, dtype=float).dropna()
    out = {f"{prefix}_{unit}": round(float(central), nd)}
    if len(v) >= MIN_RUNS:
        out |= {f"{prefix}_{unit}": round(v.median(), nd), f"{prefix}_central_{unit}": round(float(central), nd),
                f"{prefix}_mean_{unit}": round(v.mean(), nd), f"{prefix}_sd_{unit}": round(v.std(ddof=1), nd)}
    return out


def model_diameters(part, z, xy, ctx):
    """D (cm) at height z near xy on a part's model, with RCT's and measured radii; each ensemble run when there are any."""
    d, name, z_cut = ctx["dir"], ctx["name"], ctx["z_cut"]
    zf = (lambda zc: 2 * zc - z) if part == "roots" else (lambda zc: z)
    trees = rt.read_trees(d / f"{name}_trees.txt")
    radii = rt.radii_from_table(trees, pd.read_csv(d / f"{name}_segments.csv"))
    central = {"rct": rt.diameter_at(trees, [t["r"] for t in trees], zf(z_cut), xy) * 100,
               "meas": rt.diameter_at(trees, radii, zf(z_cut), xy) * 100}
    vals = {"rct": [], "meas": []}
    for _, r in ctx["runs"].iterrows():
        shift, seed = r["cut_shift_m"], int(r["seed"])
        cloud = ctx["shifted"][shift]
        sub = cloud[np.random.default_rng(seed).random(len(cloud)) < rt.ENS_FRAC]
        t = rt.read_trees(d / "ensemble" / f"{name}_{r['run']}_trees.txt")
        zc = z_cut + shift if part == "roots" else z_cut
        _, rr = rt.measured_radii(t, flip(sub, zc) if part == "roots" else sub, leaf_cap=LEAF_CAP[part])
        vals["rct"].append(rt.diameter_at(t, [x["r"] for x in t], zf(zc), xy) * 100)
        vals["meas"].append(rt.diameter_at(t, rr, zf(zc), xy) * 100)
    return stats("DBH_rct", "cm", central["rct"], vals["rct"], 1) | stats("DBH_meas", "cm", central["meas"], vals["meas"], 1)


def dbh(split, r2, z_ref, tree_xyz, stem_xyz, z_cut, ctx_stem, ctx_roots):
    """DBH of one stem at the rule's field height, or the nearest clean height. Returns columns and flags."""
    start, h_field = DBH_RULE[split]
    if start == "ground":
        pts, z0, start_xy = stem_xyz, z_cut, None
        h_target = z_ref + h_field - z_cut
    else:
        pts, z0 = tree_xyz, float(r2["z0"]) + float(r2["h_last_root"])
        start_xy, h_target = rt.stem_at_cut(stem_xyz, z_cut)[0], h_field
    out = {"DBH_h_field_m": round(z0 + h_target - z_ref, 2)}
    if h_target + DBH_WINDOW < DBH_H_MIN:
        return out, ["DBH: field height below the stem's start, not measured"]
    m = rt.clean_diameter(pts, z0, h_target, DBH_H_MIN, start_xy, window=DBH_WINDOW, max_rel_rmse=DBH_MAX_REL_RMSE)
    if m is None:
        return out, [f"DBH: no clean height within {DBH_WINDOW:.2f} m, not measured"]
    z = z0 + m["h"]
    out |= {"DBH_h_m": round(z - z_ref, 2), "DBH_points_cm": m["D_fit_cm"], "DBH_points_rmse_mm": m["fit_rmse_mm"]}
    moved = m["h"] - h_target
    flags = [f"DBH moved {'down' if moved < 0 else 'up'} {abs(moved):.2f} m (junction)"] if abs(moved) > 1e-6 else []
    if m["note"]:
        flags.append(f"DBH points: {m['note']}")
    part = "stem" if z >= z_cut else "roots" if ctx_roots is not None else None
    if part is None:
        flags.append("DBH: below the cut, no roots model")
    else:
        out |= model_diameters(part, z, m["xy"], ctx_stem if part == "stem" else ctx_roots)
    return out, flags


def model_tree(tree, r2, rule, s, ensemble, gap_max):
    """All parts of one tree. Returns the tree row and one row per stem."""
    d = s / "rct" / tree
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    z_ref, split = float(r2["z_ref"]), r2["split"]
    files = sorted((s / "stems").glob(f"{tree}_stem*.ply"), key=lambda p: int(p.stem.rsplit("_stem", 1)[1]))
    stems_xyz = [rt.read_xyz(p) for p in files]
    tree_xyz = rt.read_xyz(s / "standard" / f"{tree}.ply")
    per = r2.get("z_cut_per_stem")
    cuts = [float(v) for v in str(per).split("/")] if isinstance(per, str) else [float(r2["z_cut"])] * len(files)
    shifts = rt.ENS_SHIFTS if ensemble else []
    flags, parts, tables = [], [], []

    roots_file = s / "roots" / f"{tree}_roots.ply"
    ctx_roots = None
    if rule["roots_model"] and roots_file.exists():
        roots = rt.read_xyz(roots_file)
        z_cut = float(r2["z_cut"])
        shifted = {x: rt.shift_cut_roots(roots, stems_xyz, z_cut, x, BAND_LINK) for x in shifts}
        c, _, runs = model_part("roots", f"{tree}_roots", roots, d, shifted, flipped_at=z_cut)
        parts.append(("roots", c))
        tables.append(runs)
        if c["n_trees"] > len(files):
            flags.append("rct split")
        if c["n_trees"] < len(files):
            flags.append("stems merged")
        ctx_roots = {"dir": d, "name": f"{tree}_roots", "z_cut": z_cut, "runs": runs, "shifted": shifted}

    base_file = s / "base" / f"{tree}_base.ply"
    if base_file.exists():
        base = rt.read_xyz(base_file)
        shifted = {x: rt.shift_cut_roots(base, stems_xyz, cuts, x, BAND_LINK) for x in shifts}
        c, _, runs = model_part("base", f"{tree}_base", base, d, shifted, z_floor=z_ref)
        parts.append(("base", c))
        tables.append(runs)

    stem_rows = []
    for i, (p, xyz) in enumerate(zip(files, stems_xyz)):
        name, z_cut = p.stem, cuts[i]
        shifted = {x: rt.shift_cut(xyz, tree_xyz, stems_xyz, i, z_cut, x, BAND_LINK, sib_cuts=cuts) for x in shifts}
        c, trees, runs = model_part("stem", name, xyz, d, shifted)
        parts.append(("stem", c))
        tables.append(runs)
        fl = []
        if not rt.base_check(trees[0], xyz, z_cut, BASE_MARGIN)[2]:
            fl.append("base off stem")
        if c["n_trees"] > 1:
            fl.append("rct split")
        ctx = {"dir": d, "name": name, "z_cut": z_cut, "runs": runs, "shifted": shifted}
        dd, fd = dbh(split, r2, z_ref, tree_xyz, xyz, z_cut, ctx, ctx_roots)
        fl += fd
        flags += [f"stem{i + 1} {f}" for f in fl]
        stem_rows.append({"row": "stem", "tree": tree, "stem": i + 1, "genus": r2["genus"], "n_stems": len(files), "z_ref": z_ref,
                          "h_cut_m": round(z_cut - z_ref, 2), "height_m": round(float(xyz[:, 2].max() - z_ref), 2)}
                         | stats("vol_rct", "m3", c["vol_rct"], runs["rct"] if len(runs) else [], 4)
                         | stats("vol_meas", "m3", c["vol_meas"], runs["meas"] if len(runs) else [], 4)
                         | dd | {"n_runs": max(len(runs), 1), "flags": ";".join(fl)})

    tree_row = {"row": "tree", "tree": tree, "stem": np.nan, "genus": r2["genus"], "n_stems": len(files), "z_ref": z_ref,
                "h_cut_m": round(float(r2["z_cut"]) - z_ref, 2), "height_m": round(float(tree_xyz[:, 2].max() - z_ref), 2)}
    keep = [t for t in tables if len(t)]
    runs = None
    if keep:                                 # run by run: a band the stems gain is the band the roots lose
        ids = set.intersection(*(set(t.run) for t in keep))
        runs = {m: sum(t.set_index("run").loc[sorted(ids), m] for t in keep) for m in ("rct", "meas")}
        if len(ids) < MIN_RUNS:
            flags.append(f"only {len(ids)} ensemble runs in every part")
    for m in ("rct", "meas"):
        tree_row |= stats(f"vol_{m}", "m3", sum(c[f"vol_{m}"] for _, c in parts), runs[m] if runs else [], 4)
        for p, col in (("stem", "stems"), ("base", "base"), ("roots", "roots")):
            sel = [(q, c, t) for (q, c), t in zip(parts, tables) if q == p]
            if sel:
                r = sum(t.set_index("run")[m] for _, _, t in sel) if runs else []
                tree_row |= stats(f"vol_{col}_{m}", "m3", sum(c[f"vol_{m}"] for _, c, _ in sel), r, 4)
    ratio = tree_row["vol_rct_m3"] / tree_row["vol_meas_m3"] if tree_row["vol_meas_m3"] else np.inf
    if abs(ratio - 1) > gap_max:
        flags.append("vol gap")
    get = lambda c: [x.get(c, np.nan) for x in stem_rows]
    tree_row |= {"DBH_h_m": "/".join(f"{v:.2f}" for v in get("DBH_h_m")), "DBH_points_cm": "/".join(f"{v:.1f}" for v in get("DBH_points_cm")),
                 "DBH_rct_cm": "/".join(f"{v:.1f}" for v in get("DBH_rct_cm")), "DBH_meas_cm": "/".join(f"{v:.1f}" for v in get("DBH_meas_cm")),
                 "n_runs": len(ids) if keep else 1, "flags": ";".join(([r2["flags"]] if isinstance(r2["flags"], str) and r2["flags"] else []) + flags)}
    return tree_row, stem_rows


def run(project=".", tree=None, ensemble=False, settings_path=None, redo=False, keep=False, src=None, out=None):
    """Trees split ok and not in tree_volumes.csv; --tree T redoes T, --redo all. Trees to look at are listed at the end.
    The ensemble's tree files are deleted unless --keep. src, out: the --input route (genus from split_log.csv)."""
    global STEM, ROOT, RAY_VECTOR
    cfg = settings.load(settings_path)
    if src is None:
        part, todo = settings.part(project), settings.trees(project, cfg, tree)
        again = f"--project {project}"
    else:
        part, todo = settings.direct(cfg, src, None, out, tree, need_genus=False)
        again = f"--input {src}" + (f" --out {out}" if out else "")
    s = part / "work"
    log_p = s / "split_log.csv"
    if not log_p.exists():
        raise SystemExit(f"not found: {log_p}; run rayroots-single split first")
    STEM, ROOT, RAY_VECTOR = cfg["rct"]["stem"], cfg["rct"]["roots"], cfg["rct"]["ray_vector"]
    log = pd.read_csv(log_p).set_index("tree")
    out_p = part / FILE
    old = pd.read_csv(out_p) if out_p.exists() else pd.DataFrame(columns=["row", "tree", "flags"])
    names = [t.tree for t in todo]
    names = [t for t in names if t in log.index and log.at[t, "status"] == "ok"]
    done = set(old.tree)
    rows, failed = [], []
    for t in names:
        if t in done and tree is None and not redo:
            continue
        t0 = time.time()
        print(t, "...", flush=True)
        try:
            r, st = model_tree(t, log.loc[t], cfg["genus"][log.at[t, "genus"]], s, ensemble, cfg["gap_max"])
            rows += [r] + st
            figure.qc(s, log.loc[t].to_dict() | {"tree": t})
            if not keep:
                shutil.rmtree(s / "rct" / t / "ensemble", ignore_errors=True)
            print(f"  vol_rct {r['vol_rct_m3']}, vol_meas {r['vol_meas_m3']} m3, {time.time() - t0:.0f} s {r['flags']}", flush=True)
        except Exception as e:
            failed.append(f"{t}: {e!r}")
            print(f"  FAILED: {e!r}", flush=True)
    new = pd.DataFrame(rows)
    if len(new):
        old = old[~old.tree.isin(new.tree)].dropna(axis=1, how="all")
        tab = pd.concat([old, new], ignore_index=True).sort_values(["tree", "row", "stem"], ascending=[True, False, True])
        first = ["row", "tree", "stem", "genus", "n_stems", "height_m", "z_ref", "h_cut_m", "vol_rct_m3", "vol_meas_m3", "DBH_h_field_m",
                 "DBH_h_m", "DBH_points_cm", "DBH_points_rmse_mm", "DBH_rct_cm", "DBH_meas_cm"]
        tab = tab[[c for c in first if c in tab] + [c for c in tab if c not in first and c != "flags"] + ["flags"]]
        tab.to_csv(out_p, index=False)
        old = tab
    n_new = int((new["row"] == "tree").sum()) if len(new) else 0
    print(f"{n_new} trees modelled; {out_p}")
    see = old[(old["row"] == "tree") & old["flags"].fillna("").str.contains("rct split|stems merged|vol gap") & (old["n_runs"] == 1)]
    if len(see):
        print("Central run may be off (rct split, stems merged or vol gap); rerun these:")
        for q in see.itertuples():
            print(f"  rayroots-single model --tree {q.tree} --ensemble {again}    # {q.flags}")
    for f in failed:
        print("FAILED", f)
    return 1 if failed else 0
