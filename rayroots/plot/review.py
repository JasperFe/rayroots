"""review: the review folder of one plot, with or without root proposals. Genus comes from hints instead of
an inventory; priority is size-class coverage."""
import json
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd

from rayroots.plot.roots import run as roots_run
from rayroots.plot.roots import config, output
from rayroots.plot.roots import plotdata as cm
from rayroots.plot import page, settings
from rayroots.plot import store as sto


def root_volume(z, xyz, hag, h_cut, vox):
    """Per segment: volume (m3) of the occupied voxels of the points v1 or v2 attach to it below h_cut."""
    low = hag < h_cut
    v = np.floor(xyz[low] / vox).astype(np.int64)
    v -= v.min(0)
    key = (v[:, 0] * (v[:, 1].max() + 1) + v[:, 1]) * (v[:, 2].max() + 1) + v[:, 2]
    u = np.unique(np.c_[z["seg_id"][low].astype(np.int64), key], axis=0)
    s, n = np.unique(u[:, 0], return_counts=True)
    return {int(a): float(b) * vox ** 3 for a, b in zip(s, n)}


def read_roots(plot, npz):
    """<stem>_roots.npz with the points' xyz (float32, from the segmented file) and hag (float32)."""
    z = dict(np.load(npz))
    seg = cm.seg_path(plot)
    if str(z["segmented_file"]) != seg.name or int(z["segmented_bytes"]) != seg.stat().st_size:
        raise SystemExit(f"{plot}: {npz} was made from {z['segmented_file']}, not {seg.name}: rerun rayroots")
    mm, _ = cm.memmap(seg)
    xyz = np.empty((len(z["fi"]), 3), np.float32)
    for a in range(0, len(z["fi"]), cm.CHUNK):
        c = mm[z["fi"][a:a + cm.CHUNK]]
        xyz[a:a + cm.CHUNK] = np.c_[c["x"], c["y"], c["z"]]
    g = cm.ground(plot)
    hag = (xyz[:, 2] - g(xyz[:, 0].astype(np.float64), xyz[:, 1].astype(np.float64))).astype(np.float32)
    return z, xyz, hag


def genus_table(t, hints, plot_genus):
    """genus and genus_source per tree: genus_hints.csv (column tree or seg_id, and genus), else plot_genus."""
    t["genus"], t["genus_source"] = "unknown", "no hint"
    if plot_genus:
        t["genus"], t["genus_source"] = plot_genus, "plot_genus"
    if hints is not None:
        h = pd.read_csv(hints, dtype=str).fillna("")
        key = "tree" if "tree" in h.columns else "seg_id"
        m = h.set_index(key)["genus"]
        v = (t.tree if key == "tree" else t.seg_id.astype(str)).map(m).fillna("").astype(str)
        t.loc[v != "", "genus"] = v[v != ""]
        t.loc[v != "", "genus_source"] = "genus_hints"
    return t


def read_existing(path):
    """existing_trees.csv: genus, H_m, D_cm of trees that count as already covered."""
    if path is None:
        return pd.DataFrame({"genus": pd.Series(dtype=str), "H_m": pd.Series(dtype=float),
                             "D_cm": pd.Series(dtype=float)})
    e = pd.read_csv(path)
    bad = {"genus", "H_m", "D_cm"} - set(e.columns)
    if bad:
        raise SystemExit(f"{path}: missing columns {sorted(bad)}")
    return e


def stratum(d, h, d_bins, h_bins):
    return (int(np.digitize(d, d_bins)) if np.isfinite(d) else -1, int(np.digitize(h, h_bins)))


def priority(cand, existing, d_bins, h_bins):
    """Greedy per genus: next = the candidate whose D x H class holds the fewest trees (existing + already ranked);
    ties: larger D first. Trees flagged 'low' (implausible or missing D) come after all others of their genus."""
    cand = cand.copy()
    cand["priority"] = 0
    for genus, g in cand.groupby("genus"):
        ex = existing if genus == "unknown" else existing[existing.genus == genus]
        cnt = {}
        for r in ex.itertuples():
            k = stratum(r.D_cm, r.H_m, d_bins, h_bins)
            cnt[k] = cnt.get(k, 0) + 1
        rank = 1
        for left in (g.index[~g.low].tolist(), g.index[g.low].tolist()):
            while left:
                key = [(cnt.get(stratum(cand.at[i, "D_cm"], cand.at[i, "H_m"], d_bins, h_bins), 0),
                        -np.nan_to_num(cand.at[i, "D_cm"], nan=-1)) for i in left]
                i = left[min(range(len(left)), key=key.__getitem__)]
                cand.at[i, "priority"] = rank
                k = stratum(cand.at[i, "D_cm"], cand.at[i, "H_m"], d_bins, h_bins)
                cnt[k] = cnt.get(k, 0) + 1
                left.remove(i)
                rank += 1
    return cand


def implausible(cand, ex, d_max):
    """Flag text per tree: D above d_max, missing D, or D/H outside the range of the existing trees when given."""
    e = ex.dropna(subset=["D_cm", "H_m"])
    ratio = cand.D_cm / cand.H_m
    high = cand.D_cm > d_max
    out_range = pd.Series(False, index=cand.index)
    txt = ""
    if len(e):
        rmin, rmax = (e.D_cm / e.H_m).min(), (e.D_cm / e.H_m).max()
        out_range = np.isfinite(ratio) & ((ratio < rmin) | (ratio > rmax))
        txt = f"implausible D for H (D/H outside {rmin:.2f}-{rmax:.2f} cm/m of existing trees)"
    flags = []
    for i in cand.index:
        f = []
        if high[i]:
            f.append(f"implausible D (above {d_max:g} cm)")
        if out_range[i]:
            f.append(txt)
        if not np.isfinite(cand.at[i, "D_cm"]):
            f.append("no D")
        flags.append("; ".join(f))
    return flags, (high | out_range).to_numpy()


def suggest_genus(cand, impl, p):
    """Genus 'unknown' -> from the root index (root volume / D^2), or root volume alone without a usable D.
    Threshold calibrated on Suriname data."""
    label = "threshold calibrated on Suriname data"
    ri = cand.root_vol_m3 / (cand.D_cm / 100) ** 2
    for k, i in enumerate(cand.index):
        if cand.at[i, "genus"] != "unknown":
            continue
        if np.isfinite(ri[i]) and not impl[k]:
            r_ = ri[i] >= p["genus_index_threshold_m3_m2"]
            src = f"suggested (root volume / D2 {ri[i]:.1f}, {label})"
        else:
            r_ = cand.at[i, "root_vol_m3"] >= p["genus_volume_threshold_m3"]
            src = f"suggested (root volume {cand.at[i, 'root_vol_m3']:.2f} m3, {label})"
        cand.at[i, "genus"], cand.at[i, "genus_source"] = ("Rhizophora" if r_ else "Avicennia"), src
    return cand


def write_ply(path, rec):
    t = {"f4": "float", "u1": "uchar", "i4": "int", "f8": "double"}
    head = ["ply", "format binary_little_endian 1.0", f"element vertex {len(rec)}"]
    head += [f"property {t[rec.dtype[n].str[1:]]} {n}" for n in rec.dtype.names] + ["end_header"]
    with open(path, "wb") as f:
        f.write(("\n".join(head) + "\n").encode())
        rec.tofile(f)


def write_trees(plot, cand, paths, roots, tmp, meta, p):
    """Tree files and previews. With roots: own segment + all roots.npz rows of the segment (a contested point is in
    both trees); fields x, y, z, RGB, root_src, contested, hag. Without: x, y, z, hag."""
    rooted = roots is not None
    g = cm.ground(plot)
    if rooted:
        z, rxyz, _ = roots
        cont_all = output.contested_rows(z["fi"])
        o = np.argsort(z["seg_id"], kind="stable")                   # per segment, in file order (fi sorted)
        u, first, cnt = np.unique(z["seg_id"][o], return_index=True, return_counts=True)
        span = dict(zip(u.tolist(), zip(first.tolist(), (first + cnt).tolist())))
    bk = sto.Buckets(tmp, meta, p["bucket_m"])
    bb = pd.read_csv(tmp / "seg_bbox.csv").set_index("seg_id")
    rows = []
    for q in cand.sort_values(["x", "y"]).itertuples():
        t0 = time.time()
        s = int(q.seg_id)
        b = bb.loc[s]
        own = bk.box([b.xmin - 0.01, b.ymin - 0.01], [b.xmax + 0.01, b.ymax + 0.01])
        own = own[own["id"] == s]
        if rooted:
            a_, b_ = span.get(s, (0, 0))
            k = o[a_:b_]
            rx, src_r, cont_r = rxyz[k], z["src"][k], cont_all[k]
        else:
            rx, src_r, cont_r = np.zeros((0, 3), np.float32), np.zeros(0, np.uint8), np.zeros(0, bool)
        xyz = np.r_[np.c_[own["x"], own["y"], own["z"]], rx].astype(np.float32)
        src = np.r_[np.zeros(len(own), np.uint8), src_r]
        cont = np.r_[np.zeros(len(own), bool), cont_r]
        code = np.where(cont, page.CONTESTED, src)
        tree_p, prev_p = paths[q.tree]
        tree_p.parent.mkdir(parents=True, exist_ok=True); prev_p.parent.mkdir(parents=True, exist_ok=True)
        fields = [("x", "<f4"), ("y", "<f4"), ("z", "<f4")]
        fields += [("red", "u1"), ("green", "u1"), ("blue", "u1"), ("root_src", "u1"), ("contested", "u1")] \
            if rooted else []
        rec = np.empty(len(xyz), fields + [("hag", "<f4")])
        rec["x"], rec["y"], rec["z"] = xyz.T
        if rooted:
            colors = np.array([page.COL[i] for i in range(5)], np.uint8)
            rec["red"], rec["green"], rec["blue"] = colors[code].T
            rec["root_src"], rec["contested"] = src, cont
        rec["hag"] = xyz[:, 2] - g(xyz[:, 0].astype(np.float64), xyz[:, 1].astype(np.float64))
        write_ply(tree_p.with_suffix(".part"), rec)                  # renamed when complete
        tree_p.with_suffix(".part").replace(tree_p)
        page.preview(prev_p, xyz, code, p["preview_points"])
        rows.append({"tree": q.tree, "n_rct": int((src == 0).sum()), "n_v1_only": int((src == 1).sum()),
                     "n_v2_only": int((src == 2).sum()), "n_both": int((src == 3).sum()),
                     "n_contested": int(cont.sum()), "seconds": round(time.time() - t0, 1)})
        print(f"{q.tree}: {len(xyz):,} pts, {time.time() - t0:.1f} s", flush=True)
    return pd.DataFrame(rows)


def run(plot, project=".", seg=None, boundary=None, no_roots=False, hints=None, plot_genus=None, existing=None,
        settings_path=None, redo=False):
    t0 = time.time()
    known = settings.inputs(project, plot, segmented=seg, boundary=boundary)
    seg = settings.need(known, "segmented", "review --segmented")
    boundary = Path(known["boundary"]) if "boundary" in known else None
    roots_npz = None if no_roots else settings.need(known, "roots", "roots") if "roots" in known else None
    if not no_roots and roots_npz is None:
        print(f"{plot}: no roots yet, review without roots")
    work = settings.work_dir(project, plot)
    trees = work / "trees_in_plot.csv"
    if not trees.exists():
        raise SystemExit(f"not found: {trees}; run select first")
    out = settings.review_dir(project, plot)
    if (out / "00_gallery.html").exists() and not redo:
        print(f"{plot}: review exists, skipped (--redo rewrites it and keeps checklist.csv)")
        return 0
    if redo:                                          # the tree files follow the roots; decisions and _clean files stay
        for p in [*(out / "trees").glob("*.ply"), *(out / "previews").glob("*.png")]:
            if not p.stem.endswith("_clean"):
                p.unlink()
    cfg = settings.load(settings_path)
    p = cfg["review"]
    root_cfg = config.load()
    config.apply(root_cfg)
    out.mkdir(parents=True, exist_ok=True)
    tmp = work / "tmp_review"
    roots_run.use_tmp(tmp)
    cm.setup(seg, boundary)
    rooted = roots_npz is not None
    roots = read_roots(plot, roots_npz) if rooted else None
    st = pd.read_csv(trees)
    meta = sto.store(plot, tmp, p["bucket_m"], p["margin_m"])
    nb = sto.boundary_counts(tmp, meta, p["cell_b_m"])
    t = st.copy()
    t["boundary_points"] = [int(nb.get(int(s), 0)) for s in t.seg_id]
    trunc = t.boundary_points >= np.maximum(p["trunc_min"], p["trunc_share"] * t.n_points)
    vol = root_volume(*roots, root_cfg["h_cut"], p["root_vox_m"]) if rooted else None
    t["root_vol_m3"] = [vol.get(int(s), 0.0) for s in t.seg_id] if rooted else np.nan
    t = genus_table(t, hints, plot_genus)
    ex = read_existing(existing)
    edge = t[trunc]
    cand = t[~trunc]
    low_h = cand[cand.H_m < p["h_min_m"]].copy()
    cand = cand[cand.H_m >= p["h_min_m"]].copy()
    old_ck = pd.read_csv(out / "checklist.csv", dtype=str).fillna("").set_index("tree") \
        if (out / "checklist.csv").exists() else pd.DataFrame(columns=["genus", "status", "notes"])
    old_tr = pd.read_csv(work / "trees.csv").set_index("tree") if (work / "trees.csv").exists() else pd.DataFrame()
    cand["flags"], impl = implausible(cand, ex, cfg["d_max_cm"])
    low_h["flags"] = implausible(low_h, ex, cfg["d_max_cm"])[0]
    if p["suggest_genus"] and rooted:
        cand = suggest_genus(cand, impl, p)
    for tr in cand.tree:                              # a genus edited in the checklist is kept
        if tr in old_ck.index and tr in old_tr.index and old_ck.at[tr, "genus"] not in ("", old_tr.at[tr, "genus"]):
            cand.loc[cand.tree == tr, ["genus", "genus_source"]] = [old_ck.at[tr, "genus"], "checklist"]
    cand["rooted"] = rooted
    cand["low"] = impl | ~np.isfinite(cand.D_cm)
    cand = priority(cand, ex, p["d_bins_cm"], p["h_bins_m"])
    cand["file"] = [f"trees/{q:03d}_{t_}.ply" for q, t_ in zip(cand.priority, cand.tree)]
    cand["preview"] = [f"previews/{t_}.png" for t_ in cand.tree]
    # existing tree files: renamed when the priority changed; missing or interrupted ones (no preview) written
    have = {f.name.split("_", 1)[1][:-4]: f for f in (out / "trees").glob("*.ply")
            if not f.stem.endswith("_clean")} if (out / "trees").exists() else {}
    todo = []
    for r in cand.itertuples():
        old = have.get(r.tree)
        if old is not None and (out / r.preview).exists() and r.tree in old_tr.index \
                and bool(old_tr.at[r.tree, "rooted"]) == r.rooted:
            if old != out / r.file:
                old.rename(out / r.file)
        else:
            todo.append(r.tree)
    print(f"{plot}: {'with' if rooted else 'without'} roots, {len(cand)} trees (H >= {p['h_min_m']:g} m; {len(low_h)} "
          f"lower; {len(edge)} edge-truncated), {len(todo)} to write", flush=True)
    paths = {t_: (out / f, out / q) for t_, f, q in zip(cand.tree, cand.file, cand.preview)}
    w = write_trees(plot, cand[cand.tree.isin(todo)], paths, roots, tmp, meta, p) if todo \
        else pd.DataFrame(columns=["tree"])
    cnt_cols = ["n_rct", "n_v1_only", "n_v2_only", "n_both", "n_contested"]
    prev = old_tr[[c for c in cnt_cols if c in old_tr.columns]].reset_index() if len(old_tr) else pd.DataFrame()
    counts = pd.concat([w[["tree"] + cnt_cols] if len(w) else w,
                        prev[~prev.tree.isin(w.tree)] if len(prev) else prev])
    counts = counts.reindex(columns=["tree"] + cnt_cols)   # empty when no tree is left to review
    tab = cand.merge(counts, on="tree", how="left").sort_values(["genus", "priority"])
    tab[cnt_cols] = tab[cnt_cols].apply(pd.to_numeric).round().astype("Int64")
    tab["segment_ids"] = tab.seg_id.astype(str)
    tab.to_csv(work / "trees.part", index=False, float_format="%.3f")
    (work / "trees.part").replace(work / "trees.csv")
    low_h[["tree", "seg_id", "H_m", "D_cm", "genus", "genus_source", "flags"]].sort_values("H_m", ascending=False) \
        .to_csv(work / "trees_below_3m.csv", index=False, float_format="%.2f")
    edge[["tree", "seg_id", "H_m", "D_cm", "n_points", "boundary_points"]].to_csv(
        work / "trees_edge_truncated.csv", index=False, float_format="%.2f")
    ck = tab[["priority", "tree", "genus", "genus_source", "segment_ids", "H_m", "D_cm", "n_contested", "flags",
              "file"]].rename(columns={"H_m": "H", "D_cm": "D"}).round({"H": 1, "D": 1})
    ck["status"] = [old_ck.at[t_, "status"] if t_ in old_ck.index else "" for t_ in ck.tree]
    ck["notes"] = [old_ck.at[t_, "notes"] if t_ in old_ck.index else "" for t_ in ck.tree]
    gone = old_ck[~old_ck.index.isin(ck.tree)]
    if len(gone):                                     # never lose a decision
        g_ = gone.reset_index()
        g_["notes"] = (g_["notes"] + " [not proposed any more]").str.replace(
            r"( \[not proposed any more\])+$", " [not proposed any more]", regex=True).str.strip()
        ck = pd.concat([ck, g_[[c for c in ck.columns if c in g_.columns]]], ignore_index=True)
    ck.to_csv(out / "checklist.part", index=False)
    (out / "checklist.part").replace(out / "checklist.csv")
    (work / "run.json").write_text(json.dumps(known | {"n_selected": len(st)}, indent=1))
    header = (f'<div class="note">{"Own segment plus root proposals." if rooted else "Own segment only, no roots."} '
              f"Tree files: trees/&lt;priority&gt;_&lt;tree&gt;.ply. Not shown: {len(low_h)} trees below "
              f"{p['h_min_m']:g} m and {len(edge)} edge-truncated trees (lists in ../../work/{plot}/).</div>")
    page.gallery(tab, time.time() - t0, out / "00_gallery.html", f"Tree review {plot}", header, rooted)
    shutil.rmtree(tmp)
    print(tab.groupby(["genus", "genus_source"]).size().to_string())
    print(f"{plot}: {len(tab)} trees, {time.time() - t0:.0f} s", flush=True)
    return 0
