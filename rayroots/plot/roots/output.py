"""Outputs of a run: <stem>_roots.npz, <stem>_segmented_roots.ply, <stem>_roots.csv."""
import numpy as np
import pandas as pd

from rayroots.plot.roots import cuts
from rayroots.plot.roots import plotdata as cm


def roots_table(ids, fi, lab1, lab2):
    """(fi, seg_id, src) of every point a method attaches to a segment other than its RCT one, sorted by fi."""
    parts = []
    for lab, code in ((lab1, 1), (lab2, 2)):
        k = np.flatnonzero((lab >= 0) & (lab != ids))
        parts.append(np.c_[fi[k].astype(np.int64), lab[k].astype(np.int64), np.full(len(k), code)])
    a = np.vstack(parts)
    key, inv = np.unique(a[:, :2], axis=0, return_inverse=True)       # sorted by fi, then seg_id
    src = np.bincount(inv.ravel(), a[:, 2], len(key)).astype(np.uint8)  # 1 + 2 = 3: both
    return key[:, 0].astype(np.uint32), key[:, 1].astype(np.int32), src


def contested_rows(fi):
    """Bool per row: the row's point has rows for two segments (v1 gave it one tree, v2 another)."""
    same = fi[1:] == fi[:-1]
    return np.r_[False, same] | np.r_[same, False]


def write_npz(path, fi, sid, src, seg, config_hash, frozen):
    part = path.with_name(path.stem + ".part.npz")
    np.savez(part, fi=fi, seg_id=sid, src=src, segmented_file=seg.name, segmented_bytes=seg.stat().st_size,
             config_hash=config_hash, frozen=frozen)
    part.replace(path)


def write_ply(path, seg, fi, sid, src):
    """The segmented cloud, same point order, root points recoloured, plus field root_src. A contested point
    (rows for two segments) goes to v2's segment. Returns points per segment after recolouring {seg_id: n}."""
    cont = contested_rows(fi)
    k = ~(cont & (src == 1))             # v2 wins
    pf, ps, pc = fi[k], sid[k], src[k]
    r, g, b = cm.encode_ids(ps)
    mm, head = cm.memmap(seg)
    text = head.decode("ascii")
    if text.count("element ") != 1 or "binary_little_endian" not in text:
        raise SystemExit(f"{seg.name}: expected a binary little-endian ply with only a vertex element")
    if "root_src" in text:
        raise SystemExit(f"{seg.name}: already has root_src; give the rayextract trees output")
    dt = np.dtype(mm.dtype.descr + [("root_src", "u1")])
    n_after = np.zeros(1 << 24, np.int64)
    part = path.with_name(path.name + ".part")
    with open(part, "wb") as f:
        f.write(text.replace("end_header\n", "property uchar root_src\nend_header\n").encode("ascii"))
        for s in range(0, len(mm), cm.CHUNK):
            c = np.asarray(mm[s:s + cm.CHUNK])
            out = np.zeros(len(c), dt)
            for name in mm.dtype.names:
                out[name] = c[name]
            i0, i1 = np.searchsorted(pf, [s, s + len(c)])
            w = pf[i0:i1].astype(np.int64) - s
            out["red"][w], out["green"][w], out["blue"][w] = r[i0:i1], g[i0:i1], b[i0:i1]
            out["root_src"][w] = pc[i0:i1]
            ids = cm.rw.decode_ids(out["red"], out["green"], out["blue"])
            n_after += np.bincount(ids[ids >= 0], minlength=len(n_after))
            f.write(out.tobytes())
    part.replace(path)
    u = np.flatnonzero(n_after)
    return dict(zip(u.tolist(), n_after[u].tolist()))


def write_csv(path, plot, fi, sid, src, n_after, failed_xy):
    """One row per segment. failed_xy: list of (lo, hi) boxes of failed flip tiles."""
    st = cuts.segments(plot)
    tracks = cuts.stems(plot)
    g = cm.ground(plot)
    st["H_m"] = st.z_max - g(st.base_x.values, st.base_y.values)
    cont = contested_rows(fi)
    rows = pd.DataFrame({"seg_id": sid, "src": src, "cont": cont})
    tab = rows.groupby(["seg_id", "src"]).size().unstack(fill_value=0).reindex(columns=[1, 2, 3], fill_value=0)
    nc = rows[rows.cont].groupby("seg_id").size()
    out = pd.DataFrame({"seg_id": st.seg_id, "n_rct": st.n_points, "H_m": st.H_m, "base_x": st.base_x,
                        "base_y": st.base_y, "dist_to_boundary": st.dist_to_hull})
    for c, name in ((1, "n_v1_only"), (2, "n_v2_only"), (3, "n_both")):
        out[name] = out.seg_id.map(tab[c]).fillna(0).astype(np.int64)
    out["n_contested"] = out.seg_id.map(nc).fillna(0).astype(np.int64)
    flags = []
    for q in out.itertuples():
        f = []
        if q.dist_to_boundary > cuts.BUFFER:
            f.append("outside_boundary")
        if q.seg_id not in tracks:
            f.append("no_stem_track")
        if n_after.get(q.seg_id, 0) == 0:
            f.append("absorbed")
        if q.n_contested:
            f.append("contested")
        if any((lo <= [q.base_x, q.base_y]).all() and ([q.base_x, q.base_y] < hi).all() for lo, hi in failed_xy):
            f.append("tile_failed")
        flags.append(";".join(f))
    out["flags"] = flags
    cols = ["seg_id", "n_rct", "n_v1_only", "n_v2_only", "n_both", "n_contested", "H_m", "base_x", "base_y",
            "dist_to_boundary", "flags"]
    part = path.with_name(path.name + ".part")
    out[cols].to_csv(part, index=False, float_format="%.3f")
    part.replace(path)
