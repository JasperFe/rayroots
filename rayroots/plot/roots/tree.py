"""rayroots-plot tree: one file per segment: <stem>_tree_<seg_id>.ply with its RCT points and ALL its roots.npz rows. A contested point
(v1 gave it one tree, v2 another) is in both trees. All fields of the plot file, plus root_src and contested;
RGB shows the source (COL below)."""
import numpy as np
from numpy.lib.recfunctions import repack_fields

from rayroots.plot.roots import ground as rw
from rayroots.plot.roots import output
from rayroots.plot.roots import plotdata as cm

COL = {0: (127, 127, 127), 3: (0, 150, 60), 1: (235, 104, 52), 2: (42, 120, 214)}   # RCT, both, v1 only, v2 only
COL_CONTESTED = (220, 0, 200)


def write_tree(path, head, parts):
    a = np.concatenate(parts)
    a = a[np.argsort(a["fi"], kind="stable")]
    names = [n for n in a.dtype.names if n != "fi"]
    for code, col in COL.items():
        k = a["root_src"] == code
        a["red"][k], a["green"][k], a["blue"][k] = col
    k = a["contested"] == 1
    a["red"][k], a["green"][k], a["blue"][k] = COL_CONTESTED
    lines = head.decode("ascii").splitlines(keepends=True)
    lines = [f"element vertex {len(a)}\n" if s.startswith("element vertex") else s for s in lines]
    i = lines.index("end_header\n")
    lines.insert(i, "property uchar contested\n")
    with open(path, "wb") as f:
        f.write("".join(lines).encode("ascii"))
        f.write(repack_fields(a[names]).tobytes())


def run(roots_ply, seg_ids, npz, out):
    stem = roots_ply.stem.removesuffix("_segmented_roots")
    z = np.load(npz)
    fi, sid, src = z["fi"], z["seg_id"], z["src"]
    cont = output.contested_rows(fi)
    want = sorted(set(seg_ids))
    rows = {s: np.flatnonzero(sid == s) for s in want}
    mm, head = cm.memmap(roots_ply)
    if "root_src" not in mm.dtype.names:
        raise SystemExit(f"{roots_ply.name} has no root_src: give the rayroots output, not the rayextract file")
    dt = np.dtype(mm.dtype.descr + [("contested", "u1"), ("fi", "<i8")])
    parts = {s: [] for s in want}

    def add(s, recs, index, src_v, cont_v):
        o = np.zeros(len(recs), dt)
        for n in mm.dtype.names:
            o[n] = recs[n]
        o["root_src"], o["contested"], o["fi"] = src_v, cont_v, index
        parts[s].append(o)

    for s0 in range(0, len(mm), cm.CHUNK):
        c = np.asarray(mm[s0:s0 + cm.CHUNK])
        ids = rw.decode_ids(c["red"], c["green"], c["blue"])
        for s in want:
            k = np.flatnonzero((ids == s) & (c["root_src"] == 0))       # its RCT points
            if len(k):
                add(s, c[k], k + s0, 0, 0)
            r = rows[s]
            r = r[(fi[r] >= s0) & (fi[r] < s0 + len(c))]                # its roots, contested ones included
            if len(r):
                add(s, c[fi[r].astype(np.int64) - s0], fi[r], src[r], cont[r].astype(np.uint8))
    out.mkdir(parents=True, exist_ok=True)
    for s in want:
        if not parts[s]:
            print(f"segment {s}: no points")
            continue
        path = out / f"{stem}_tree_{s}.ply"
        write_tree(path, head, parts[s])
        n = np.concatenate(parts[s])
        print(f"{path.name}: {len(n):,} points: RCT {(n['root_src'] == 0).sum():,}, v1 only {(n['root_src'] == 1).sum():,}, "
              f"v2 only {(n['root_src'] == 2).sum():,}, both {(n['root_src'] == 3).sum():,}, "
              f"contested {(n['contested'] == 1).sum():,}", flush=True)
