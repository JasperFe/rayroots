"""finish (A5): reads a plot's checklist.csv.
  ok    copy the tree .ply (xyz only) to finished_trees/<genus>/<tree>.ply
  fix   copy trees/<tree>_clean.ply there once it exists, else pending
  drop  log it
An existing finished_trees file is never overwritten. Every action goes to work/<plot>/review_log.csv; trees
already copied are skipped on a rerun. Copied and dropped trees go to tree_register.csv, one row per tree. finished_trees
and the register are the input of the single-tree steps."""
import datetime
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from rayroots.plot import settings

LOG_COLS = ["date", "tree", "status", "action", "target", "note"]
REG_COLS = ["plot", "tree", "segment_ids", "genus", "status", "date"]


def read_ply(path):
    """Vertex records of a binary little-endian PLY (as written by review / CloudCompare)."""
    types = {"float": "<f4", "float32": "<f4", "double": "<f8", "float64": "<f8", "uchar": "u1", "uint8": "u1",
             "char": "i1", "int8": "i1", "ushort": "<u2", "uint16": "<u2", "short": "<i2", "int16": "<i2",
             "int": "<i4", "int32": "<i4", "uint": "<u4", "uint32": "<u4"}
    with open(path, "rb") as f:
        props, n, inv = [], 0, False
        while True:
            line = f.readline().decode("ascii", "replace").strip()
            if line.startswith("format") and "binary_little_endian" not in line:
                raise ValueError(f"{path}: only binary little-endian PLY")
            if line.startswith("element vertex"):
                n = int(line.split()[2]); inv = True
            elif line.startswith("element"):
                inv = False
            elif line.startswith("property") and inv:
                _, t, name = line.split()
                props.append((name, types[t]))
            if line == "end_header":
                break
        return np.fromfile(f, np.dtype(props), count=n)


def write_xyz_ply(path, rec):
    out = np.empty(len(rec), [("x", "<f4"), ("y", "<f4"), ("z", "<f4")])
    out["x"], out["y"], out["z"] = rec["x"], rec["y"], rec["z"]
    part = path.with_suffix(".part")
    with open(part, "wb") as f:
        f.write(f"ply\nformat binary_little_endian 1.0\nelement vertex {len(out)}\nproperty float x\n"
                f"property float y\nproperty float z\nend_header\n".encode())
        out.tofile(f)
    part.replace(path)


def plan(plot, project):
    """One row per checklist decision not handled yet: action ('copy', 'copy clean', 'dropped', 'skipped', 'pending',
    'unknown status'), source and target."""
    rv, work = settings.review_dir(project, plot), settings.work_dir(project, plot)
    fl = pd.read_csv(work / "trees.csv").set_index("tree")["file"]
    ck = pd.read_csv(rv / "checklist.csv", dtype=str).fillna("")
    log_p = work / "review_log.csv"
    log = pd.read_csv(log_p, dtype=str).fillna("") if log_p.exists() else pd.DataFrame(columns=LOG_COLS)
    done = set(log.loc[log.action == "copied", "tree"])
    print(ck.status.str.strip().str.lower().replace("", "(empty)").value_counts().to_string())
    rows = []
    for r in ck.itertuples():
        status = r.status.strip().lower()
        if not status or r.tree in done:
            continue
        target_p = settings.part(project) / "finished_trees" / r.genus / f"{r.tree}.ply"
        clean_p = rv / "trees" / f"{r.tree}_clean.ply"
        row = {"tree": r.tree, "status": status, "genus": r.genus, "source": "", "target": "", "note": r.notes,
               "plot": plot, "segment_ids": r.segment_ids, "n_contested": r.n_contested}
        if status == "drop":
            row["action"] = "dropped"
        elif status not in ("ok", "fix"):
            row["action"], row["note"] = "unknown status", f"status '{r.status}'"
        elif r.genus in ("", "unknown"):
            row["action"], row["note"] = "pending", "set the genus in checklist.csv"
        elif target_p.exists():
            row["action"], row["target"], row["note"] = "skipped", str(target_p), "finished_trees file exists (not overwritten)"
        elif status == "ok" and r.tree not in fl.index:
            row["action"], row["note"] = "pending", "tree no longer proposed (not in trees.csv)"
        elif status == "ok":
            row["action"], row["source"], row["target"] = "copy", str(rv / fl[r.tree]), str(target_p)
        elif clean_p.exists():
            row["action"], row["source"], row["target"] = "copy clean", str(clean_p), str(target_p)
        else:
            row["action"], row["note"] = "pending", f"{r.tree}_clean.ply not there yet"
        rows.append(row)
    return work, log, pd.DataFrame(rows, columns=["tree", "status", "genus", "action", "source", "target", "note",
                                                "plot", "segment_ids", "n_contested"])


def apply(work, log, todo, project):
    """Copy and log; pending / dropped / skipped / unknown rows are logged again only when their state changes."""
    today = datetime.date.today().isoformat()
    new = []
    for r in todo.itertuples():
        action = r.action
        if action in ("copy", "copy clean"):
            Path(r.target).parent.mkdir(parents=True, exist_ok=True)
        if action == "copy":
            write_xyz_ply(Path(r.target), read_ply(r.source)); action = "copied"
        elif action == "copy clean":
            shutil.copy2(r.source, r.target); action = "copied"
        new.append({"date": today, "tree": r.tree, "status": r.status, "action": action, "target": r.target,
                    "note": r.note})
    new = pd.DataFrame(new, columns=LOG_COLS)
    last = log.drop_duplicates("tree", keep="last").set_index("tree") if len(log) else None
    if last is not None and len(new):
        keep = [not (t in last.index and last.at[t, "action"] == a and a in ("pending", "dropped", "skipped",
                                                                             "unknown status"))
                for t, a in zip(new.tree, new.action)]
        new = new[keep]
    pd.concat([log, new], ignore_index=True).to_csv(work / "review_log.csv", index=False)
    register(settings.part(project) / "tree_register.csv",
             todo[todo.tree.isin(new.loc[new.action.isin(["copied", "dropped"]), "tree"])], today)
    print(new.groupby("action").size().to_string() if len(new) else "nothing new")


def register(path, done, today):
    """Upsert the copied / dropped trees (one row per tree; a later decision replaces it)."""
    if not len(done):
        return
    reg = pd.read_csv(path, dtype=str).fillna("") if path.exists() else pd.DataFrame(columns=REG_COLS)
    rows = pd.DataFrame({"plot": done["plot"], "tree": done.tree, "segment_ids": done.segment_ids, "genus": done.genus,
                         "status": done.status, "date": today}, columns=REG_COLS)
    reg = pd.concat([reg[~reg.tree.isin(rows.tree)], rows], ignore_index=True).sort_values(["plot", "tree"])
    reg.to_csv(path.with_suffix(".part"), index=False)
    path.with_suffix(".part").replace(path)
    print(f"tree_register.csv: {len(rows)} trees registered ({len(reg)} in total)")


def run(plot, project=".", do_apply=False):
    if not (settings.review_dir(project, plot) / "checklist.csv").exists():
        print(f"no checklist of {plot} in {settings.review_dir(project, plot)}")
        return 2
    work, log, todo = plan(plot, project)
    print(todo.groupby("action").size().to_string() if len(todo) else "nothing to do")
    pend = todo[todo.action == "pending"]
    print("pending:", "; ".join(f"{a} ({b})" for a, b in zip(pend.tree, pend.note)) if len(pend) else "none")
    c = todo[todo.action.isin(["copy", "copy clean"]) & (pd.to_numeric(todo.n_contested, errors="coerce") > 0)]
    if len(c):
        print(f"WARNING: {len(c)} trees hold contested points that may also be in a neighbouring tree: "
              + "; ".join(f"{a} ({b})" for a, b in zip(c.tree, c.n_contested)))
    if do_apply:
        apply(work, log, todo, project)
    else:
        print("dry run: nothing written; --apply copies")
    return 0
