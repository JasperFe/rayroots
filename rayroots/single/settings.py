"""Settings, the finished trees of a project, and overrides.csv."""
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import yaml


def load(path=None):
    cfg = yaml.safe_load(open(Path(__file__).with_name("single.yaml")))
    if path:
        new = yaml.safe_load(open(path)) or {}
        for k, v in new.pop("rct", {}).items():
            cfg["rct"][k] = cfg["rct"][k] | v if isinstance(v, dict) else v
        cfg |= new
    return cfg


def part(project):
    return Path(project) / "B_trees"


def work(project):
    return part(project) / "work"


def trees(project, cfg, name=None):
    """Finished trees (status ok or fix) of tree_register.csv with a known genus; other genera skipped with a message."""
    p = Path(project) / "A_plot" / "tree_register.csv"
    if not p.exists():
        raise SystemExit(f"not found: {p} (rayroots-plot finish writes it)")
    reg = pd.read_csv(p, dtype=str).fillna("")
    reg = reg[reg.status.isin(["ok", "fix"])]
    if name:
        reg = reg[reg.tree == name]
        if reg.empty:
            raise SystemExit(f"{name}: not a finished tree in {p}")
    out = []
    for r in reg.itertuples():
        if r.genus not in cfg["genus"]:
            print(f"{r.tree}: genus {r.genus} skipped (no rule in single.yaml)")
            continue
        path = Path(project) / "A_plot" / "finished_trees" / r.genus / f"{r.tree}.ply"
        if not path.exists():
            print(f"{r.tree}: {path} not found, skipped")
            continue
        out.append(SimpleNamespace(tree=r.tree, genus=r.genus, path=path))
    return out


def overrides(project):
    """overrides.csv in the project: tree, h_cut_above_z0 (roots), n_stems (dbh), reason. Returns {tree: {key: value}}."""
    p = part(project) / "overrides.csv"
    if not p.exists():
        return {}
    t = pd.read_csv(p)
    if "h_cut" in t:
        raise SystemExit(f"{p}: column h_cut is now h_cut_above_z0 (z0 = 0.1th percentile of z, in split_log.csv)")
    return {r["tree"]: {k: r[k] for k in ("h_cut_above_z0", "n_stems") if k in r and pd.notna(r[k])} for _, r in t.iterrows()}
