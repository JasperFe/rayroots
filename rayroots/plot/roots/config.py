"""roots.yaml -> the constants of the frozen modules. The modules keep their literal values as written in the
source; tests/test_config.py checks that the defaults here equal them."""
import hashlib
import json
from pathlib import Path

import numpy as np
import yaml

from rayroots.plot.roots import connect as cn
from rayroots.plot.roots import cuts, flip, ground, v2

DEFAULTS = Path(__file__).with_name("roots.yaml")


def load(path=None):
    """Defaults, overlaid with the keys of `path`. An unknown key is an error."""
    cfg = yaml.safe_load(DEFAULTS.read_text())
    if path is None:
        return cfg
    for k, v in (yaml.safe_load(Path(path).read_text()) or {}).items():
        if k not in cfg:
            raise SystemExit(f"{path}: unknown key {k}")
        if isinstance(cfg[k], dict):
            bad = set(v) - set(cfg[k])
            if bad:
                raise SystemExit(f"{path}: unknown keys in {k}: {sorted(bad)}")
            cfg[k].update(v)
        else:
            cfg[k] = v
    return cfg


def method_part(cfg):
    return {k: v for k, v in cfg.items() if k != "ground"}


def digest(cfg):
    return hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:12]


def is_frozen(cfg):
    """True when the method settings are the frozen ones and the ground is the grid."""
    return cfg["ground"] == "grid" and method_part(cfg) == method_part(load())


def apply(cfg):
    """Set the module constants, derived ones included."""
    for k, v in cfg["grid"].items():
        setattr(ground, k, v)
    c = dict(cfg["cuts"])
    region_pad, assign_pad = c.pop("REGION_PAD"), c.pop("ASSIGN_PAD")
    for k, v in c.items():
        setattr(cuts, k, v)
    cuts.WINDOW = cuts.BUFFER + cuts.R_MAX
    cuts.REGION = cuts.WINDOW + region_pad
    cuts.ASSIGN_REGION = cuts.BUFFER + cuts.R_MAX + assign_pad
    cuts.HEIGHTS = np.round(np.arange(cuts.H_START, cuts.S0 + 1e-9, cuts.STEP), 4)
    for k, v in cfg["flip"].items():
        setattr(flip, k, v)
    flip.H_CUT, flip.CYL_TOL, flip.Z_MIN = cfg["h_cut"], cfg["cyl_tol"], cuts.Z_MIN
    for k, v in cfg["connect"].items():
        setattr(cn, k, v)
    cn.R_REACH = cuts.R_MAX
    w = cfg["v2"]
    v2.H_CUT, v2.CORE, v2.COL_PAD = cfg["h_cut"], w["CORE"], w["COL_PAD"]
    v2.PAD = cn.R_REACH + w["PAD_EXTRA"]
    v2.DEFAULT = (w["TOL"], w["TURN_MAX"], w["BRANCH_MAX"])
