"""plot.yaml, overlaid with --settings. Project folders, and run.json: the inputs of a plot, given once."""
import json
from pathlib import Path

import yaml

DEFAULTS = Path(__file__).with_name("plot.yaml")


def load(path=None):
    cfg = yaml.safe_load(DEFAULTS.read_text())
    if path is None:
        return cfg
    for k, v in (yaml.safe_load(Path(path).read_text()) or {}).items():
        if k not in cfg:
            raise SystemExit(f"{path}: unknown key {k}")
        if k == "rct":                                 # passed through to RCT: any flag
            bad = set(v) - set(cfg[k])
            if bad:
                raise SystemExit(f"{path}: unknown keys in rct: {sorted(bad)}")
            for tool, flags in v.items():
                cfg[k][tool].update(flags or {})
        elif isinstance(cfg[k], dict):
            bad = set(v) - set(cfg[k])
            if bad:
                raise SystemExit(f"{path}: unknown keys in {k}: {sorted(bad)}")
            cfg[k].update(v)
        else:
            cfg[k] = v
    return cfg


def part(project):
    return Path(project) / "A_plot"


def review_dir(project, plot):
    return part(project) / "review" / plot


def work_dir(project, plot):
    return part(project) / "work" / plot


def inputs(project, plot, **given):
    """run.json of a plot, updated with the paths given (None: keep what is there). Returns its content."""
    w = work_dir(project, plot)
    w.mkdir(parents=True, exist_ok=True)
    p = w / "run.json"
    r = json.loads(p.read_text()) if p.exists() else {"plot": plot}
    r |= {k: str(Path(v).resolve()) for k, v in given.items() if v is not None}
    p.write_text(json.dumps(r, indent=1))
    return r


def need(r, key, step):
    """A path from run.json as a Path."""
    if key not in r:
        raise SystemExit(f"{r['plot']}: no {key} known; give it to {step} or run the earlier step")
    p = Path(r[key])
    if not p.is_file():
        raise SystemExit(f"not found: {p} (moved? give it again to {step})")
    return p
