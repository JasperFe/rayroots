# Frozen method (v1 + v2). Do not edit here: a change is a new version.
"""Run rayimport / rayextract trees."""
import shutil
import subprocess
from pathlib import Path

GROUND_OFFSET = 0.001   # m; the flat ground goes this far above the cloud's lowest point


def rct(cmd, cwd=None, timeout=3600):
    """Run an RCT tool; returns '' on success, else the last line of its output."""
    try:
        r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=timeout, cwd=cwd)
    except subprocess.TimeoutExpired:
        return f"timeout {timeout} s"
    if r.returncode != 0:
        out = (r.stderr.strip() or r.stdout.strip()).splitlines()
        return f"exit {r.returncode}: {out[-1] if out else ''}"
    return ""


def rayimport(ply, raycloud, ray_vector="0,0,-0.1"):
    """rayimport writes <name>_raycloud.ply next to its input; move it to `raycloud`. Returns error or ''."""
    made = Path(ply).with_name(Path(ply).stem + "_raycloud.ply")
    made.unlink(missing_ok=True)
    err = rct(["rayimport", ply, "ray", ray_vector, "--max_intensity", "0"])
    if err or not made.exists():
        return err or "rayimport: no output"
    shutil.move(made, raycloud)
    return ""


def rayextract_trees(raycloud, mesh, settings, out_prefix, keep=("trees.txt",)):
    """rayextract trees on `raycloud` with ground `mesh`; outputs (written next to the ray cloud) are
    renamed to <out_prefix>_<suffix> for the suffixes in `keep` and the rest deleted.
    Returns error or ''."""
    raycloud = Path(raycloud)
    stem = raycloud.stem
    made = {s: raycloud.with_name(f"{stem}_{s}") for s in ("trees.txt", "trees_mesh.ply", "segmented.ply")}
    for p in made.values():
        p.unlink(missing_ok=True)
    args = [a for k, v in settings.items() for a in (f"--{k}", str(v))]
    err = rct(["rayextract", "trees", raycloud, mesh, *args])
    ok = made["trees.txt"].exists()
    for s, p in made.items():
        if p.exists():
            if s in keep and ok:
                p.rename(f"{out_prefix}_{s}")
            else:
                p.unlink()
    return "" if ok else (err or "rayextract trees: no trees file")
