"""segment: plot LAS -> rayimport -> rayextract terrain -> rayextract trees, in A_plot/work/<plot>/. Settings: `rct` in
plot.yaml, passed through as flags. A step whose output exists is skipped. The LAS copy, the ray cloud and RCT's tree mesh
are deleted at the end, unless --keep. Own rayextract output: skip this step, give --segmented to `roots`."""
import shutil
import struct
import subprocess
from datetime import datetime
from pathlib import Path

from rayroots.plot import settings


def flags(d):
    """{key: value} -> --key value; True -> --key; None or False -> left out."""
    out = []
    for k, v in d.items():
        if v is None or v is False:
            continue
        out += [f"--{k}"] + ([] if v is True else [str(v)])
    return out


def header(path):
    """LAS 1.x public header: point count, scale, the exact offset."""
    b = open(path, "rb").read(227)
    return (f"las_points: {struct.unpack('<I', b[107:111])[0]}\n"
            "las_scale: %r %r %r\n" % struct.unpack("<3d", b[131:155]) +
            "las_offset_exact: %r %r %r\n" % struct.unpack("<3d", b[155:179]))


def call(cmd, w, log):
    """Run in w; the command line and the output go to the log (with time -v: memory peak, wall clock)."""
    print(f"== {datetime.now():%H:%M:%S} {' '.join(cmd)}", flush=True)
    timed = ["/usr/bin/time", "-v"] + cmd if Path("/usr/bin/time").exists() else cmd
    with open(log, "a") as f:
        f.write(f"command: {' '.join(cmd)}\n")
        f.flush()
        return subprocess.run(timed, cwd=w, stdout=f, stderr=subprocess.STDOUT).returncode


def run(las, plot=None, project=".", settings_path=None, redo=False, keep=False):
    las = Path(las).resolve()
    plot = plot or las.stem
    rct = settings.load(settings_path)["rct"]
    w = settings.work_dir(project, plot)
    w.mkdir(parents=True, exist_ok=True)
    seg, mesh, ray = w / f"{plot}_segmented.ply", w / f"{plot}_mesh.ply", w / f"{plot}.ply"
    settings.inputs(project, plot, las=las)
    if redo:
        for p in (seg, mesh, ray):
            p.unlink(missing_ok=True)
    if seg.exists():
        print(f"{plot}: segmented cloud exists, skipped (--redo runs it again)")
        settings.inputs(project, plot, segmented=seg, mesh=mesh if mesh.exists() else None)
        return 0
    for t in ("rayimport", "rayextract"):
        if not shutil.which(t):
            print(f"{t} not on PATH (RCT installed?)")
            return 2
    log = w / f"{plot}.rct.log"
    local = w / f"{plot}.las"
    if not ray.exists() and not local.exists():
        try:
            local.symlink_to(las)
        except OSError:                                # Windows drives refuse symlinks
            shutil.copy(las, local)
    with open(log, "a") as f:
        f.write(f"date: {datetime.now().isoformat()}\nrayextract: {Path(shutil.which('rayextract')).resolve()}\n")
        if local.exists():
            f.write(header(local))
    imp = dict(rct["rayimport"])
    vec = imp.pop("ray_vector")
    steps = [(ray, ["rayimport", local.name, "ray", vec] + flags(imp)),
             (mesh, ["rayextract", "terrain", ray.name] + flags(rct["terrain"])),
             (seg, ["rayextract", "trees", ray.name, mesh.name] + flags(rct["trees"]))]
    for out, cmd in steps:
        if not out.exists():
            rc = call(cmd, w, log)
            if rc or not out.exists():
                print(f"{plot}: {cmd[0]} failed (exit {rc}), see {log}")
                return 1
    settings.inputs(project, plot, segmented=seg, mesh=mesh)
    if not keep:
        for p in (local, ray, w / f"{plot}_trees_mesh.ply"):
            p.unlink(missing_ok=True)
    print(f"{plot}: segmented -> {seg}")
    return 0
