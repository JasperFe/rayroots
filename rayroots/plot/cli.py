"""rayroots-plot: from a plot cloud to finished trees. Each step has its own inputs and outputs, in <project>/A_plot/.
A step whose output exists is skipped; --redo runs it again. A path given to one step is remembered for the next.
Exit codes: 0 done; 1 failed; 2 bad input; 3 memory guard (roots: rerun, it resumes); 5 outputs written, but tiles failed."""
import argparse
import sys
from pathlib import Path


def add(sub, name, help_, *args):
    p = sub.add_parser(name, help=help_)
    p.add_argument("--plot", required=name != "segment" and name != "overview", help="plot name")
    p.add_argument("--project", type=Path, default=Path("."))
    for a in args:
        if a == "settings":
            p.add_argument("--settings", type=Path, help="yaml with the settings to change (see plot.yaml)")
        elif a == "redo":
            p.add_argument("--redo", action="store_true", help="run this step again")
        elif a == "keep":
            p.add_argument("--keep", action="store_true", help="keep the large intermediate files")
        elif a == "segmented":
            p.add_argument("--segmented", type=Path, help="segmented .ply, e.g. your own rayextract output")
        elif a == "boundary":
            p.add_argument("--boundary", type=Path, help="csv x,y in the segmented frame; its convex hull is the plot")
    return p


def main():
    a = argparse.ArgumentParser(prog="rayroots-plot", description=__doc__)
    sub = a.add_subparsers(dest="step", required=True)
    p = add(sub, "segment", "LAS -> RCT segmented cloud (rct in plot.yaml)", "settings", "redo", "keep")
    p.add_argument("las", type=Path)
    p = add(sub, "roots", "root proposals: roots.npz, roots.csv", "segmented", "boundary", "redo", "keep")
    p.add_argument("--mesh", type=Path, help="terrain mesh, used with ground = mesh")
    p.add_argument("--config", type=Path, help="yaml with the method settings to change (see roots.yaml)")
    p.add_argument("--ground", choices=["grid", "mesh"])
    add(sub, "select", "trees of the plot", "segmented", "boundary", "settings", "redo")
    p = add(sub, "review", "gallery, checklist, tree files", "segmented", "boundary", "settings", "redo")
    p.add_argument("--no-roots", action="store_true", help="own segment only")
    p.add_argument("--genus-hints", type=Path, help="csv: tree or seg_id, genus")
    p.add_argument("--plot-genus", help="genus of every tree without a hint")
    p.add_argument("--existing-trees", type=Path, help="csv: genus, H_m, D_cm of trees that count as covered")
    p = add(sub, "context", "surroundings of trees marked fix, or of the trees named", "segmented", "settings")
    p.add_argument("trees", nargs="*", help="tree names or segment ids; default: all marked fix")
    p.add_argument("--radius", type=float)
    p = add(sub, "tree", "a tree with root_src and contested fields; needs roots --keep")
    p.add_argument("seg_ids", type=int, nargs="+")
    p = add(sub, "finish", "dry run; --apply copies to finished_trees/<genus>/ and tree_register.csv")
    p.add_argument("--apply", action="store_true")
    p = add(sub, "overview", "overview.html", "settings")
    p.add_argument("--plots", type=Path, help="plots.csv: also lists plots without a review")
    x = a.parse_args()
    for f in ("segmented", "boundary", "mesh", "config", "genus_hints", "existing_trees", "settings", "las"):
        v = getattr(x, f, None)
        if v is not None and not v.is_file():
            print(f"not found: {v}")
            sys.exit(2)
    if x.step == "segment":
        from rayroots.plot import segment
        sys.exit(segment.run(x.las, x.plot, x.project, x.settings, x.redo, x.keep))
    if x.step == "roots":
        from rayroots.plot.roots import run
        sys.exit(run.step(x.plot, x.project, x.segmented, x.mesh, x.boundary, x.config, x.ground, x.redo, x.keep))
    if x.step == "select":
        from rayroots.plot import select
        sys.exit(select.run(x.plot, x.project, x.segmented, x.boundary, x.settings, x.redo))
    if x.step == "review":
        from rayroots.plot import review
        sys.exit(review.run(x.plot, x.project, x.segmented, x.boundary, x.no_roots, x.genus_hints, x.plot_genus,
                            x.existing_trees, x.settings, x.redo))
    if x.step == "context":
        from rayroots.plot import context
        sys.exit(context.run(x.plot, x.trees or None, x.project, x.radius, x.segmented, x.settings))
    if x.step == "tree":
        from rayroots.plot import settings
        from rayroots.plot.roots import tree
        r = settings.inputs(x.project, x.plot)
        npz = settings.need(r, "roots", "roots")
        w = settings.work_dir(x.project, x.plot)
        ply = w / npz.name.replace("_roots.npz", "_segmented_roots.ply")
        if not ply.is_file():
            print(f"not found: {ply}; run roots with --redo --keep")
            sys.exit(2)
        sys.exit(tree.run(ply, x.seg_ids, npz, w / "tree"))
    if x.step == "finish":
        from rayroots.plot import finish
        sys.exit(finish.run(x.plot, x.project, x.apply))
    from rayroots.plot import overview
    sys.exit(overview.run(x.project, x.plots, x.settings))


if __name__ == "__main__":
    main()
