"""rayroots-single split | model: from finished trees to volumes. Not validated yet, see README.
Trees from a project's tree_register.csv (--project), or straight from files (--input; no register needed).
Exit codes: 0 done; 1 a tree failed."""
import argparse
import sys
from pathlib import Path


def main():
    a = argparse.ArgumentParser(prog="rayroots-single", description=__doc__)
    sub = a.add_subparsers(dest="step", required=True)
    s = sub.add_parser("split", help="find the cut; stems and roots files, one QC figure per tree")
    m = sub.add_parser("model", help="RCT per part; tree_volumes.csv")
    m.add_argument("--ensemble", action="store_true", help="15 runs instead of the central run only")
    m.add_argument("--keep", action="store_true", help="keep the ensemble's tree files")
    s.add_argument("--genus", help="genus of all --input files (default: the folder name, e.g. Rhizophora)")
    for p in (s, m):
        p.add_argument("--project", type=Path, help="project with A_plot/tree_register.csv (default: current folder)")
        p.add_argument("--input", type=Path, help="one tree .ply, or a folder of them; tree = file name without _clean")
        p.add_argument("--out", type=Path, help="with --input: output folder (default: <input folder>/rayroots_single)")
        p.add_argument("--tree", help="this tree only; it is redone")
        p.add_argument("--redo", action="store_true", help="all trees again")
        p.add_argument("--settings", type=Path)
    x = a.parse_args()
    if x.input and x.project:
        a.error("give --project or --input, not both")
    if not x.input and (x.out or getattr(x, "genus", None)):
        a.error("--out and --genus go with --input")
    project = x.project or Path(".")
    if x.step == "split":
        from rayroots.single import split
        sys.exit(split.run(project, x.tree, x.settings, x.redo, x.input, x.genus, x.out))
    from rayroots.single import model
    sys.exit(model.run(project, x.tree, x.ensemble, x.settings, x.redo, x.keep, x.input, x.out))


if __name__ == "__main__":
    main()
