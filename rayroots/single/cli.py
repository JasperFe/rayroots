"""rayroots-single split | model: from finished trees to volumes. Not validated yet, see README.
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
    for p in (s, m):
        p.add_argument("--project", type=Path, default=Path("."))
        p.add_argument("--tree", help="this tree only; it is redone")
        p.add_argument("--redo", action="store_true", help="all trees again")
        p.add_argument("--settings", type=Path)
    x = a.parse_args()
    if x.step == "split":
        from rayroots.single import split
        sys.exit(split.run(x.project, x.tree, x.settings, x.redo))
    from rayroots.single import model
    sys.exit(model.run(x.project, x.tree, x.ensemble, x.settings, x.redo, x.keep))


if __name__ == "__main__":
    main()
