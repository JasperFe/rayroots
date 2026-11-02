"""End-to-end check of the crop: tree count, points per tree file, checksum of the roots.
Usage: crop_check.py <A_plot/work/<plot>> <stem> <expected.json> [--write]. Needs a run with --keep (raysplit trees/)."""
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np

run, stem, exp = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
z = np.load(run / f"{stem}_roots.npz")
sha = hashlib.sha256(b"".join(z[k].tobytes() for k in ("fi", "seg_id", "src"))).hexdigest()
trees = {}
for p in sorted((run / "trees").glob("*.ply")):
    with open(p, "rb") as f:
        head = b""
        while not head.endswith(b"end_header\n"):
            head += f.readline()
    trees[p.stem.rsplit("_", 1)[1]] = int(re.search(rb"element vertex (\d+)", head).group(1))
got = {"tree_count": len(trees), "points_per_tree": trees, "root_rows": len(z["fi"]), "roots_sha256": sha}
if "--write" in sys.argv:
    exp.write_text(json.dumps(got, indent=1) + "\n")
    print("written", exp)
    sys.exit(0)
want = json.loads(exp.read_text())
bad = [k for k in want if want[k] != got[k]]
print(f"trees {got['tree_count']} (expected {want['tree_count']}), root rows {got['root_rows']:,}, sha {sha[:16]}")
print("CROP CHECK", "PASS" if not bad else f"FAIL: {bad}")
sys.exit(1 if bad else 0)
