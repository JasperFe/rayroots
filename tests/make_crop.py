"""Crop a plot LAS at full density, same header offset and scale, so the crop stays in the frame of the full plot.
Box in the segmented frame (LAS minus header offset).
Usage: make_crop.py <in.las> <out.las> xmin xmax ymin ymax"""
import sys

import laspy
import numpy as np

src, dst = sys.argv[1:3]
x0, x1, y0, y1 = (float(v) for v in sys.argv[3:7])
with laspy.open(src) as f:
    h = f.header
    out = laspy.LasHeader(version=h.version, point_format=h.point_format)
    out.offsets, out.scales = h.offsets, h.scales
    n = 0
    with laspy.open(dst, mode="w", header=out) as w:
        for p in f.chunk_iterator(5_000_000):
            x, y = np.asarray(p.x) - h.offsets[0], np.asarray(p.y) - h.offsets[1]
            k = (x >= x0) & (x < x1) & (y >= y0) & (y < y1)
            if k.any():
                w.write_points(p[k])
                n += int(k.sum())
print(f"{n:,} points -> {dst}")
