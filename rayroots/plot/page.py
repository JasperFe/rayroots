"""Previews and gallery of a review folder."""
import html
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

CONTESTED = 4         # preview code of a contested point
COL = {0: (127, 127, 127), 3: (0, 150, 60), 1: (235, 104, 52), 2: (42, 120, 214), CONTESTED: (220, 0, 200)}
NAME = {0: "RCT tree", 3: "roots from both", 1: "v1 only", 2: "v2 only", CONTESTED: "contested (also in another tree)"}
CTX_COL = (215, 215, 215)


def preview(path, tree_xyz, code, n_points):
    rng = np.random.default_rng(0)
    k = rng.permutation(len(tree_xyz))[:n_points]
    x, s = tree_xyz[k], code[k]
    ext = np.ptp(x, 0)
    w = [ext[0], ext[1], ext[0]]
    fig, axs = plt.subplots(1, 3, figsize=(15, 5.5), gridspec_kw={"width_ratios": [max(v, 1) for v in w]})
    views = [(0, 2), (1, 2), (0, 1)]           # side (x-z), side at 90 deg (y-z), top (x-y)
    for a, (u, v) in zip(axs, views):
        for c in (0, 3, 1, 2, CONTESTED):
            m = s == c
            a.scatter(x[m, u], x[m, v], s=0.3 if c == 0 else 0.6, c=[np.array(COL[c]) / 255], lw=0, rasterized=True)
        a.set_aspect("equal"); a.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.03)
    fig.savefig(path, dpi=70); plt.close(fig)


def gallery(trees, seconds, out, title, header, rooted):
    """Gallery page; legend: the colours the previews use (without roots: RCT tree only)."""
    sec = []
    n = trees.genus.value_counts()
    order = [g for g in n.index if g != "unknown"] + (["unknown"] if "unknown" in n.index else [])
    for genus in order:
        g = trees[trees.genus == genus].sort_values("priority")
        cards = []
        for r in g.itertuples():
            roots = (f"<br>roots: both {r.n_both:,} &middot; v1 only {r.n_v1_only:,} &middot; v2 only {r.n_v2_only:,}"
                     f" &middot; contested {r.n_contested:,}" if rooted else "<br>own segment only")
            gsrc = f'<br>genus: <b>{html.escape(genus)}</b> ({html.escape(r.genus_source)})'
            flags = f'<div class="flag">{html.escape(r.flags)}</div>' if isinstance(r.flags, str) and r.flags else ""
            d = f"{r.D_cm:.0f} cm" if np.isfinite(r.D_cm) else "?"
            cards.append(
                f'<div class="card"><a href="{r.preview}"><img loading="lazy" src="{r.preview}" '
                f'alt="{html.escape(r.tree)}"></a><div class="meta"><b>{r.priority:03d} {html.escape(r.tree)}</b>'
                f'<br>H {r.H_m:.1f} m &middot; D {d}{gsrc}{roots}{flags}</div></div>')
        sec.append(f"<h2>{html.escape(genus)} ({len(g)})</h2><div class='grid'>{''.join(cards)}</div>")
    sw = lambda c: f'<span class="sw" style="background:rgb{c}"></span>'
    legend = " ".join(f"{sw(COL[k])}{NAME[k]}" for k in ((0, 3, 1, 2, CONTESTED) if rooted else (0,)))
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>{title}</title>
<style>
:root {{ --bg:#fafaf9; --card:#ffffff; --ink:#1c1c1a; --muted:#6b6b66; --line:#e2e2de; --flag:#a8461b; }}
@media (prefers-color-scheme: dark) {{ :root {{ --bg:#161615; --card:#21211f; --ink:#ececea; --muted:#a4a49e;
  --line:#3a3a37; --flag:#f08a5d; }} }}
body {{ margin:0; padding:16px; background:var(--bg); color:var(--ink); font:14px/1.4 system-ui, sans-serif; }}
h1 {{ font-size:20px; margin:0 0 4px; }} h2 {{ font-size:17px; margin:24px 0 8px; }}
.legend, .note {{ color:var(--muted); margin:4px 0; }}
.sw {{ display:inline-block; width:12px; height:12px; border-radius:2px; margin:0 4px 0 10px; vertical-align:-1px; }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fill, minmax(340px, 1fr)); gap:12px; }}
.card {{ background:var(--card); border:1px solid var(--line); border-radius:6px; overflow:hidden; }}
.card img {{ width:100%; display:block; background:#fff; }}
.meta {{ padding:8px 10px; }} .flag {{ color:var(--flag); margin-top:4px; }}
</style></head><body>
<h1>{title}: {len(trees)} trees</h1>
{header}<div class="legend">{legend}</div>
<div class="note">Sorted by priority per genus (trees filling thin D x H classes first). Views: side, side at 90&deg;,
top. Fill in checklist.csv: status ok / fix / drop, and genus. For a tree marked fix: rayroots-plot context, open its .ply and
context/&lt;tree&gt;_context.ply together in CloudCompare (keep the proposed shift), save as trees/&lt;tree&gt;_clean.ply. Built {time.strftime('%Y-%m-%d %H:%M')} in {seconds:.0f} s.</div>
{''.join(sec)}
</body></html>"""
    out.write_text(page, encoding="utf-8")
