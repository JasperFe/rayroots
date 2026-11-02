"""overview across plots: A_plot/overview.html (one table); in work/: overview.csv, covered_trees.csv (finished trees,
usable as --existing-trees for a next review)."""
import html
import json
from pathlib import Path

import numpy as np
import pandas as pd

from rayroots.plot import settings


def class_names(bins):
    return [f"D {a:g}-{b:g}" if np.isfinite(b) else f"D {a:g}+" for a, b in zip(bins[:-1], bins[1:])]


def run(project=".", plots_csv=None, settings_path=None):
    a = settings.part(project)
    (a / "work").mkdir(parents=True, exist_ok=True)
    bins = settings.load(settings_path)["review"]["d_bins_cm"]
    names = class_names(bins)
    plots = sorted(p.name for p in (a / "review").glob("*") if (p / "checklist.csv").exists())
    if plots_csv:
        listed = pd.read_csv(plots_csv, dtype=str).fillna("")["plot"].tolist()
        plots = listed + [p for p in plots if p not in listed]
    rows, covered = [], []
    for plot in plots:
        rv, work = settings.review_dir(project, plot), settings.work_dir(project, plot)
        row = dict.fromkeys(["trees_selected", "trees_review", *names, "no D", "open", "ok", "fix", "drop", "finished"], 0)
        row = {"plot": plot, **row}
        if not (rv / "checklist.csv").exists():
            rows.append(row)
            continue
        row["trees_selected"] = json.load(open(work / "run.json"))["n_selected"]
        tr = pd.read_csv(work / "trees.csv")
        row["trees_review"] = len(tr)
        for i, c in enumerate(names):
            row[c] = int(((tr.D_cm >= bins[i]) & (tr.D_cm < bins[i + 1])).sum())
        row["no D"] = int(tr.D_cm.isna().sum())
        ck = pd.read_csv(rv / "checklist.csv", dtype=str).fillna("")
        st = ck.status.str.strip().str.lower()
        gone = (st == "") & ck.notes.str.contains("[not proposed any more]", regex=False)
        row["open"] = int(((st == "") & ~gone).sum())
        for s in ("ok", "fix", "drop"):
            row[s] = int((st == s).sum())
        log_p = work / "review_log.csv"
        if log_p.exists():
            log = pd.read_csv(log_p, dtype=str).fillna("")
            done = log.loc[log.action == "copied", "tree"].unique()
            row["finished"] = len(done)
            ck_ = ck.set_index("tree")
            have = tr[tr.tree.isin(done)]
            covered.append(pd.DataFrame({"plot": plot, "tree": have.tree, "genus": ck_.genus.reindex(have.tree).values,
                                         "H_m": have.H_m, "D_cm": have.D_cm}))
        rows.append(row)
    t = pd.DataFrame(rows)
    total = {"plot": "total", **{c: int(t[c].sum()) if len(t) else 0 for c in t.columns if c != "plot"}}
    t = pd.concat([t, pd.DataFrame([total])], ignore_index=True)
    t.to_csv(a / "work" / "overview.csv", index=False)
    cov = pd.concat(covered) if covered else pd.DataFrame(columns=["plot", "tree", "genus", "H_m", "D_cm"])
    cov.to_csv(a / "work" / "covered_trees.csv", index=False, float_format="%.3f")
    head = "".join(f"<th>{html.escape(c)}</th>" for c in t.columns)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(v))}</td>" for v in r) + "</tr>"
                   for r in t.itertuples(index=False))
    (a / "overview.html").write_text(f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>Plot overview</title><style>
:root {{ --bg:#fafaf9; --ink:#1c1c1a; --line:#e2e2de; }}
@media (prefers-color-scheme: dark) {{ :root {{ --bg:#161615; --ink:#ececea; --line:#3a3a37; }} }}
body {{ margin:0; padding:16px; background:var(--bg); color:var(--ink); font:14px/1.4 system-ui, sans-serif; }}
.wrap {{ overflow-x:auto; }} table {{ border-collapse:collapse; }}
th, td {{ border-bottom:1px solid var(--line); padding:4px 10px; text-align:right; white-space:nowrap; }}
th:first-child, td:first-child {{ text-align:left; }} tr:last-child td {{ font-weight:600; }}
</style></head><body><h1 style="font-size:20px">Plot overview</h1>
<div class="wrap"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div></body></html>
""", encoding="utf-8")
    print(t.to_string(index=False))
    return 0
