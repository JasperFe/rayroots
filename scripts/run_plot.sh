#!/usr/bin/env bash
# Part A in one go: chains the rayroots-plot steps. Every step can also be run on its own (rayroots-plot -h).
# run_plot.sh <plot>.las|<own>_segmented.ply --project P [--plot NAME] [--boundary b.csv] [--mesh m.ply] [--no-roots]
#             [--plot-genus G] [--genus-hints h.csv] [--existing-trees e.csv] [--settings plot.yaml] [--config roots.yaml]
#             [--ground grid|mesh] [--keep]
# run_plot.sh --plots plots.csv --project P [--settings plot.yaml] [--keep]
#   plots.csv: columns plot, cloud, boundary, roots (yes/no), optional plot_genus; paths relative to plots.csv.
#   A failed plot is logged and skipped; the overview comes last.
# Rerun after a stop: finished steps are skipped. Roots exit 3 (memory guard) restarts from the tile cache, at most 3 times.
# Run it detached: setsid nohup bash run_plot.sh <plot>.las --project P > run.log 2>&1 &
set -euo pipefail
SELF=$(realpath "$0")
USAGE="usage: run_plot.sh <plot>.las|<own>_segmented.ply --project P [options] | run_plot.sh --plots plots.csv --project P"
INPUT=; PROJECT=; PLOT=; PLOTS=; ROOTS=yes; BOUND=(); MESH=(); CONFIG=(); REVIEW=(); SET=(); KEEP=()
[ $# -gt 0 ] || { echo "$USAGE"; exit 1; }
if [ "$1" != "--plots" ]; then INPUT=$(realpath "$1"); shift; fi
while [ $# -gt 0 ]; do
  case $1 in
    --boundary) BOUND=("$1" "$(realpath "$2")"); shift 2 ;;
    --mesh) MESH=("$1" "$(realpath "$2")"); shift 2 ;;
    --config) CONFIG+=("$1" "$(realpath "$2")"); shift 2 ;;
    --ground) CONFIG+=("$1" "$2"); shift 2 ;;
    --project) PROJECT=$(realpath -m "$2"); shift 2 ;;
    --plot) PLOT=$2; shift 2 ;;
    --plots) PLOTS=$(realpath "$2"); shift 2 ;;
    --no-roots) ROOTS=no; shift ;;
    --plot-genus) REVIEW+=("$1" "$2"); shift 2 ;;
    --genus-hints|--existing-trees) REVIEW+=("$1" "$(realpath "$2")"); shift 2 ;;
    --settings) SET=("$1" "$(realpath "$2")"); shift 2 ;;
    --keep) KEEP=(--keep); shift ;;
    *) echo "unknown option $1"; exit 1 ;;
  esac
done
[ -n "$PROJECT" ] || { echo "$USAGE"; exit 1; }
command -v rayroots-plot > /dev/null || { echo "rayroots-plot not on PATH (conda activate rayroots?)"; exit 1; }

if [ -n "$PLOTS" ]; then
  failed=0
  while IFS='|' read -r p cloud bnd roots pg; do
    args=("$cloud" --project "$PROJECT" --plot "$p" ${SET[@]+"${SET[@]}"} ${KEEP[@]+"${KEEP[@]}"})
    [ -n "$bnd" ] && args+=(--boundary "$bnd")
    [ "$roots" = no ] && args+=(--no-roots)
    [ -n "$pg" ] && args+=(--plot-genus "$pg")
    echo "=== $p"
    if ! bash "$SELF" "${args[@]}" < /dev/null; then echo "FAILED: $p"; failed=$((failed + 1)); fi
  done < <(python - "$PLOTS" <<'PY'
import sys
from pathlib import Path
import pandas as pd
f = Path(sys.argv[1])
t = pd.read_csv(f, dtype=str).fillna("")
miss = {"plot", "cloud", "boundary", "roots"} - set(t.columns)
if miss:
    sys.exit(f"{f}: missing columns {sorted(miss)}")
for r in t.to_dict("records"):
    path = lambda v: str((f.parent / v).resolve()) if v else ""
    roots = "no" if r["roots"].strip().lower() == "no" else "yes"
    print("|".join([r["plot"], path(r["cloud"]), path(r["boundary"]), roots, r.get("plot_genus", "")]))
PY
)
  rayroots-plot overview --project "$PROJECT" --plots "$PLOTS" ${SET[@]+"${SET[@]}"}
  echo "$failed plots failed"
  [ "$failed" = 0 ] && exit 0 || exit 1
fi

[ -n "$INPUT" ] || { echo "$USAGE"; exit 1; }
SEG=()
case $INPUT in
  *_segmented.ply) [ -n "$PLOT" ] || PLOT=$(basename "$INPUT" _segmented.ply); SEG=(--segmented "$INPUT") ;;
  *) [ -n "$PLOT" ] || PLOT=$(basename "$INPUT" .las)
     rayroots-plot segment "$INPUT" --plot "$PLOT" --project "$PROJECT" ${SET[@]+"${SET[@]}"} ${KEEP[@]+"${KEEP[@]}"} ;;
esac
COMMON=(--plot "$PLOT" --project "$PROJECT")
NO=()
if [ $ROOTS = yes ]; then
  attempt=1
  while true; do
    rc=0
    rayroots-plot roots "${COMMON[@]}" ${SEG[@]+"${SEG[@]}"} ${MESH[@]+"${MESH[@]}"} ${BOUND[@]+"${BOUND[@]}"} \
      ${CONFIG[@]+"${CONFIG[@]}"} ${KEEP[@]+"${KEEP[@]}"} || rc=$?
    [ $rc = 0 ] && break
    if [ $rc = 3 ] && [ $attempt -le 3 ]; then attempt=$((attempt + 1)); sleep 60; continue; fi
    echo "STOPPED: roots exit $rc (5 = failed tiles, see the WARNING lines; 2 = bad input)"; exit "$rc"
  done
else
  NO=(--no-roots)
fi
rayroots-plot select "${COMMON[@]}" ${SEG[@]+"${SEG[@]}"} ${BOUND[@]+"${BOUND[@]}"} ${SET[@]+"${SET[@]}"}
rayroots-plot review "${COMMON[@]}" ${NO[@]+"${NO[@]}"} ${REVIEW[@]+"${REVIEW[@]}"} ${SET[@]+"${SET[@]}"}
echo "done: $(date -Is)"
