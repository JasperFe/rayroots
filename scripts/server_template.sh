#!/usr/bin/env bash
# Template: run_plot.sh on a remote machine. Fill in the three names, then run the lines by hand from your PC.
# rayextract needs far more memory than the laptop has on a full plot; RCT and the rayroots env must be on the server.
HOST=user@server          # ssh host (add -p PORT to ssh and -P PORT to scp where needed)
DIR=/path/on/server       # work folder there
LAS=/path/to/<stem>.las   # here
STEM=$(basename "$LAS" .las)

# 1. Copy in, check the copy
scp "$LAS" "$HOST:$DIR/"
sha256sum "$LAS"
ssh "$HOST" "sha256sum $DIR/$STEM.las"        # must match

# 2. Run detached (survives the ssh session)
ssh "$HOST" "cd $DIR && conda run -n rayroots setsid nohup bash /path/to/rayroots/scripts/run_plot.sh $STEM.las > run.out 2>&1 &"

# 3. Later: finished when $STEM.rct.log ends with "done:"; run.out says why it stopped otherwise
ssh "$HOST" "tail -1 $DIR/$STEM.rct.log; tail -3 $DIR/run.out"

# 4. Copy back; check the targets do not exist first. The ray cloud $STEM.ply is not needed back.
scp "$HOST:$DIR/${STEM}_segmented.ply" "$HOST:$DIR/${STEM}_mesh.ply" "$HOST:$DIR/${STEM}_trees.txt" .
scp "$HOST:$DIR/${STEM}_segmented_roots.ply" "$HOST:$DIR/${STEM}_roots.npz" "$HOST:$DIR/${STEM}_roots.csv" .
scp "$HOST:$DIR/$STEM.rct.log" "$HOST:$DIR/$STEM.rayimport.log" .
scp -r "$HOST:$DIR/trees" .
