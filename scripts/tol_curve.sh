#!/bin/bash
# Umbilicus tolerance curve (pre-registered 2026-09-29 19:36 CDT). PHerc0826 [8136, 9136), arm A, N = 110.
# Each fit = the umb_sens.sh / pilot_fit.sh recipe; only the umbilicus shift (dx, dy) or the seed changes.
# Usage: tol_curve.sh            (launches 2 lanes and waits)
#        tol_curve.sh one TAG DX DY SEED
SC=PHerc0826; VID=20250821151701; Z0=8136; Z1=9136; N=110
BASE=${DATA:-$HOME/data}/${SC}_${VID}; OUT=${R3:-$HOME/r3}/fits/tol

one() {
  TAG=$1; DX=$2; DY=$3; SEED=$4
  D=${DATA:-$HOME/data}/${SC}_${VID}_tol_$TAG; mkdir -p $D $OUT
  ln -sfn $BASE/lasagna_inputs $D/lasagna_inputs; ln -sfn ${R3:-$HOME/r3}/bands/${SC}_${Z0}/A $D/tracks; cp $BASE/spiral-scroll.json $D/
  python3 - "$BASE/umbilicus.json" "$D/umbilicus.json" $DX $DY <<'PY'
import json, sys
u = json.load(open(sys.argv[1]))
for c in u["control_points"]:
    c["x"] += int(sys.argv[3]); c["y"] += int(sys.argv[4])
json.dump(u, open(sys.argv[2], "w"))
PY
  export FIT_SPIRAL_OUT_DIR=$OUT/${SC}_${Z0}_$TAG; rm -rf $FIT_SPIRAL_OUT_DIR; mkdir -p $FIT_SPIRAL_OUT_DIR
  export FIT_SPIRAL_CONFIG_OVERRIDES="{\"z_begin\": $Z0, \"z_end\": $Z1, \"optimizer_num_training_steps\": 30000,
   \"input_use_tracks\": true, \"input_disable_patches\": true, \"input_use_outer_shell\": false,
   \"dense_spacing_mode\": \"grad_mag\", \"loss_weight_dense_spacing\": 0, \"optimizer_random_seed\": $SEED,
   \"shell_outer_winding_idx\": $N, \"model_gap_expander_num_windings\": $N, \"model_gap_expander_capacity_windings\": $((N+4))}"
  cd ${R3:-$HOME/r3}/spiral-fitting
  for attempt in 1 2; do
    if .venv/bin/python fit_spiral.py --dataset $D > $FIT_SPIRAL_OUT_DIR.log 2>&1; then echo "TOL_FIT_DONE $TAG"; return; fi
    echo "TOL_FIT_FAIL $TAG attempt $attempt: $(grep -a -m1 -E 'Error' $FIT_SPIRAL_OUT_DIR.log)"
  done
}

if [ "$1" = one ]; then shift; one "$@"; exit; fi
mkdir -p $OUT
lane1() { one seed3 0 0 3; one y50 0 50 1; one y100 0 100 1; one y250 0 250 1; one y500 0 500 1; }
lane2() { one seed4 0 0 4; one x50 50 0 1; one x100 100 0 1; one x500 500 0 1; one xm250 -250 0 1; }
lane1 > $OUT/lane1.log 2>&1 &
sleep 120   # stagger the single-core loading phase
lane2 > $OUT/lane2.log 2>&1 &
wait
echo TOL_ALL_DONE > $OUT/ALL_DONE
