#!/bin/bash
# Umbilicus sensitivity fits (pre-registered). Waits for full-run part 1 to finish, then runs S (shifted) and R (seed 2).
W="$(cd "$(dirname "$0")" && pwd)"   # this scripts/ folder

while pgrep -f "full_scroll.sh PHerc0826" > /dev/null; do sleep 60; done
SC=PHerc0826; VID=20250821151701; Z0=8136; Z1=9136; N=110
BASE=${DATA:-$HOME/data}/${SC}_${VID}; cd ${R3:-$HOME/r3}/spiral-fitting
for V in S R; do
  D=${DATA:-$HOME/data}/${SC}_${VID}_sens$V; mkdir -p $D
  ln -sfn $BASE/lasagna_inputs $D/lasagna_inputs; ln -sfn ${R3:-$HOME/r3}/bands/${SC}_${Z0}/A $D/tracks; cp $BASE/spiral-scroll.json $D/
  SEED=1; [ $V = R ] && SEED=2
  python3 - "$BASE/umbilicus.json" "$D/umbilicus.json" $V <<'PY'
import json, sys
u = json.load(open(sys.argv[1]))
if sys.argv[3] == "S":
    for c in u["control_points"]: c["x"] += 250
json.dump(u, open(sys.argv[2], "w"))
PY
  export FIT_SPIRAL_OUT_DIR=${R3:-$HOME/r3}/fits/sens/${SC}_${Z0}_$V; rm -rf $FIT_SPIRAL_OUT_DIR; mkdir -p $FIT_SPIRAL_OUT_DIR
  export FIT_SPIRAL_CONFIG_OVERRIDES="{\"z_begin\": $Z0, \"z_end\": $Z1, \"optimizer_num_training_steps\": 30000,
   \"input_use_tracks\": true, \"input_disable_patches\": true, \"input_use_outer_shell\": false,
   \"dense_spacing_mode\": \"grad_mag\", \"loss_weight_dense_spacing\": 0, \"optimizer_random_seed\": $SEED,
   \"shell_outer_winding_idx\": $N, \"model_gap_expander_num_windings\": $N, \"model_gap_expander_capacity_windings\": $((N+4))}"
  .venv/bin/python fit_spiral.py --dataset $D > $FIT_SPIRAL_OUT_DIR.log 2>&1 || .venv/bin/python fit_spiral.py --dataset $D > $FIT_SPIRAL_OUT_DIR.log 2>&1
  echo "SENS_FIT_DONE $V"
done
.venv/bin/python $W/umb_sens_compare.py \
  $(ls -d ${R3:-$HOME/r3}/fits/full/${SC}_${Z0}_A/*/ | head -1) $(ls -d ${R3:-$HOME/r3}/fits/sens/${SC}_${Z0}_S/*/ | head -1) $(ls -d ${R3:-$HOME/r3}/fits/sens/${SC}_${Z0}_R/*/ | head -1)
echo SENS_DONE
