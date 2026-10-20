#!/bin/bash
# Pilot fit for one scroll band and one arm (pre-registered recipe; only the tracks DBM differs between arms).
# Usage: pilot_fit.sh SCROLL VOLID Z0 Z1 N ARM
W="$(cd "$(dirname "$0")" && pwd)"   # this scripts/ folder
set -e
SC=$1; VID=$2; Z0=$3; Z1=$4; N=$5; ARM=$6
BASE=${DATA:-$HOME/data}/${SC}_${VID}
D=${DATA:-$HOME/data}/${SC}_${VID}_arm${ARM}
BAND=${R3:-$HOME/r3}/bands/${SC}_${Z0}/${ARM}
mkdir -p $D
ln -sfn $BASE/umbilicus.json $D/umbilicus.json
ln -sfn $BASE/lasagna_inputs $D/lasagna_inputs
ln -sfn $BAND $D/tracks
cp $BASE/spiral-scroll.json $D/spiral-scroll.json
cd ${R3:-$HOME/r3}/spiral-fitting
[ -e $BASE/lasagna_inputs/las_008_nx.ome.zarr.respool_g2_pair/meta.json ] || .venv/bin/python pack_resident_pools.py $BASE/lasagna_inputs --what normals --normal-group 2
export FIT_SPIRAL_OUT_DIR=${R3:-$HOME/r3}/fits/${FITROOT:-pilot}/${SC}_${Z0}_${ARM}
rm -rf $FIT_SPIRAL_OUT_DIR; mkdir -p $FIT_SPIRAL_OUT_DIR
export FIT_SPIRAL_CONFIG_OVERRIDES="{\"z_begin\": $Z0, \"z_end\": $Z1, \"optimizer_num_training_steps\": 30000,
 \"input_use_tracks\": true, \"input_disable_patches\": true, \"input_use_outer_shell\": false,
 \"dense_spacing_mode\": \"grad_mag\", \"loss_weight_dense_spacing\": 0, \"optimizer_random_seed\": ${SEED:-1},
 \"shell_outer_winding_idx\": $N, \"model_gap_expander_num_windings\": $N,
 \"model_gap_expander_capacity_windings\": $((N+4))}"
for attempt in 1 2; do
  if .venv/bin/python fit_spiral.py --dataset $D > $FIT_SPIRAL_OUT_DIR.log 2>&1; then break; fi
  echo "fit attempt $attempt failed: $(grep -a -m1 -E 'Error' $FIT_SPIRAL_OUT_DIR.log)"
done
grep -a -E "satisfied_tracks|satisfied_track_points|save_mesh fitted" $FIT_SPIRAL_OUT_DIR.log | tail -3
echo FIT_DONE $SC $ARM
