#!/bin/bash
# SPIRAL-ACCURACY-0139 fit driver v2 (two cards: one copy per card, per-fit lock) = the gpu-share tenant job (prereg/2026-10-13_partner_spiral_accuracy.md sections 3 and 8).
# Usage: run_fits.sh PLAN          PLAN = smoke | main   (fits in the frozen priority order; finished fits are skipped)
# Job contract (GPU_COORDINATION.md): SIGTERM -> the fit is stopped (its own autosave every 1,000 iterations is the checkpoint)
# and the script exits within ~60 s; on restart an unfinished fit resumes from checkpoint_fitted.ckpt; exit 0 only when done.
set -u
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
W=$ROOT/work/spiral_accuracy
PY=$ROOT/work/rocm_ports/spiral_venv_rocm/bin/python
FITTER=$ROOT/work/villa_main/spiral-fitting
BASE=$W/data/PHerc0139_20250728140407
Z0=5706; Z1=6706
N=$(python3 -c "import json;print(json.load(open('$W/out/wind_count.json'))['N'])")
CAP_H=27
export ROCR_VISIBLE_DEVICES=0,1 HIP_VISIBLE_DEVICES=${HIP_VISIBLE_DEVICES:-0} ROCBLAS_USE_HIPBLASLT=0
export WANDB_MODE=disabled FIT_SPIRAL_CACHE_DIR=$W/data/fitcache
PLAN=${1:-main}
CARD=$HIP_VISIBLE_DEVICES
mkdir -p $W/fits/locks
if [ "$PLAN" = smoke ]; then
  FITS="smoke_base_s1:0:0:1:4000"
else
  FITS="base_s1:0:0:1 y50_s1:0:50:1 x50_s1:50:0:1 base_s2:0:0:2 y25_s1:0:25:1 y10_s1:0:10:1 base_s3:0:0:3 \
y50_s2:0:50:2 y25_s2:0:25:2 y10_s2:0:10:2 x25_s1:25:0:1 x10_s1:10:0:1 x50_s2:50:0:2 x25_s2:25:0:2 x10_s2:10:0:2"
fi
CHILD=; T0=; TAG=
on_term() {
  echo "$(date -Is) SIGTERM: stopping fit $TAG pid $CHILD" >> $W/logs/driver.log
  [ -n "$T0" ] && echo "{\"tag\": \"$TAG\", \"start\": $T0, \"end\": $(date +%s), \"hours\": $(python3 -c "import time;print(round((time.time()-$T0)/3600,4))"), \"rc\": \"sigterm\"}" >> $W/logs/gpu_hours.jsonl
  if [ -n "$CHILD" ]; then kill -TERM $CHILD 2>/dev/null; for i in $(seq 60); do kill -0 $CHILD 2>/dev/null || break; sleep 1; done; kill -KILL $CHILD 2>/dev/null; fi
  exit 143
}
trap on_term TERM INT
gpu_hours() { python3 -c "
import json,os
p='$W/logs/gpu_hours.jsonl'; t=0
if os.path.exists(p):
  for l in open(p): t+=json.loads(l)['hours']
print(round(t,3))"; }

for spec in $FITS; do
  IFS=: read TAG DX DY SEED STEPS <<< "$spec"; STEPS=${STEPS:-30000}
  RUN=$W/fits/$TAG
  [ -e $RUN/DONE ] && continue
  if [ -d $W/fits/locks/$TAG ] && [ "$(cat $W/fits/locks/$TAG/card 2>/dev/null)" != "$CARD" ]; then continue; fi
  OTHER=0; for l in $W/fits/locks/*/; do [ -d "$l" ] || continue; t=$(basename $l); [ "$(cat $l/card)" != "$CARD" ] && [ ! -e $W/fits/$t/DONE ] && OTHER=$((OTHER+1)); done
  MINE=$([ -d $W/fits/locks/$TAG ] && echo 1 || echo 0)
  if [ $MINE = 0 ] && python3 -c "import sys; sys.exit(0 if $(gpu_hours) + 2.4 * $OTHER >= $CAP_H else 1)"; then
    echo "$(date -Is) cap $CAP_H GPU-h reached, $TAG not started" >> $W/logs/driver.log; continue
  fi
  if [ $MINE = 0 ]; then mkdir $W/fits/locks/$TAG 2>/dev/null || continue; echo $CARD > $W/fits/locks/$TAG/card; fi
  DS=$W/data/fit_ds/$TAG; mkdir -p $DS $RUN
  ln -sfn $BASE/lasagna_inputs $DS/lasagna_inputs; ln -sfn $BASE/tracks $DS/tracks; cp $BASE/spiral-scroll.json $DS/
  python3 - $BASE/umbilicus.json $DS/umbilicus.json $DX $DY <<'PY'
import json, sys
u = json.load(open(sys.argv[1]))
for c in u["control_points"]:
    c["x"] += float(sys.argv[3]); c["y"] += float(sys.argv[4])
u.setdefault("metadata", {})["spiral_accuracy_shift_xy_vx"] = [float(sys.argv[3]), float(sys.argv[4])]
json.dump(u, open(sys.argv[2], "w"))
PY
  export FIT_SPIRAL_RUN_DIR=$RUN
  unset FIT_SPIRAL_RESUME_PATH
  [ -e $RUN/checkpoint_fitted.ckpt ] && export FIT_SPIRAL_RESUME_PATH=$RUN/checkpoint_fitted.ckpt
  export FIT_SPIRAL_CONFIG_OVERRIDES="{\"z_begin\": $Z0, \"z_end\": $Z1, \"optimizer_num_training_steps\": $STEPS,
 \"input_use_tracks\": true, \"input_disable_patches\": true, \"input_use_outer_shell\": false,
 \"dense_spacing_mode\": \"grad_mag\", \"loss_weight_dense_spacing\": 0, \"optimizer_random_seed\": $SEED,
 \"shell_outer_winding_idx\": $N, \"model_gap_expander_num_windings\": $N, \"model_gap_expander_capacity_windings\": $((N+4))}"
  LOG=$W/logs/fit_$TAG.log; echo "$(date -Is) card $CARD starts $TAG" >> $W/logs/driver.log
  { echo "=== $(date -Is) start $TAG dx=$DX dy=$DY seed=$SEED steps=$STEPS N=$N resume=${FIT_SPIRAL_RESUME_PATH:-none} HIP=$HIP_VISIBLE_DEVICES"
    echo "dataset root:"; ls -la $DS; echo "overrides: $FIT_SPIRAL_CONFIG_OVERRIDES"; } >> $LOG
  T0=$(date +%s)
  ( cd $FITTER && exec $PY fit_spiral.py --dataset $DS ) >> $LOG 2>&1 &
  CHILD=$!
  wait $CHILD; RC=$?
  T1=$(date +%s); CHILD=
  echo "{\"tag\": \"$TAG\", \"start\": $T0, \"end\": $T1, \"hours\": $(python3 -c "print(round(($T1-$T0)/3600,4))"), \"rc\": $RC}" >> $W/logs/gpu_hours.jsonl
  echo "=== $(date -Is) end $TAG rc=$RC" >> $LOG; T0=
  if [ $RC -eq 0 ] && ls $RUN/meshes/fitted/w* >/dev/null 2>&1; then
    date -Is > $RUN/DONE
    rm -rf $DS
    ( cd $ROOT && $ROOT/.venv/bin/python -I $W/scripts/score_accuracy.py $W/out/accuracy.json $TAG=$RUN ) >> $W/logs/score.log 2>&1 || true
  else
    echo "$(date -Is) card $CARD $TAG failed rc=$RC" >> $W/logs/driver.log; exit 1
  fi
done
LIVE=0; for l in $W/fits/locks/*/; do [ -d "$l" ] || continue; [ -e $W/fits/$(basename $l)/DONE ] || LIVE=1; done
if [ $LIVE = 1 ]; then echo "$(date -Is) card $CARD: no fit left for this card; the other card finishes the plan" >> $W/logs/driver.log; exit 0; fi
echo "$(date -Is) plan $PLAN finished" >> $W/logs/driver.log
if [ "$PLAN" = main ]; then
  bash $W/scripts/h2.sh >> $W/logs/h2.log 2>&1 || true
  ( cd $ROOT && $ROOT/.venv/bin/python -I $W/scripts/summarize.py ) >> $W/logs/summarize.log 2>&1 || true
fi
exit 0
