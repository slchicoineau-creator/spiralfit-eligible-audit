#!/bin/bash
# H2 (descriptive, prereg section 6): the package's instrument lanes/A_geometry/tools/fit_compare.py, unchanged, on every base pair
# and every shifted fit vs every base fit (finished fits only). Angles around the unshifted mapped umbilicus. CPU only.
ROOT=$(cd "$(dirname "$0")/../../.." && pwd); W=$ROOT/work/spiral_accuracy; F=$W/fits
PAIRS=()
B=(); for b in base_s1 base_s2 base_s3; do [ -e $F/$b/DONE ] && B+=($b); done
for ((i=0;i<${#B[@]};i++)); do for ((j=i+1;j<${#B[@]};j++)); do PAIRS+=("${B[i]}|${B[j]}=$F/${B[i]}:$F/${B[j]}"); done; done
for o in $(ls $F); do case $o in base_*|smoke_*) continue;; esac; [ -e $F/$o/DONE ] || continue
  for b in "${B[@]}"; do PAIRS+=("$b|$o=$F/$b:$F/$o"); done; done
cd $ROOT && $ROOT/.venv/bin/python -I lanes/A_geometry/tools/fit_compare.py $W/out/h2_fit_compare.json \
  $W/data/PHerc0139_20250728140407/umbilicus.json "${PAIRS[@]}"
