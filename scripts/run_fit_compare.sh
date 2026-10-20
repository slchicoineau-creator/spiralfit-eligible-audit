#!/bin/bash
# Runs scripts/fit_compare.py on every pair of fits in this note -> results/fit_compare.json.
# FITS = the fit output root used by pilot_fit.sh / umb_sens.sh (pilot, seed2, seed3, full, sens); DATA = dataset root.
W="$(cd "$(dirname "$0")" && pwd)"; R="$(dirname "$W")/results"
FITS=${FITS:-${R3:-$HOME/r3}/fits}; DATA=${DATA:-$HOME/data}; PY=${PY:-python}
d() { ls -d $FITS/$1/*/ | head -1; }
U826=$DATA/PHerc0826_20250821151701/umbilicus.json; U211=$DATA/PHerc0211_20250821151803/umbilicus.json
OUT=$R/fit_compare.json
$PY $W/fit_compare.py $OUT $U826 \
  "0826_8136 seed1 vs seed2=$(d full/PHerc0826_8136_A):$(d sens/PHerc0826_8136_R)" \
  "0826_8136 umbilicus +250x=$(d full/PHerc0826_8136_A):$(d sens/PHerc0826_8136_S)" \
  "0826_8136 A vs Astar=$(d full/PHerc0826_8136_A):$(d full/PHerc0826_8136_Astar)" \
  "0826_8136 A vs B=$(d full/PHerc0826_8136_A):$(d full/PHerc0826_8136_B)" \
  "0826_8136 A vs C=$(d full/PHerc0826_8136_A):$(d full/PHerc0826_8136_C)" \
  "0826_10208 A seed1 vs seed2=$(d pilot/PHerc0826_10208_A):$(d seed2/PHerc0826_10208_A)" \
  "0826_10208 A seed1 vs seed3=$(d pilot/PHerc0826_10208_A):$(d seed3/PHerc0826_10208_A)" \
  "0826_10208 A seed2 vs seed3=$(d seed2/PHerc0826_10208_A):$(d seed3/PHerc0826_10208_A)" \
  "0826_10208 C seed1 vs seed2=$(d pilot/PHerc0826_10208_C):$(d seed2/PHerc0826_10208_C)" \
  "0826_10208 C seed1 vs seed3=$(d pilot/PHerc0826_10208_C):$(d seed3/PHerc0826_10208_C)" \
  "0826_10208 seed1 A vs C=$(d pilot/PHerc0826_10208_A):$(d pilot/PHerc0826_10208_C)"
$PY $W/fit_compare.py $OUT $U211 \
  "0211_8328 A vs Astar=$(d full/PHerc0211_8328_A):$(d full/PHerc0211_8328_Astar)" \
  "0211_8328 A vs B=$(d full/PHerc0211_8328_A):$(d full/PHerc0211_8328_B)" \
  "0211_8328 A vs C=$(d full/PHerc0211_8328_A):$(d full/PHerc0211_8328_C)" \
  "0211_10496 A seed1 vs seed2=$(d pilot/PHerc0211_10496_A):$(d seed2/PHerc0211_10496_A)" \
  "0211_10496 A seed1 vs seed3=$(d pilot/PHerc0211_10496_A):$(d seed3/PHerc0211_10496_A)" \
  "0211_10496 A seed2 vs seed3=$(d seed2/PHerc0211_10496_A):$(d seed3/PHerc0211_10496_A)" \
  "0211_10496 C seed1 vs seed2=$(d pilot/PHerc0211_10496_C):$(d seed2/PHerc0211_10496_C)" \
  "0211_10496 C seed1 vs seed3=$(d pilot/PHerc0211_10496_C):$(d seed3/PHerc0211_10496_C)" \
  "0211_10496 seed1 A vs C=$(d pilot/PHerc0211_10496_A):$(d pilot/PHerc0211_10496_C)"
echo FIT_COMPARE_DONE
