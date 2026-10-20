#!/bin/bash
# Build the three pilot arm track sets for one scroll band, then crossings + packed stores per arm (WSL).
# Usage: run_band.sh SCROLL VOLID Z0 Z1
W="$(cd "$(dirname "$0")" && pwd)"   # this scripts/ folder
set -e
SC=$1; VID=$2; Z0=$3; Z1=$4
OUT=${R3:-$HOME/r3}/bands/${SC}_${Z0}
cd ${R3:-$HOME/r3}/spiral-fitting
.venv/bin/python $W/band_tracks.py $SC $VID $Z0 $Z1 $OUT
for A in A Astar B C; do
  DBM=$OUT/$A/${SC}_${VID}_surface_m7_L0_th0.2.dbm
  .venv/bin/python convert_track_store.py $DBM --force --quiet
  .venv/bin/python build_track_crossings.py $DBM --z-min $((Z0-64)) --z-max $((Z1+64)) --force
done
cp $OUT/band_summary.json $W/../results/band_summary_${SC}_${Z0}.json
echo BAND_DONE $SC
