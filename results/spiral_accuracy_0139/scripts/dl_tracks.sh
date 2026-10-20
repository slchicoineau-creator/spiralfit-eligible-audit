#!/bin/bash
# Download the published m7 track store of one scroll volume (public, dl.ash2txt.org). Usage: dl_tracks.sh SCROLL VOLID
# Paths relative to the repo root (run from anywhere). Adapted from work/rocm_ports/scripts/dl_tracks.sh.
set -e
SC=$1; VID=$2
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
D=$ROOT/work/spiral_accuracy/data/${SC}_${VID}
mkdir -p $D/tracks; cd $D/tracks
U=https://dl.ash2txt.org/datasets/spiral_datasets/$SC/$VID/tracks
B=${SC}_${VID}_surface_m7_L0_th0.2
curl -sS -O "$U/$B.extract.json"
mkdir -p $B.dbm.vctracks
for f in header.bin metadata.json arclengths.f64 family_codes.i8 offsets.i64 source_ids.u64 tortuosities.f64 z_bounds.i32 coordinates.i32; do
  curl -sS -C - --retry 10 -o $B.dbm.vctracks/$f "$U/$B.dbm.vctracks/$f"
done
curl -sS -C - --retry 10 -O "$U/$B.dbm.crossings.npz"
curl -sS -C - --retry 10 -O "$U/$B.dbm"
echo DL_DONE
