#!/bin/bash
# villa's vc_gen_umbilicus (post-#1823) run densely: every 50 z, 21 repeats per slice (coordinatewise median), seed 1.
# Needs a volume-cartographer build with the vc_gen_umbilicus target (here: the vc3d-deps Docker image with the build in
# the Docker volume vc3d_build at /build/out) and the xy normal-grids fetched by fetch_grids.py into ./ng/<scroll>/.
cd "$(dirname "$0")"
IMG=ghcr.io/scrollprize/vc3d-deps/linux:sha-0c371b1d472c5281b703d65517e980d945da693f
for spec in PHerc0125:2150:20350 PHerc0211:1750:13200 PHerc0826:1950:16250 PHerc0139:1450:20400; do
  IFS=: read sc z0 z1 <<< "$spec"
  MSYS_NO_PATHCONV=1 docker run --rm -v vc3d_build:/build -v "$PWD:/w" $IMG /build/out/bin/vc_gen_umbilicus \
    -i /w/ng/$sc -o /w/strong_$sc.json --slices $z0:$z1:50 --seed 1 --repeats 21 --threads 16 > strong_$sc.log 2>&1
done
