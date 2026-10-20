#!/bin/bash
# Download the published reference umbilici (vesuvius-challenge-open-data, CC BY-NC 4.0) into ./truth/ next to this script.
cd "$(dirname "$0")" && mkdir -p truth && cd truth
B=https://vesuvius-challenge-open-data.s3.amazonaws.com
for k in PHerc0125/representations/umbilicus/20250821151825-umbilicus-20260808111524.json \
         PHerc0139/representations/umbilicus/20260102150214-umbilicus-20260826112529.json \
         PHerc0211/representations/umbilicus/20250821151803-umbilicus-20260808112626.json \
         PHerc0826/representations/umbilicus/20250821151701-umbilicus-20260808113303.json; do
  curl -sS -o "${k%%/*}_$(basename $k)" "$B/$k"
done
ls
