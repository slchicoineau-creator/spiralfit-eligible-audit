"""Download published m7 xy normal-grid slices for a truth volume. Usage: fetch_grids.py SCROLL [truth|z0:z1:step]
'truth' = the slice nearest each truth control point (9.362 um frame). Output: ng/<SCROLL>/{metadata.json,xy/NNNNNN.grid}."""
import os, sys, urllib.request, numpy as np
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import truth
S3 = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
SC = sys.argv[1]; VID = truth.VOL[SC]
PRE = f"{S3}/{SC}/representations/predictions/surfaces/{VID}-surface-20260413222639-surface-m7-L0-th0.2.normal-grids"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ng", SC); os.makedirs(os.path.join(OUT, "xy"), exist_ok=True)
spec = sys.argv[2] if len(sys.argv) > 2 else "truth"
zs = sorted(set(np.rint(truth.points(SC)[:, 2]).astype(int))) if spec == "truth" else list(range(*map(int, spec.split(":"))))
if not os.path.exists(os.path.join(OUT, "metadata.json")):
    urllib.request.urlretrieve(f"{PRE}/metadata.json", os.path.join(OUT, "metadata.json"))


def get(z):
    p = os.path.join(OUT, "xy", f"{z:06d}.grid")
    if os.path.exists(p) and os.path.getsize(p) > 0:
        return 0
    for a in range(4):
        try:
            urllib.request.urlretrieve(f"{PRE}/xy/{z:06d}.grid", p + ".part"); os.replace(p + ".part", p); return 1
        except Exception as e:
            err = e
    print("FAIL", z, err, flush=True); return 0


with ThreadPoolExecutor(16) as ex:
    n = sum(ex.map(get, zs))
print(SC, "slices", len(zs), "downloaded", n, flush=True)
