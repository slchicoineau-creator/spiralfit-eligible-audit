"""Winding count N for the band, by the package's frozen rule (work/r3/w2/wind_count.py, unchanged logic; portable paths):
"ceil(1.2 x max sheet crossings of the published m7 mask along 8 radial rays from the published umbilicus at the band's centre
slice. Rising edges (0->1) with >= 3 vx of background between sheets." Umbilicus = the mapped one (prereg 3.3).
Run from the repo root: .venv/bin/python -I work/spiral_accuracy/scripts/wind_count.py"""
import json, math, os, re, sys, urllib.request
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
os.environ.setdefault("SB_CACHE", os.path.join(ROOT, "work/spiral_accuracy/data/zio_cache"))
sys.path.insert(0, os.path.join(ROOT, "lib"))
from zio import RemoteArray  # noqa
SC, VID, ZC = "PHerc0139", "20250728140407", 6206
D = os.path.join(ROOT, f"work/spiral_accuracy/data/{SC}_{VID}")
cp = json.load(open(f"{D}/umbilicus.json"))["control_points"]
zs = np.array([c["z"] for c in cp]); ux = float(np.interp(ZC, zs, [c["x"] for c in cp])); uy = float(np.interp(ZC, zs, [c["y"] for c in cp]))
S3 = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
ej = json.load(open(f"{D}/tracks/{SC}_{VID}_surface_m7_L0_th0.2.extract.json"))
url = ej["array_url"].replace("s3://vesuvius-challenge-open-data", S3)
arr = RemoteArray(url.rsplit("/", 1)[0], url.rsplit("/", 1)[1])
sl = arr.read(ZC, ZC + 1, 0, arr.shape[1], 0, arr.shape[2], workers=48)[0] > 0
counts = []
for k in range(8):
    a = 2 * math.pi * k / 8; r = np.arange(0, 6000); x = (ux + r * math.cos(a)).astype(int); y = (uy + r * math.sin(a)).astype(int)
    ok = (x >= 0) & (x < sl.shape[1]) & (y >= 0) & (y < sl.shape[0]); v = sl[y[ok], x[ok]].astype(np.int8)
    n, gap = 0, 99
    for val in v:
        if val:
            if gap >= 3:
                n += 1
            gap = 0
        else:
            gap += 1
    counts.append(n)
N = math.ceil(1.2 * max(counts))
res = dict(scroll=SC, zc=ZC, umbilicus_xy=[ux, uy], crossings=counts, N=N, array_url=ej["array_url"])
json.dump(res, open(os.path.join(ROOT, "work/spiral_accuracy/out/wind_count.json"), "w"), indent=1)
print(json.dumps(res))
