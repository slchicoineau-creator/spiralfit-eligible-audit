"""Pre-registered winding count N for a pilot band: ceil(1.2 x max sheet crossings of the published m7 mask along 8
radial rays from the published umbilicus at the band's centre slice). Rising edges (0->1) with >= 3 vx of background
between sheets. Runs in WSL (shares the band builder's mask cache). Usage: python wind_count.py SCROLL VOLID ZC"""
import json, math, os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
os.environ.setdefault("SB_CACHE", os.path.join(os.environ.get("R3", os.path.expanduser("~/r3")), "/maskcache".lstrip("/")))
from zio import RemoteArray  # noqa
SC, VID, ZC = sys.argv[1], sys.argv[2], int(sys.argv[3])
D = os.path.join(os.environ.get("DATA", os.path.expanduser("~/data")), f"{SC}_{VID}")
cp = json.load(open(f"{D}/umbilicus.json"))["control_points"]
zs = np.array([c["z"] for c in cp]); ux = float(np.interp(ZC, zs, [c["x"] for c in cp])); uy = float(np.interp(ZC, zs, [c["y"] for c in cp]))
S3 = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
import re, urllib.request
base = f"https://dl.ash2txt.org/datasets/spiral_datasets/{SC}/{VID}/tracks/"
ej = re.findall(r'href="([^"]*extract\.json)"', urllib.request.urlopen(base, timeout=60).read().decode())[0]
url = json.loads(urllib.request.urlopen(base + ej, timeout=60).read())["array_url"].replace("s3://vesuvius-challenge-open-data", S3)
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
print(json.dumps(dict(scroll=SC, zc=ZC, umbilicus_xy=[ux, uy], crossings=counts, N=N)))
