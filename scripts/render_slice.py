"""Render one CT z-slice (downsampled pyramid level) with a labelled full-resolution coordinate grid, for visual
umbilicus annotation. Optional zoom box (full-res x0 y0 size) renders at a finer level.
Usage: python render_slice.py SCROLL VOLID Z LEVEL OUT.png [X0 Y0 SIZE]"""
import os, re, sys, urllib.request, numpy as np
from PIL import Image, ImageDraw
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib")); os.environ.setdefault("SB_CACHE", os.path.join(REPO, "cache", "slices"))
from zio import RemoteArray  # noqa
S3 = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
SC, VID, Z, LV, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
lst = urllib.request.urlopen(f"{S3}/?list-type=2&delimiter=/&prefix={SC}/volumes/", timeout=60).read().decode()
vol = sorted([q for q in re.findall(r"<Prefix>([^<]*)</Prefix>", lst) if VID in q and q.endswith(".zarr/")], key=lambda q: "masked" not in q)[0].rstrip("/")
ra = RemoteArray(f"{S3}/{vol}", LV); f = 2 ** LV
if len(sys.argv) > 6:
    X0, Y0, SZ = (int(v) for v in sys.argv[6:9])
else:
    X0, Y0, SZ = 0, 0, max(ra.shape[1], ra.shape[2]) * f
x0, y0, n = X0 // f, Y0 // f, SZ // f
img = ra.read(Z // f, Z // f + 1, y0, min(y0 + n, ra.shape[1]), x0, min(x0 + n, ra.shape[2]), workers=32)[0].astype(np.float32)
lo, hi = np.percentile(img[img > 0], [1, 99.5]) if (img > 0).any() else (0, 1)
im = Image.fromarray(np.clip((img - lo) / (hi - lo + 1e-9) * 255, 0, 255).astype(np.uint8)).convert("RGB")
s = max(1, 1000 // max(im.size)); im = im.resize((im.size[0] * s, im.size[1] * s), Image.NEAREST) if s > 1 else im
sc = im.size[0] / img.shape[1]; d = ImageDraw.Draw(im)
step = 1000 if SZ > 4000 else (250 if SZ > 1000 else 50)
for v in range((X0 // step + 1) * step, X0 + SZ, step):
    px = (v - X0) / f * sc; d.line([(px, 0), (px, im.size[1])], fill=(255, 60, 60), width=1); d.text((px + 2, 2), str(v), fill=(255, 255, 0))
for v in range((Y0 // step + 1) * step, Y0 + SZ, step):
    py = (v - Y0) / f * sc; d.line([(0, py), (im.size[0], py)], fill=(60, 160, 255), width=1); d.text((2, py + 2), str(v), fill=(0, 255, 255))
d.text((4, im.size[1] - 14), f"{SC} z={Z} level {LV}; red lines = x, blue lines = y (full-res voxels)", fill=(255, 255, 255))
im.save(OUT); print(OUT, im.size, "full-res box", X0, Y0, SZ)
