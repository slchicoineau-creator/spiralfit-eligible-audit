"""Pilot renders (pre-registered): for one fitted run, 12 crops of 18x18 vertices (341x341 px at 1 vx/px, the GrowGate /
rater-quiz scale) chosen uniformly at random among fully valid windows of all non-spliced fitted windings
(meshes/fitted/w*), rng seed = scroll_index*10 + arm_index; CT from the scroll's 9.362 um volume (nearest voxel after
bilinear vertex interpolation, exactly as build_quiz_small.windows). Runs in WSL.
Usage: python render_fit.py FIT_DIR SCROLL VOLID SEED OUT_DIR"""
import glob, json, os, re, sys, urllib.request, numpy as np, tifffile
from PIL import Image
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
os.environ.setdefault("SB_CACHE", os.path.join(os.environ.get("R3", os.path.expanduser("~/r3")), "/ctcache".lstrip("/")))
from zio import RemoteArray  # noqa
FIT, SC, VID, SEED, OUT = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
S3 = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
lst = urllib.request.urlopen(f"{S3}/?list-type=2&delimiter=/&prefix={SC}/volumes/", timeout=60).read().decode()
vol = sorted([q for q in re.findall(r"<Prefix>([^<]*)</Prefix>", lst) if VID in q and q.endswith(".zarr/")], key=lambda q: "masked" not in q)[0].rstrip("/")
ra = RemoteArray(f"{S3}/{vol}", 0)
os.makedirs(OUT, exist_ok=True); rng = np.random.default_rng(SEED); cands = []
ws = sorted(w for w in glob.glob(os.path.join(FIT, "meshes", "fitted", "w*")) if not w.endswith("_spliced"))
grids = {}
for w in ws:
    X, Y, Z = (tifffile.imread(os.path.join(w, a + ".tif")).astype(np.float64) for a in "xyz")
    ok = (X > 0) & (Y > 0) & (Z > 0) & np.isfinite(X) & np.isfinite(Y) & np.isfinite(Z)
    if ok.shape[0] < 18 or ok.shape[1] < 18:
        continue
    good = ok[:-17, :-17].copy()
    for dr in range(18):
        for dc in range(18):
            good &= ok[dr:dr + good.shape[0], dc:dc + good.shape[1]]
    grids[w] = (X, Y, Z)
    cands += [(w, int(r), int(c)) for r, c in np.argwhere(good)]
NCROP = int(os.environ.get("NCROP", "12"))
pick = [cands[i] for i in rng.choice(len(cands), min(NCROP, len(cands)), replace=False)] if cands else []
rows = []
for k, (w, r0, c0) in enumerate(pick):
    X, Y, Z = grids[w]
    t = np.linspace(0, 17, 17 * 20 + 1); R, C = np.meshgrid(t, t, indexing="ij"); i0, j0 = np.floor(R).astype(int), np.floor(C).astype(int)
    i1, j1 = np.minimum(i0 + 1, 17), np.minimum(j0 + 1, 17); fr, fc = R - i0, C - j0
    bil = lambda G: (lambda g: g[i0, j0] * (1 - fr) * (1 - fc) + g[i1, j0] * fr * (1 - fc) + g[i0, j1] * (1 - fr) * fc + g[i1, j1] * fr * fc)(G[r0:r0 + 18, c0:c0 + 18])
    px, py, pz = (np.rint(bil(G)).astype(int) for G in (X, Y, Z))
    box = ra.read(pz.min(), pz.max() + 1, py.min(), py.max() + 1, px.min(), px.max() + 1, workers=16)
    img = box[pz - pz.min(), py - py.min(), px - px.min()].astype(np.uint8)
    Image.fromarray(img).save(os.path.join(OUT, f"crop{k:02d}.png"))
    rows.append(dict(k=k, winding=os.path.basename(w), r0=r0, c0=c0))
json.dump(dict(fit=FIT, scroll=SC, seed=SEED, n_candidates=len(cands), crops=rows), open(os.path.join(OUT, "crops.json"), "w"), indent=1)
print(f"{SC} seed {SEED}: {len(cands)} candidate windows over {len(grids)} windings, rendered {len(rows)} -> {OUT}")
