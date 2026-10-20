"""Ink-map sensitivity to surface error on PHerc1447 (public data; §8 of this note). Usage: ink_surface_error.py infer | score

Data (put under DATA=./data_1447):
  - HF dataset YoussefMoNader/ink-8um-pherc1447-surfaces (CC BY-NC 4.0): w062/predictions/v8in.png, w062/tifxyz/{x,y,z}.tif
  - Vesuvius open data PHerc1447/segments/20250702235910-auto_grown_20250702235910292/: surface-volumes/8.64um-...-20250521151220.tifs/*.tif
    -> tierb/<segment>/layers/, mesh/tifxyz/{x,y,z}.tif -> tierb/<segment>/tifxyz/
  - model: hf download YoussefMoNader/ink-8um-v8in --local-dir models/v8in (MIT)
infer: v8-in on the auto-grown segment, 24-layer windows s = 0..7 of 31 (normal offset s - 3.5 vx), stride 32
score: (a) offset sweep on points where the auto surface is <= 1 vx from refined w062; (b) r by surface distance
"""
import os, sys, glob, json, numpy as np, tifffile, PIL.Image as I
from scipy.ndimage import gaussian_filter
I.MAX_IMAGE_PIXELS = None; H = os.environ.get("DATA", "./data_1447"); A = os.path.join(H, "tierb", "20250702235910-auto_grown_20250702235910292")


def up(g, k):
    import torch, torch.nn.functional as F
    t = torch.from_numpy(g.astype(np.float32))[None, None]
    return F.interpolate(t, scale_factor=k, mode="bilinear", align_corners=True)[0, 0].numpy()


def infer():
    import torch
    sys.path.insert(0, os.path.join(H, "models", "v8in"))
    from ink8um import InkDetector
    from ink8um.inference import predict_stack
    f = sorted(glob.glob(os.path.join(A, "layers", "*.tif"))); full = np.stack([np.clip(tifffile.imread(p), 0, 200).astype(np.uint8) for p in f], 2)
    m = InkDetector.from_pretrained(os.path.join(H, "models", "v8in")).cuda().eval()
    for s in range(8):
        o = os.path.join(H, "tierb", f"ob_s{s}.npy")
        if not os.path.exists(o):
            p = predict_stack(m, np.ascontiguousarray(full[:, :, s:s + 24]), reverse=False, stride=32, batch_size=8, num_workers=0)
            np.save(o, p.astype(np.float16)); print("s", s, float(p.mean()), flush=True); torch.cuda.empty_cache()


def matches():
    from scipy.spatial import cKDTree
    ax, ay, az = (tifffile.imread(os.path.join(A, "tifxyz", c + ".tif")) for c in "xyz"); va = (ax > 0) & (ay > 0) & (az > 0)
    X, Y, Z = (up(np.where(va, g, 0), 20) for g in (ax, ay, az)); V = up(va.astype(np.float32), 20) > 0.999
    rr, cc = np.nonzero(V); t = cKDTree(np.stack([X[rr, cc], Y[rr, cc], Z[rr, cc]], 1))
    hx, hy, hz = (tifffile.imread(os.path.join(H, "w062", "tifxyz", c + ".tif")) for c in "xyz"); hv = (hx > 0) & (hy > 0) & (hz > 0)
    gi, gj = np.nonzero(hv); d, k = t.query(np.stack([hx[gi, gj], hy[gi, gj], hz[gi, gj]], 1), distance_upper_bound=20); s = np.isfinite(d)
    return gi[s] * 4, gj[s] * 4, d[s], rr[k[s]], cc[k[s]]


def score():
    r, c, d, a1, a2 = matches()
    ref = gaussian_filter(np.asarray(I.open(next(p for p in (os.path.join(H, "w062", "predictions", "v8in.png"), os.path.join(H, "w062", "v8in.png")) if os.path.exists(p))).convert("L")).astype(np.float32) / 255, 4)[r, c]
    maps = {s: gaussian_filter(np.load(os.path.join(H, "tierb", f"ob_s{s}.npy")).astype(np.float32), 4)[a1, a2] for s in range(8)}
    near = d <= 1
    sweep = {f"{s - 3.5:+.1f}": round(float(np.corrcoef(ref[near & (maps[s] > 0) & (ref > 0)], maps[s][near & (maps[s] > 0) & (ref > 0)])[0, 1]), 3) for s in range(8)}
    bins = {}
    for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 5), (5, 8), (8, 12), (12, 20)):
        mk = (d >= lo) & (d < hi) & (maps[3] > 0) & (ref > 0)
        if mk.sum() > 500:
            bins[f"{lo}-{hi}"] = dict(n=int(mk.sum()), r=round(float(np.corrcoef(ref[mk], maps[3][mk])[0, 1]), 3))
    mx = max(sweep.values()); ends = min(sweep["-3.5"], sweep["+3.5"])
    reading = "offset tolerated" if all(v >= mx - 0.05 for v in sweep.values()) else "offset matters" if mx - ends >= 0.15 else "partial"
    res = dict(n_near=int(near.sum()), offset_sweep_r=sweep, reading=reading, distance_bins_s3=bins)
    json.dump(res, open("results/ink_vs_surface_error.json", "w"), indent=1); print(json.dumps(res, indent=1))


if __name__ == "__main__":
    {"infer": infer, "score": score}[sys.argv[1]]()
