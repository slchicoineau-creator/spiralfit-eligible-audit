"""Envelope-centroid baseline ENV: area centroid of the filled papyrus envelope (Otsu threshold, closing, holes filled,
largest component) at pyramid level 5, at each scored truth point's own z (rules in score.py); PHerc0139 on its
9.362 um volume at the affine-mapped z. Writes env_<scroll>.npy (z, x, y) and env_scores.json next to this file."""
import json, os, sys, numpy as np
from scipy import ndimage as ndi
H = os.path.dirname(os.path.abspath(__file__)); sys.path[:0] = [H, os.path.join(H, "..", "..", "lib")]
import score
from zio import RemoteArray, S3
VOL = {s: f"{s}/volumes/{v}-9.362um-1.2m-113keV-masked.zarr" for s, v in
       [("PHerc0125", "20250821151825"), ("PHerc0211", "20250821151803"), ("PHerc0826", "20250821151701"), ("PHerc0139", "20250728140407")]}
LV = 5


def otsu(v):
    h, _ = np.histogram(v, 256, (0, 256)); p = h / max(h.sum(), 1); w = np.cumsum(p); m = np.cumsum(p * np.arange(256))
    return int(np.argmax((m[-1] * w - m) ** 2 / (w * (1 - w) + 1e-12)))


def env_centroid(sl):
    nz = sl[sl > 0]
    if nz.size < 500:
        return None
    yy, xx = np.mgrid[-4:5, -4:5]
    env = ndi.binary_fill_holes(ndi.binary_closing(sl > otsu(nz), xx ** 2 + yy ** 2 <= 16))
    lab, n = ndi.label(env)
    if n == 0:
        return None
    ys, xs = np.nonzero(lab == np.argmax(ndi.sum(env, lab, range(1, n + 1))) + 1)
    return np.array([xs.mean(), ys.mean()])


out = {}
for sc in score.SCROLLS:
    ra = RemoteArray(f"{S3}/{VOL[sc]}", LV); f = 2 ** LV; rows = []
    for x, y, z in score.scored_points(sc):
        zl = int(round(z)) // f
        if zl < ra.shape[0]:
            c = env_centroid(ra.read(zl, zl + 1, 0, ra.shape[1], 0, ra.shape[2], workers=16)[0])
            if c is not None:
                e = (c + 0.5) * f - 0.5; rows.append((z, e[0], e[1]))
    A = np.array(rows); np.save(os.path.join(H, f"env_{sc}.npy"), A)
    out[sc] = score.score(sc, A); print(sc, out[sc], flush=True)
json.dump(out, open(os.path.join(H, "env_scores.json"), "w"), indent=1)
