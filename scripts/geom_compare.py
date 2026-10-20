"""EXPLORATORY (not pre-registered; run after the full-run renders failed coverage): geometric distance of each arm's
fit from the arm-A fit on the same band (median nearest-vertex distance and share within 4 vx), next to the seed-noise
floor (A seed 2 vs A seed 1, 0826 band only). Answers: do the track changes move the fitted surfaces more than
re-running the fitter with another seed does?"""
import glob, json, os, numpy as np, tifffile
from scipy.spatial import cKDTree
F = os.path.join(os.environ.get("R3", os.path.expanduser("~/r3")), "/fits".lstrip("/"))


def verts(d):
    P = []
    for w in sorted(x for x in glob.glob(os.path.join(d, "meshes", "fitted", "w*")) if not x.endswith("_spliced")):
        X, Y, Z = (tifffile.imread(os.path.join(w, a + ".tif")).astype(np.float64) for a in "xyz")
        ok = (X > 0) & (Y > 0) & (Z > 0) & np.isfinite(X) & np.isfinite(Y) & np.isfinite(Z)
        P.append(np.stack([Z[ok], Y[ok], X[ok]], 1))
    return np.concatenate(P)


run = lambda p: glob.glob(os.path.join(F, p, "*/"))[0]
out = {}
for band, extra in (("PHerc0826_8136", {"A_seed2": "sens/PHerc0826_8136_R"}), ("PHerc0211_8328", {})):
    base = cKDTree(verts(run(f"full/{band}_A")))
    arms = {a: f"full/{band}_{a}" for a in ("Astar", "B", "C")}; arms.update(extra)
    for a, p in arms.items():
        d, _ = base.query(verts(run(p)))
        out[f"{band}:{a}"] = dict(median=round(float(np.median(d)), 2), within4=round(float((d <= 4).mean()), 3))
        print(band, a, out[f"{band}:{a}"], flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "geom_compare.json"), "w"), indent=1)
