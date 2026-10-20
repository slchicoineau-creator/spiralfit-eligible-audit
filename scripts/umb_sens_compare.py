"""Umbilicus sensitivity comparison (pre-registered): distance from each valid vertex of fits S and R to the nearest
vertex of the baseline fit (non-spliced fitted windings). Usage: python umb_sens_compare.py BASE_DIR S_DIR R_DIR"""
import glob, json, os, sys, numpy as np, tifffile
from scipy.spatial import cKDTree


def verts(d):
    P = []
    for w in sorted(x for x in glob.glob(os.path.join(d, "meshes", "fitted", "w*")) if not x.endswith("_spliced")):
        X, Y, Z = (tifffile.imread(os.path.join(w, a + ".tif")).astype(np.float64) for a in "xyz")
        ok = (X > 0) & (Y > 0) & (Z > 0) & np.isfinite(X) & np.isfinite(Y) & np.isfinite(Z)
        P.append(np.stack([Z[ok], Y[ok], X[ok]], 1))
    return np.concatenate(P)


B = cKDTree(verts(sys.argv[1])); out = {}
for name, d in (("S_shifted", sys.argv[2]), ("R_seed2", sys.argv[3])):
    dist, _ = B.query(verts(d)); out[name] = dict(median=float(np.median(dist)), share_within_4=float((dist <= 4).mean()), n=len(dist))
    print(name, out[name])
r = out["S_shifted"]["share_within_4"] / max(out["R_seed2"]["share_within_4"], 1e-9)
out["ratio"] = r; out["reading"] = "insensitive" if r >= 0.9 else ("sensitive" if r < 0.5 else "moderate")
print(f"ratio S/R = {r:.3f} -> {out['reading']}")
json.dump(out, open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "umb_sens_result.json"), "w"), indent=1)
