"""Rebuild the PHerc0139 2.399 um (20260102150214) -> 9.362 um (20250728140407) frame map.
Same-segment tifxyz meshes exist in both frames; vertices are paired by fractional grid index (bilinear in the
2.4 um grid), and a 3D affine (12 params) is fitted by iteratively trimmed least squares.
Acceptance (round 3): median residual about 4.1 vx at 9.362 um. Writes affine_0139.json."""
import io, json, os, urllib.request, numpy as np, tifffile
from scipy.ndimage import map_coordinates
S3 = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
HERE = os.path.dirname(os.path.abspath(__file__)); CACHE = os.path.join(HERE, "meshcache"); os.makedirs(CACHE, exist_ok=True)
SEGS = ["20250108000000-w025_2025010863", "20250108000003-w028_2025010836", "20250831000000-w040_2025083102",
        "20251226000000-w055_2025122611", "20250223000000-w059_2025022312"]
F9, F2 = "20250728140407-9.362um", "20260102150214-2.399um"


def grid(seg, frame):
    sid = seg.split("-")[0]; out = []
    for a in "xyz":
        p = os.path.join(CACHE, f"{sid}_{frame}_{a}.tif")
        if not os.path.exists(p):
            urllib.request.urlretrieve(f"{S3}/PHerc0139/segments/{seg}/mesh/{sid}-on-{frame}.tifxyz/{a}.tif", p)
        out.append(tifffile.imread(p).astype(np.float64))
    X, Y, Z = out
    return np.stack([X, Y, Z], -1), (X > 0) & (Y > 0) & (Z > 0)


P9, P2 = [], []
for s in SEGS:
    G9, V9 = grid(s, F9); G2, V2 = grid(s, F2)
    H9, W9 = V9.shape; H2, W2 = V2.shape
    ii, jj = np.nonzero(V9)
    fi, fj = ii * (H2 - 1) / (H9 - 1), jj * (W2 - 1) / (W9 - 1)
    ok = map_coordinates(V2.astype(np.float64), [fi, fj], order=1) > 0.999      # all 4 neighbours valid
    q2 = np.stack([map_coordinates(G2[..., k], [fi, fj], order=1) for k in range(3)], -1)
    P9.append(G9[ii, jj][ok]); P2.append(q2[ok]); print(s, "pairs", ok.sum(), "grids", V9.shape, V2.shape, flush=True)
P9, P2 = np.concatenate(P9), np.concatenate(P2)
A = np.c_[P2, np.ones(len(P2))]; keep = np.ones(len(P2), bool)
for it in range(5):
    M, *_ = np.linalg.lstsq(A[keep], P9[keep], rcond=None)
    r = np.linalg.norm(A @ M - P9, axis=1); keep = r < max(3 * np.median(r), 10)
r = np.linalg.norm(A @ M - P9, axis=1)
res = dict(frame_from=F2, frame_to=F9, M_xyz1_to_xyz=M.T.tolist(), n_pairs=int(len(r)), n_kept=int(keep.sum()),
           resid_median_vx=float(np.median(r)), resid_p95_vx=float(np.percentile(r, 95)),
           linear_part_scale=np.linalg.svd(M[:3].T, compute_uv=False).tolist())
json.dump(res, open(os.path.join(HERE, "affine_0139.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "M_xyz1_to_xyz"}, indent=1))
