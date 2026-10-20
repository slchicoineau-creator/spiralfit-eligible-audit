"""Pre-registered scoring (work log 2026-09-29 19:40 CDT). An estimate is a set of (z, x, y) samples in the 9.362 um
frame; it is linearly interpolated in z (the fitter's convention) to every scored truth point.
  score(sc, zxy)        -> dict(n, median_mm, share_le_094, median_mm_smooth_truth)
  gate(medians, strong) -> bool, per the rule: median <= 0.94 on >= 3 of 4 AND below the strong incumbent on >= 3 of 4.
CLI: score.py EST_PATTERN   (pattern with {sc}; .json control_points or .npy (N,3) z,x,y)"""
import json, os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import truth
SCROLLS = ["PHerc0125", "PHerc0211", "PHerc0826", "PHerc0139"]
EXCL = {"PHerc0125": [(10016, 10516)], "PHerc0139": [(14610, 15110)]}
ZMAX = {"PHerc0211": 13200}
MM = truth.VOX_UM / 1000


def scored_points(sc):
    P = truth.points(sc); keep = np.ones(len(P), bool)
    for a, b in EXCL.get(sc, []):
        keep &= ~((P[:, 2] >= a) & (P[:, 2] <= b))
    if sc in ZMAX:
        keep &= P[:, 2] < ZMAX[sc]
    return P[keep]


def smooth_truth(sc, P):
    """Secondary: running median of the (unexcluded, in-range) truth over +-2 mm of z."""
    A = truth.points(sc)
    if sc in ZMAX:
        A = A[A[:, 2] < ZMAX[sc]]
    w = 2.0 / MM; out = []
    for p in P:
        m = np.abs(A[:, 2] - p[2]) <= w; out.append(np.median(A[m, :2], 0))
    return np.array(out)


def score(sc, zxy):
    zxy = np.asarray(zxy, float); zxy = zxy[np.argsort(zxy[:, 0])]
    P = scored_points(sc)
    ex = np.interp(P[:, 2], zxy[:, 0], zxy[:, 1]); ey = np.interp(P[:, 2], zxy[:, 0], zxy[:, 2])
    E = np.c_[ex, ey]; e = np.linalg.norm(E - P[:, :2], axis=1) * MM
    es = np.linalg.norm(E - smooth_truth(sc, P), axis=1) * MM
    return dict(n=int(len(e)), median_mm=float(np.median(e)), share_le_094=float((e <= 0.94).mean()),
                median_mm_smooth_truth=float(np.median(es)))


def gate(med, strong):
    ok = sum(med[s] <= 0.94 for s in SCROLLS); better = sum(med[s] < strong[s] for s in SCROLLS)
    return dict(n_le_094=int(ok), n_better_than_strong=int(better), PASS=bool(ok >= 3 and better >= 3))


def load(path):
    if path.endswith(".npy"):
        return np.load(path)
    cp = json.load(open(path))["control_points"]
    return np.array([(p["z"], p["x"], p["y"]) for p in cp], float)


if __name__ == "__main__":
    for sc in SCROLLS:
        f = sys.argv[1].format(sc=sc)
        if os.path.exists(f):
            print(sc, {k: round(v, 3) if isinstance(v, float) else v for k, v in score(sc, load(f)).items()})
