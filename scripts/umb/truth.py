"""Truth umbilici in the 9.362 um frame of each 9 um truth volume, and scoring helpers.
points(sc) -> (N,3) array of (x, y, z) full-res voxels in the 9.362 um frame, sorted by z.
at(sc, z) -> (x, y) linear interpolation (the fitter's convention), NaN outside the truth z range."""
import glob, json, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
VOL = {"PHerc0125": "20250821151825", "PHerc0211": "20250821151803", "PHerc0826": "20250821151701", "PHerc0139": "20250728140407"}
VOX_UM = 9.362


def _raw(sc):
    f = glob.glob(os.path.join(HERE, "truth", f"{sc}_*.json"))[0]
    return np.array([(p["x"], p["y"], p["z"]) for p in json.load(open(f))["control_points"]], float)


def points(sc):
    P = _raw(sc)
    if sc == "PHerc0139":                                     # 2.399 um frame -> 9.362 um frame
        M = np.array(json.load(open(os.path.join(HERE, "affine_0139.json")))["M_xyz1_to_xyz"])
        P = np.c_[P, np.ones(len(P))] @ M.T
    return P[np.argsort(P[:, 2])]


def at(sc, z):
    P = points(sc); z = np.asarray(z, float)
    x = np.interp(z, P[:, 2], P[:, 0], left=np.nan, right=np.nan); y = np.interp(z, P[:, 2], P[:, 1], left=np.nan, right=np.nan)
    return np.stack([x, y], -1)


def err_mm(sc, z, xy):
    """xy error (mm) of estimates xy (N,2) at slices z (N,) against the interpolated truth."""
    return np.linalg.norm(np.asarray(xy, float) - at(sc, z), axis=-1) * VOX_UM / 1000


if __name__ == "__main__":
    for sc in VOL:
        P = points(sc); d = np.diff(P[:, 2]); r = []
        for i in range(1, len(P) - 1):                        # roughness: distance of each point from the chord of its neighbours
            t = (P[i, 2] - P[i - 1, 2]) / (P[i + 1, 2] - P[i - 1, 2]); r.append(np.linalg.norm(P[i, :2] - (P[i - 1, :2] + t * (P[i + 1, :2] - P[i - 1, :2]))))
        print(f"{sc}: n={len(P)} z {P[0,2]:.0f}-{P[-1,2]:.0f} median dz {np.median(d):.0f}  chord-residual median {np.median(r):.1f} vx p90 {np.percentile(r,90):.1f} vx")
