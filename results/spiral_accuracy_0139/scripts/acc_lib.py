"""Shared helpers for SPIRAL-ACCURACY-0139 (prereg/2026-10-13_partner_spiral_accuracy.md). Paths are relative to the repo root.

Meshes are tifxyz grids (x.tif, y.tif, z.tif; invalid <= 0). Each mesh carries one integer label: a truth wrap wNNN (k = NNN) or a
fit winding wNNN (W = NNN). The turn coordinate tau (prereg section 5) is computed the same way for both."""
import glob, json, os, re
import numpy as np
import tifffile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
D = os.path.join(ROOT, "work/spiral_accuracy/data")
USED = [k for k in range(23, 60) if k not in (30, 45, 46)]        # 34 wraps (prereg section 2)
BAND = (5706, 6706)
SCORE_Z = (5756, 6656)
TWO_PI = 2 * np.pi


def load_umbilicus(path):
    cp = sorted(json.load(open(path))["control_points"], key=lambda c: c["z"])
    z = np.array([c["z"] for c in cp], float)
    return z, np.array([c["y"] for c in cp], float), np.array([c["x"] for c in cp], float)


def umb_at(umb, z):
    return np.interp(z, umb[0], umb[1]), np.interp(z, umb[0], umb[2])     # y, x


def load_mesh(d):
    """-> dict(C=(R,C,3) zyx with NaN invalid, ax=around-axis index (0 rows / 1 cols))"""
    X, Y, Z = (tifffile.imread(os.path.join(d, a + ".tif")).astype(np.float64) for a in "xyz")
    ok = (X > 0) & (Y > 0) & (Z > 0) & np.isfinite(X) & np.isfinite(Y) & np.isfinite(Z)
    C = np.stack([Z, Y, X], -1); C[~ok] = np.nan
    dz0 = np.nanmedian(np.abs(np.diff(C[..., 0], axis=0)))
    dz1 = np.nanmedian(np.abs(np.diff(C[..., 0], axis=1)))
    ax = 1 if dz0 > dz1 else 0       # around-the-scroll axis = the one along which z changes least
    return dict(C=C, ax=ax)


def truth_meshes(frame="9um"):
    out = {}
    for d in sorted(glob.glob(os.path.join(D, f"labels_{frame}", "w*"))):
        k = int(os.path.basename(d)[1:4])
        if k in USED and os.path.exists(os.path.join(d, "x.tif")):
            out[k] = d
    return out


def fit_meshes(run_dir):
    ws = sorted(x for x in glob.glob(os.path.join(run_dir, "meshes", "fitted", "w*")) if not x.endswith("_spliced"))
    return {int(re.match(r"w(\d+)", os.path.basename(w)).group(1)): w for w in ws}


def lines(m):
    """Iterate the mesh's lines along the around axis: yields (P (n,3) zyx of valid vertices in order, index tuple list)."""
    C, ax = m["C"], m["ax"]
    n_lines = C.shape[1 - ax]
    for i in range(n_lines):
        L = C[:, i] if ax == 0 else C[i]
        ok = np.isfinite(L).all(1)
        if ok.sum() >= 3:
            idx = np.nonzero(ok)[0]
            yield L[ok], [((j, i) if ax == 0 else (i, j)) for j in idx]


def raw_angles(P, umb):
    uy, ux = umb_at(umb, P[:, 0])
    return np.arctan2(P[:, 1] - uy, P[:, 2] - ux), np.hypot(P[:, 1] - uy, P[:, 2] - ux)


def line_sA(P, umb, s):
    th, r = raw_angles(P, umb)
    # unwrap along the line, but break at grid holes: consecutive valid vertices far apart are still unwrapped in order
    return s * np.unwrap(th), r


def spiral_sense(meshes, umb, zr=BAND):
    """s = +1 if going along a wrap with increasing atan2 angle moves OUTWARD, else -1.
    Main estimator (continuity): the end of a line at its largest raw angle is nearer to wrap k+1 (s = +1) or to wrap k-1
    (s = -1); majority over lines of wraps with both neighbours. Secondary: sign vote of the r-vs-angle slope per line."""
    from scipy.spatial import cKDTree
    V = {}
    for k, d in meshes.items():
        C = load_mesh(d)["C"].reshape(-1, 3); V[k] = C[np.isfinite(C).all(1)]
    trees = {k: cKDTree(v) for k, v in V.items()}
    cont, votes = [], []
    for k, d in meshes.items():
        for P, _ in lines(load_mesh(d)):
            if zr is not None and not (zr[0] <= np.median(P[:, 0]) < zr[1]):
                continue
            th, r = raw_angles(P, umb); A = np.unwrap(th)
            if A.max() - A.min() < 1.8 * np.pi:
                continue
            votes.append(np.sign(np.polyfit(A, r, 1)[0]))
            if k - 1 in trees and k + 1 in trees:
                Ep = P[np.argmax(A)]
                cont.append(1.0 if trees[k + 1].query(Ep)[0] < trees[k - 1].query(Ep)[0] else -1.0)
    c, v = np.array(cont), np.array(votes)
    return (1.0 if c.mean() > 0 else -1.0), dict(continuity_lines=int(len(c)), continuity_share_plus=round(float((c > 0).mean()), 4),
                                                 slope_lines=int(len(v)), slope_share_plus=round(float((v > 0).mean()), 4))


def seam_angles(meshes, umb, s, zr=None):
    """Seam (turn start) angle, mod 2 pi, of every line: min of s*A along the line."""
    out = []
    for d in meshes.values():
        m = load_mesh(d)
        for P, _ in lines(m):
            if zr is not None and not (zr[0] <= np.median(P[:, 0]) < zr[1]):
                continue
            A, _ = line_sA(P, umb, s)
            out.append(float(np.min(A) % TWO_PI))
    return np.array(out)


def theta_ref_from(seams):
    a = np.sort(np.asarray(seams) % TWO_PI)
    if len(a) == 0:
        return 0.0, TWO_PI
    gaps = np.diff(np.concatenate([a, [a[0] + TWO_PI]]))
    i = int(np.argmax(gaps))
    return float((a[i] + gaps[i] / 2) % TWO_PI), float(gaps[i])


def mesh_tau(d, label, umb, s, theta_ref, zr=None):
    """-> P (n,3) zyx, tau (n,), N (n,3) normals (zyx) for valid vertices of the mesh (optionally only z in zr)."""
    m = load_mesh(d); C = m["C"]
    nrm = np.cross(np.gradient(C, axis=0), np.gradient(C, axis=1))
    nrm /= np.linalg.norm(nrm, axis=-1, keepdims=True) + 1e-12
    Ps, Ts, Ns = [], [], []
    for P, ij in lines(m):
        A, _ = line_sA(P, umb, s)
        a0 = np.min(A)
        shift = theta_ref + ((a0 - theta_ref) % TWO_PI) - a0          # move the line's start into [theta_ref, theta_ref + 2 pi)
        tau = label + (A + shift - theta_ref) / TWO_PI
        N = np.array([nrm[i, j] for i, j in ij])
        keep = np.isfinite(N).all(1)
        if zr is not None:
            keep &= (P[:, 0] >= zr[0]) & (P[:, 0] < zr[1])
        Ps.append(P[keep]); Ts.append(tau[keep]); Ns.append(N[keep])
    if not Ps:
        return np.zeros((0, 3)), np.zeros(0), np.zeros((0, 3))
    return np.concatenate(Ps), np.concatenate(Ts), np.concatenate(Ns)


def map_umbilicus(src_json, matrix, out_json):
    """Map published control points (xyz, source frame) through a 3x4 matrix; keep score; sort by z."""
    M = np.asarray(matrix, float)
    raw = json.load(open(src_json))
    cps = []
    for c in raw["control_points"]:
        q = M @ np.array([c["x"], c["y"], c["z"], 1.0])
        cps.append(dict(x=float(q[0]), y=float(q[1]), z=float(q[2]), score=c.get("score", 100)))
    cps.sort(key=lambda c: c["z"])
    json.dump(dict(control_points=cps, metadata=dict(mapped_from=os.path.basename(src_json), matrix=M.tolist())),
              open(out_json, "w"), indent=1)
    return cps
