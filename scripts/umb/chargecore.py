"""ChargeCore: umbilicus from the topological charge (Poincare index) of the sheet-orientation field of one xy slice's
published m7 normal-grid polylines. Spec: 04_Work_Log.md 2026-09-29 20:02 CDT (judge 15_Judge_Gate_UMB.md sect. C) plus
the 20:06 clarifications. Coordinates: full-res 9.362 um voxels (x, y). Parameters are fixed by the pre-registration."""
import numpy as np
from scipy import ndimage as ndi

CELL = 8            # vx per field cell (C3.1)
SIGMA = 48.0        # vx, Gaussian of the normalised convolution (C3.2)
VALID_FRAC = 0.02   # valid cell: W_s >= 0.02 * p90(W_s) (C3.3)
CLOSE_VX, ERODE_VX = 64, 96   # body mask (C3.4)
ANNIH_VX = 192.0    # (+1/2, -1/2) annihilation radius (C4.1)
CORE_VX = 16.0      # two surviving +1/2 closer than this = one core candidate (clarification 1)
PAIR_MAX_VX = 3000.0  # (C4.3)
DOT_MIN = 0.5       # (C4.3)
SAMPLE_VX = 4.0     # rasterisation step (clarification 3)


def _disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return x * x + y * y <= r * r


def rasterise(paths, W, H):
    """paths: list of (N,2) int arrays (x, y). Returns Qx, Qy, Wt on the cell grid (rows = y)."""
    nx, ny = int(np.ceil(W / CELL)), int(np.ceil(H / CELL))
    xs, ys, c2, s2, ws = [], [], [], [], []
    for P in paths:
        if len(P) < 2:
            continue
        P = P.astype(np.float64); d = np.diff(P, axis=0); L = np.hypot(d[:, 0], d[:, 1]); ok = L > 0
        P0, d, L = P[:-1][ok], d[ok], L[ok]
        phi = np.arctan2(d[:, 1], d[:, 0]); n = np.maximum(1, np.rint(L / SAMPLE_VX)).astype(int)
        rep = np.repeat(np.arange(len(L)), n); k = np.concatenate([(np.arange(m) + 0.5) / m for m in n])
        pts = P0[rep] + d[rep] * k[:, None]
        xs.append(pts[:, 0]); ys.append(pts[:, 1]); c2.append(np.cos(2 * phi[rep])); s2.append(np.sin(2 * phi[rep]))
        ws.append((L / n)[rep])
    if not xs:
        z = np.zeros((ny, nx)); return z, z.copy(), z.copy()
    xs, ys, c2, s2, ws = (np.concatenate(a) for a in (xs, ys, c2, s2, ws))
    ix = np.clip((xs // CELL).astype(int), 0, nx - 1); iy = np.clip((ys // CELL).astype(int), 0, ny - 1)
    flat = iy * nx + ix; sz = nx * ny
    Qx = np.bincount(flat, ws * c2, sz).reshape(ny, nx); Qy = np.bincount(flat, ws * s2, sz).reshape(ny, nx)
    Wt = np.bincount(flat, ws, sz).reshape(ny, nx)
    return Qx, Qy, Wt


def detect(paths, W, H, sigma=SIGMA, return_field=False):
    Qx, Qy, Wt = rasterise(paths, W, H)
    s = sigma / CELL
    Qxs, Qys, Ws = (ndi.gaussian_filter(a, s, mode="constant") for a in (Qx, Qy, Wt))
    th = np.arctan2(Qys, Qxs)                                         # = 2*phi
    valid = Ws >= VALID_FRAC * np.percentile(Ws, 90) if Ws.max() > 0 else np.zeros_like(Ws, bool)
    body = ndi.binary_closing(Wt > 0, _disk(CLOSE_VX // CELL))
    body = ndi.binary_fill_holes(body)
    body = ndi.binary_erosion(body, _disk(ERODE_VX // CELL))
    ok = valid & body
    # plaquette charge: loop (i,j) -> (i,j+1) -> (i+1,j+1) -> (i+1,j) -> (i,j)
    a, b, c, d = th[:-1, :-1], th[:-1, 1:], th[1:, 1:], th[1:, :-1]
    w = lambda x: (x + np.pi) % (2 * np.pi) - np.pi
    q = (w(b - a) + w(c - b) + w(d - c) + w(a - d)) / (4 * np.pi)
    okp = ok[:-1, :-1] & ok[:-1, 1:] & ok[1:, 1:] & ok[1:, :-1]
    k = np.rint(q * 2) / 2 * okp                                     # multiples of 1/2
    pi_, pj = np.nonzero(k != 0)
    plus, minus = [], []
    for i, j in zip(pi_, pj):
        xy = ((j + 1) * CELL, (i + 1) * CELL); m = int(round(abs(k[i, j]) * 2))
        (plus if k[i, j] > 0 else minus).extend([xy] * m)
    plus, minus = np.array(plus, float).reshape(-1, 2), np.array(minus, float).reshape(-1, 2)
    # C4.1 greedy nearest (+,-) annihilation within ANNIH_VX
    used_p, used_m = np.zeros(len(plus), bool), np.zeros(len(minus), bool)
    if len(plus) and len(minus):
        D = np.linalg.norm(plus[:, None] - minus[None], axis=-1)
        for f in np.argsort(D, axis=None):
            ip, im = np.unravel_index(f, D.shape)
            if D[ip, im] > ANNIH_VX:
                break
            if not used_p[ip] and not used_m[im]:
                used_p[ip] = used_m[im] = True
    surv = plus[~used_p]
    n_plus = int(len(surv))
    # clarification 1: merge +1/2 closer than CORE_VX into core candidates
    cores, halves, used = [], [], np.zeros(len(surv), bool)
    if len(surv) > 1:
        D = np.linalg.norm(surv[:, None] - surv[None], axis=-1); np.fill_diagonal(D, np.inf)
        for f in np.argsort(D, axis=None):
            i, j = np.unravel_index(f, D.shape)
            if D[i, j] >= CORE_VX:
                break
            if not used[i] and not used[j]:
                used[i] = used[j] = True; cores.append((surv[i] + surv[j]) / 2)
    halves = surv[~used]
    # C4.2 comet axis on the normalised director (clarification 2)
    nrm = np.hypot(Qxs, Qys) + 1e-12; cx, sx = Qxs / nrm, Qys / nrm
    dcx_dy, dcx_dx = np.gradient(cx); dsx_dy, dsx_dx = np.gradient(sx)
    px_f, py_f = dcx_dx + dsx_dy, dsx_dx - dcx_dy
    pol = []
    for x, y in halves:
        j, i = int(x // CELL), int(y // CELL)
        v = np.array([px_f[i - 1:i + 1, j - 1:j + 1].mean(), py_f[i - 1:i + 1, j - 1:j + 1].mean()])
        pol.append(v / (np.linalg.norm(v) + 1e-12))
    cands = [("core", c) for c in cores]
    for i in range(len(halves)):
        for j in range(i + 1, len(halves)):
            u = halves[j] - halves[i]; dist = np.linalg.norm(u)
            if dist > PAIR_MAX_VX or dist == 0:
                continue
            u = u / dist; ai, aj = pol[i] @ u, pol[j] @ (-u)
            if abs(ai) >= DOT_MIN and abs(aj) >= DOT_MIN and np.sign(ai) == np.sign(aj):
                cands.append(("pair", (halves[i] + halves[j]) / 2))
    by, bx = np.nonzero(body)
    cen = np.array([(bx.mean() + 0.5) * CELL, (by.mean() + 0.5) * CELL]) if len(bx) else np.array([W / 2, H / 2])
    est = None
    if cands:
        best = min(cands, key=lambda t: np.linalg.norm(t[1] - cen)); est = (float(best[1][0]), float(best[1][1]), best[0])
    out = dict(est=est, n_plus=n_plus, n_minus_surv=int((~used_m).sum()), n_cands=len(cands), n_cores=len(cores),
               body_centroid=cen.tolist(), plus_surv=surv.tolist())
    if return_field:
        out.update(theta=th, valid=valid, body=body, Ws=Ws, minus_surv=minus[~used_m].tolist())
    return out
