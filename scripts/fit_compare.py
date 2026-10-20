"""Compare two fit_spiral runs of the same band: where the fitted windings sit, and which turn each is assigned to.

Why not nearest-vertex distance: fit_spiral writes each winding on a grid of `output_step_size` (20 vx) vertices, so
the nearest *vertex* of an identical surface is already a median 5-8 vx away (the grid floor, printed below). This
script measures instead
  1. point-to-surface distance: for each vertex of OTHER, the smallest |normal offset| to BASE among its 8 nearest
     BASE vertices whose tangential offset is <= 15 vx (normals from the BASE grid);
  2. turn agreement: each vertex gets a continuous turn number t = winding + position along the winding's theta axis;
     for OTHER vertices within 4 vx of BASE, dt = t_BASE - t_OTHER. A seam-free comparison: a constant dt is a harmless
     relabel, a dt that changes with angle means the fits disagree on which sheet is which turn;
  3. the modal dt in each 45-degree sector (around the umbilicus) of the inner / middle / outer third of the windings,
     and the outer-third sector range (max - min of the sector modes; 0 for a pure relabel).
The nearest-vertex numbers of the first version of this note are kept as `nearest_vertex` for comparison.

Usage: python fit_compare.py OUT.json UMBILICUS.json NAME=BASE_FIT_DIR:OTHER_FIT_DIR [NAME=...]
(a FIT_DIR is a fit_spiral run directory containing meshes/fitted/w*; non-spliced windings are used; with
input_disable_patches the _spliced meshes are identical anyway.)"""
import glob, json, os, sys, numpy as np, tifffile
from scipy.spatial import cKDTree


def load(run_dir):
    ws = sorted(x for x in glob.glob(os.path.join(run_dir, "meshes", "fitted", "w*")) if not x.endswith("_spliced"))
    shapes = [tifffile.imread(os.path.join(w, "x.tif")).shape for w in ws]
    ax = 0 if len({s[0] for s in shapes}) > len({s[1] for s in shapes}) else 1  # theta axis: its length varies by winding
    step = 1.0 / json.load(open(os.path.join(ws[0], "meta.json")))["scale"][0]
    allv, P, N, T, W = [], [], [], [], []
    for w in ws:
        X, Y, Z = (tifffile.imread(os.path.join(w, a + ".tif")).astype(np.float64) for a in "xyz")
        ok = (X > 0) & (Y > 0) & (Z > 0) & np.isfinite(X) & np.isfinite(Y) & np.isfinite(Z)
        if not ok.any():
            continue
        allv.append(np.stack([Z[ok], Y[ok], X[ok]], 1))
        C = np.stack([Z, Y, X], -1); C[~ok] = np.nan
        n = np.cross(np.gradient(C, axis=0), np.gradient(C, axis=1)); n /= np.linalg.norm(n, axis=-1, keepdims=True) + 1e-12
        good = ok & np.isfinite(n).all(-1)
        wi = int(os.path.basename(w)[1:4]); frac = (np.indices(ok.shape)[ax] + 0.5) / ok.shape[ax]
        P.append(C[good]); N.append(n[good]); T.append(wi + frac[good]); W.append(np.full(good.sum(), wi))
    return dict(all=np.concatenate(allv), P=np.concatenate(P), N=np.concatenate(N), T=np.concatenate(T),
                W=np.concatenate(W), step=step)


def surface_match(B, Q, k=8, tmax=15.0):
    d, idx = cKDTree(B["P"]).query(Q, k=k)
    v = Q[:, None, :] - B["P"][idx]; dn = np.abs((v * B["N"][idx]).sum(-1))
    dt = np.sqrt(np.maximum((v ** 2).sum(-1) - dn ** 2, 0)); dn = np.where(dt <= tmax, dn, np.inf)
    j = np.argmin(dn, 1); r = np.arange(len(Q)); best = dn[r, j]
    return np.where(np.isfinite(best), best, d[:, 0]), idx[r, j]


def dist_stats(x):
    return dict(median=round(float(np.median(x)), 2), within2=round(float((x <= 2).mean()), 3),
                within4=round(float((x <= 4).mean()), 3), p90=round(float(np.quantile(x, .9)), 2))


def mode(a):
    u, c = np.unique(a, return_counts=True); k = np.argmax(c)
    return int(u[k]), round(float(c[k] / len(a)), 3)


def grid_floor(step, n=200000, seed=0):
    """Nearest-vertex distance between two IDENTICAL surfaces whose vertex grids are offset at random."""
    r = np.random.default_rng(seed)
    p2 = r.uniform(0, step, (n, 2)); d2 = np.linalg.norm(np.minimum(p2, step - p2), axis=1)
    p1 = r.uniform(0, step, n); d1 = np.minimum(p1, step - p1)
    return {"2d_offset": dist_stats(d2), "1d_offset": dist_stats(d1)}


def compare(base_dir, other_dir, umb):
    B, Q = load(base_dir), load(other_dir)
    nv, _ = cKDTree(B["all"]).query(Q["all"])
    ds, j = surface_match(B, Q["P"])
    on = ds <= 4
    dt = np.round(B["T"][j[on]] - Q["T"][on]).astype(int); wb = B["W"][j[on]]; P = Q["P"][on]
    ang = np.degrees(np.arctan2(P[:, 1] - np.interp(P[:, 0], umb[0], umb[1]), P[:, 2] - np.interp(P[:, 0], umb[0], umb[2]))) % 360
    lo, hi = int(B["W"].min()), int(B["W"].max()) + 1; cut = [lo + (hi - lo) / 3, lo + 2 * (hi - lo) / 3]
    bands = {"inner": wb < cut[0], "middle": (wb >= cut[0]) & (wb < cut[1]), "outer": wb >= cut[1]}
    sectors = {}
    for bn, m in bands.items():
        sectors[bn] = {}
        for s in range(8):
            mm = m & (ang >= 45 * s) & (ang < 45 * (s + 1))
            if mm.sum() >= 500:
                sectors[bn][f"{45 * s:03d}"] = mode(dt[mm])[0]
    om = list(sectors["outer"].values())
    m0, share = mode(dt)
    u, c = np.unique(dt, return_counts=True); top = {int(u[i]): round(float(c[i] / len(dt)), 3) for i in np.argsort(-c)[:5]}
    return dict(nearest_vertex=dist_stats(nv), point_to_surface=dist_stats(ds), grid_step_vx=B["step"],
                turn_modal_offset=m0, turn_share_at_mode=share, turn_offsets_top=top,
                outer_windings=f"w{int(np.ceil(cut[1]))}-{hi - 1}", sector_modal_offset=sectors,
                outer_sector_range=(max(om) - min(om)) if om else None,
                n_other=int(len(Q["all"])), n_colocated=int(on.sum()))


if __name__ == "__main__":
    out_path, umb_path, pairs = sys.argv[1], sys.argv[2], sys.argv[3:]
    cp = sorted(json.load(open(umb_path))["control_points"], key=lambda c: c["z"])
    umb = tuple(np.array([c[k] for c in cp], float) for k in ("z", "y", "x"))
    res = json.load(open(out_path)) if os.path.exists(out_path) else {}
    for p in pairs:
        name, dirs = p.split("=", 1); b, o = dirs.split(":")
        res[name] = compare(b, o, umb)
        if "grid_floor" not in res:
            res["grid_floor"] = grid_floor(res[name]["grid_step_vx"])
        r = res[name]
        print(f"{name}: nearest-vertex {r['nearest_vertex']['median']} | point-to-surface {r['point_to_surface']['median']} "
              f"({r['point_to_surface']['within4']} <= 4 vx) | turn offset {r['turn_modal_offset']:+d} on {r['turn_share_at_mode']} "
              f"| outer sector range {r['outer_sector_range']}", flush=True)
        json.dump(res, open(out_path, "w"), indent=1)
