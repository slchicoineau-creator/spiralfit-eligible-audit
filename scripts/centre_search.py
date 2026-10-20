"""Centre search (pre-registered 2026-09-30, 17_Plan_CentreSearch.md): does the fitter's own track agreement peak at
a centre position? PHerc0826 [8136, 9136), arm A. Reads satisfaction_metrics_fitted.json of every probe fit
(published centre with seeds 1-4, the tolerance-curve shifts, the centre-search shifts) and fits a quadratic surface.
Q1: the -x 250 signal replicates iff the seed-2 -x 250 fit >= bar (seed band max + range).
Q2: clear peak iff the quadratic is concave, its maximum lies within +-500 vx on both axes, and its value there >= bar.
PASS iff Q1 and Q2. Writes results/centre_search.json and figs/centre_search.png.
Usage: python centre_search.py FITS_ROOT"""
import glob, json, os, sys
import numpy as np
F = sys.argv[1]; REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); MM = 9.362e-3
PROBES = {  # tag: (dx, dy, seed, run dir relative to FITS_ROOT)
    "published_s1": (0, 0, 1, "full/PHerc0826_8136_A"), "published_s2": (0, 0, 2, "sens/PHerc0826_8136_R"),
    "published_s3": (0, 0, 3, "tol/PHerc0826_8136_seed3"), "published_s4": (0, 0, 4, "tol/PHerc0826_8136_seed4"),
    "x50": (50, 0, 1, "tol/PHerc0826_8136_x50"), "x100": (100, 0, 1, "tol/PHerc0826_8136_x100"),
    "x250": (250, 0, 1, "sens/PHerc0826_8136_S"), "x500": (500, 0, 1, "tol/PHerc0826_8136_x500"),
    "y50": (0, 50, 1, "tol/PHerc0826_8136_y50"), "y100": (0, 100, 1, "tol/PHerc0826_8136_y100"),
    "y250": (0, 250, 1, "tol/PHerc0826_8136_y250"), "y500": (0, 500, 1, "tol/PHerc0826_8136_y500"),
    "xm250": (-250, 0, 1, "tol/PHerc0826_8136_xm250"),
    "cs_xm250_s2": (-250, 0, 2, "tol/PHerc0826_8136_cs_xm250_s2"), "cs_xm100": (-100, 0, 1, "tol/PHerc0826_8136_cs_xm100"),
    "cs_xm500": (-500, 0, 1, "tol/PHerc0826_8136_cs_xm500"), "cs_ym100": (0, -100, 1, "tol/PHerc0826_8136_cs_ym100"),
    "cs_ym250": (0, -250, 1, "tol/PHerc0826_8136_cs_ym250"), "cs_ym500": (0, -500, 1, "tol/PHerc0826_8136_cs_ym500"),
    "cs_xm250_ym250": (-250, -250, 1, "tol/PHerc0826_8136_cs_xm250_ym250"),
    "cs_xm250_yp250": (-250, 250, 1, "tol/PHerc0826_8136_cs_xm250_yp250"),
}
rows = {}
for tag, (dx, dy, seed, p) in PROBES.items():
    f = glob.glob(os.path.join(F, p, "*", "satisfaction_metrics_fitted.json"))
    if not f:
        print(f"{tag}: missing"); continue
    s = json.load(open(f[0]))["summary"]
    rows[tag] = dict(dx=dx, dy=dy, seed=seed, points=round(s["satisfied_track_points_fraction"], 4),
                     tracks=round(s["satisfied_tracks_fraction"], 4))


def analyse(key):
    band = [r[key] for t, r in rows.items() if t.startswith("published_")]
    bar = max(band) + (max(band) - min(band))
    q1v = rows.get("cs_xm250_s2", {}).get(key)
    P = [r for r in rows.values() if abs(r["dx"]) <= 500 and abs(r["dy"]) <= 500]
    X = np.array([[r["dx"] * MM, r["dy"] * MM] for r in P]); y = np.array([r[key] for r in P])
    A = np.c_[np.ones(len(X)), X[:, 0], X[:, 1], X[:, 0] ** 2, X[:, 1] ** 2, X[:, 0] * X[:, 1]]
    a, b, c, d, e, f = np.linalg.lstsq(A, y, rcond=None)[0]
    H = np.array([[2 * d, f], [f, 2 * e]]); eig = np.linalg.eigvalsh(H)
    concave = bool((eig < 0).all())
    peak = np.linalg.solve(H, -np.array([b, c])) if abs(np.linalg.det(H)) > 1e-12 else np.array([np.nan, np.nan])
    val = float(a + b * peak[0] + c * peak[1] + d * peak[0] ** 2 + e * peak[1] ** 2 + f * peak[0] * peak[1])
    inside = bool(np.all(np.abs(peak) <= 500 * MM))
    resid = y - A @ np.array([a, b, c, d, e, f])
    q1 = q1v is not None and q1v >= bar; q2 = concave and inside and val >= bar
    return dict(seed_band=[min(band), max(band)], bar=round(bar, 4), q1_value=q1v, q1_replicates=bool(q1),
                quad_coef=[float(v) for v in (a, b, c, d, e, f)], hessian_eigs=[float(v) for v in eig], concave=concave,
                peak_mm=[round(float(v), 3) for v in peak], peak_inside=inside, peak_value=round(val, 4),
                rms_residual=round(float(np.sqrt(np.mean(resid ** 2))), 4), n_probes=len(P), q2_clear_peak=bool(q2),
                verdict="PASS" if (q1 and q2) else "FAIL")


out = dict(probes=rows, primary_points=analyse("points"), secondary_tracks=analyse("tracks"))
json.dump(out, open(os.path.join(REPO, "results", "centre_search.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "probes"}, indent=1))
try:
    import matplotlib
    matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    P = list(rows.values()); sc = ax.scatter([r["dx"] * MM for r in P], [r["dy"] * MM for r in P], c=[r["points"] for r in P],
                                             cmap="viridis", s=90, edgecolor="k")
    for r in P:
        ax.annotate(f"{r['points']:.3f}", (r["dx"] * MM, r["dy"] * MM), fontsize=6, xytext=(4, 4), textcoords="offset points")
    pk = out["primary_points"]["peak_mm"]
    if np.all(np.isfinite(pk)):
        ax.plot(*pk, "r*", ms=14, label=f"quadratic peak ({out['primary_points']['verdict']})"); ax.legend(fontsize=7)
    ax.axhline(0, c="#999", lw=.5); ax.axvline(0, c="#999", lw=.5); plt.colorbar(sc, label="satisfied track-point fraction")
    ax.set_xlabel("centre shift x (mm)"); ax.set_ylabel("centre shift y (mm)")
    ax.set_title("PHerc0826 [8136, 9136): track agreement vs centre position", fontsize=9)
    plt.tight_layout(); plt.savefig(os.path.join(REPO, "figs", "centre_search.png"), dpi=130); print("wrote figs/centre_search.png")
except ImportError:
    print("matplotlib not available: figure skipped")
