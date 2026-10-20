"""Stage 0 (prereg sections 3.3 and 4): map the umbilicus, then gates S0-a (labels ordered + tau branch-consistent),
S0-c (w023 surrounds the umbilicus with > 60 vx clearance) and S0-d (catalog transform vs mesh-fitted affine <= 5 vx).
CPU only, labels only (no fit exists). Run from the repo root:
  .venv/bin/python -I work/spiral_accuracy/scripts/s0_labels.py"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from acc_lib import (D, ROOT, SCORE_Z, BAND, TWO_PI, USED, load_umbilicus, umb_at, load_mesh, lines, line_sA, raw_angles,
                     truth_meshes, spiral_sense, seam_angles, theta_ref_from, map_umbilicus)

OUT = os.path.join(ROOT, "work/spiral_accuracy/out")
DS = os.path.join(D, "PHerc0139_20250728140407")
os.makedirs(OUT, exist_ok=True); os.makedirs(DS, exist_ok=True)
cat = json.load(open(os.path.join(D, "catalog_0139.json")))
T = cat["properties"]["20260102150214"]["transforms"][0]
assert T["to_volume_id"] == "20250728140407"
M = np.array(T["transformation_matrix"], float)
res = dict(transform_source="catalog properties[20260102150214].transforms[0]", matrix=M.tolist())

# --- umbilicus mapped into the fit frame (the dataset root's umbilicus.json) ---
map_umbilicus(os.path.join(D, "umbilicus_published_2399.json"), M, os.path.join(DS, "umbilicus.json"))
umb = load_umbilicus(os.path.join(DS, "umbilicus.json"))

# --- S0-d: catalog vs the team's mesh-fitted affine on the band's umbilicus points ---
Ma = np.array(json.load(open(os.path.join(ROOT, "work/umb/affine_0139.json")))["M_xyz1_to_xyz"], float)
src = json.load(open(os.path.join(D, "umbilicus_published_2399.json")))["control_points"]
P = np.array([[c["x"], c["y"], c["z"], 1.0] for c in src])
qc, qa = P @ M.T, P @ Ma.T
inb = (qc[:, 2] >= BAND[0]) & (qc[:, 2] < BAND[1])
d = np.linalg.norm(qc - qa, axis=1)[inb]
res["S0_d"] = dict(n_points=int(inb.sum()), median_vx=round(float(np.median(d)), 2), max_vx=round(float(d.max()), 2),
                   PASS=bool(d.max() <= 5.0))

tm = truth_meshes("9um")
res["truth_wraps"] = sorted(tm)
assert len(tm) == 34, len(tm)
s, sv = spiral_sense(tm, umb)
res["spiral_sense"] = dict(s=s, **sv)
seams = seam_angles(tm, umb, s, zr=BAND)
theta_ref, gap = theta_ref_from(seams)
res["truth_seams"] = dict(n=int(len(seams)), theta_ref_deg=round(float(np.degrees(theta_ref)), 1), largest_gap_deg=round(float(np.degrees(gap)), 1))

# --- S0-a: rays every 5 deg on 10 slices; main crossings (tau - k < 1) of every used wrap's nearest line ---
zs = np.linspace(SCORE_Z[0] + 25, SCORE_Z[1] - 25, 10)
cache = {k: list(lines(load_mesh(dd))) for k, dd in tm.items()}
order_ok = order_n = br_ok = br_n = 0
n_cross = []
for zj in zs:
    for phi in np.radians(np.arange(0, 360, 5)):
        cr = []
        for k, L in cache.items():
            best = None
            for Pl, _ in L:
                mz = np.median(Pl[:, 0])
                if abs(mz - zj) <= 10 and (best is None or abs(mz - zj) < abs(np.median(best[:, 0]) - zj)):
                    best = Pl
            if best is None:
                continue
            A, r = line_sA(best, umb, s)
            a0 = np.min(A); A = A + theta_ref + ((a0 - theta_ref) % TWO_PI) - a0
            tgt = theta_ref + ((s * phi - theta_ref) % TWO_PI)          # main revolution only
            for i in range(len(A) - 1):
                lo, hi = sorted((A[i], A[i + 1]))
                if lo <= tgt < hi and hi - lo < 0.2:
                    f = (tgt - A[i]) / (A[i + 1] - A[i])
                    rr = r[i] + f * (r[i + 1] - r[i])
                    cr.append((rr, k, k + (tgt - theta_ref) / TWO_PI))
        cr.sort()
        n_cross.append(len(cr))
        for (r1, k1, t1), (r2, k2, t2) in zip(cr, cr[1:]):
            if k1 == k2:
                continue
            br_n += 1; br_ok += abs((t2 - t1) - (k2 - k1)) < 0.1
            if r2 - r1 >= 4:
                order_n += 1; order_ok += k2 > k1
res["S0_a"] = dict(rays=int(len(zs) * 72), median_crossings_per_ray=float(np.median(n_cross)),
                   order_pairs=order_n, order_share=round(order_ok / max(order_n, 1), 4),
                   branch_pairs=br_n, branch_share=round(br_ok / max(br_n, 1), 4))
res["S0_a"]["PASS"] = bool(res["S0_a"]["order_share"] >= 0.95 and res["S0_a"]["branch_share"] >= 0.95)

# --- S0-c: w023 surrounds the umbilicus in every 20-slice bin; clearance > 60 vx ---
m = load_mesh(tm[23]); C = m["C"].reshape(-1, 3); C = C[np.isfinite(C).all(1)]
bins = []
for z0 in range(BAND[0], BAND[1], 20):
    sel = C[(C[:, 0] >= z0) & (C[:, 0] < z0 + 20)]
    if len(sel) < 10:
        bins.append(dict(z0=z0, n=int(len(sel)), ok=False)); continue
    th, r = raw_angles(sel, umb)
    a = np.sort(th % TWO_PI); g = np.diff(np.concatenate([a, [a[0] + TWO_PI]])).max()
    bins.append(dict(z0=z0, n=int(len(sel)), max_gap_deg=round(float(np.degrees(g)), 1), min_r_vx=round(float(r.min()), 1),
                     ok=bool(np.degrees(g) < 45 and r.min() > 60)))
res["S0_c"] = dict(bins=len(bins), bins_ok=int(sum(b["ok"] for b in bins)),
                   worst_gap_deg=max(b.get("max_gap_deg", 360) for b in bins),
                   min_r_vx=min(b.get("min_r_vx", 0) for b in bins))
res["S0_c"]["PASS"] = res["S0_c"]["bins_ok"] == res["S0_c"]["bins"]
uy, ux = umb_at(umb, np.array([6206.0]))
res["umbilicus_at_band_centre_xy"] = [round(float(ux[0]), 1), round(float(uy[0]), 1)]
json.dump(res, open(os.path.join(OUT, "s0_labels.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k not in ("matrix", "truth_wraps")}, indent=1))
