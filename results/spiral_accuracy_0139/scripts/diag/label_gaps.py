"""Diagnostic (after the smoke fit; labels only, no fit): are the label gaps w029->w031 and w044->w047 really 2 and 3 pitches?
Along the S0-a rays: radial gap between consecutive used labels / median one-step gap on the same ray. -> out/diag_label_gaps.json"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from acc_lib import D, ROOT, SCORE_Z, TWO_PI, load_umbilicus, load_mesh, lines, line_sA, truth_meshes, spiral_sense, seam_angles, theta_ref_from, BAND
umb = load_umbilicus(os.path.join(D, "PHerc0139_20250728140407/umbilicus.json"))
tm = truth_meshes("9um"); s, _ = spiral_sense(tm, umb); th_ref, _ = theta_ref_from(seam_angles(tm, umb, s, zr=BAND))
cache = {k: list(lines(load_mesh(d))) for k, d in tm.items()}
zs = np.linspace(SCORE_Z[0] + 25, SCORE_Z[1] - 25, 10); ratios = {}
for zj in zs:
    for phi in np.radians(np.arange(0, 360, 5)):
        cr = {}
        for k, L in cache.items():
            best = min((P for P, _ in L if abs(np.median(P[:, 0]) - zj) <= 10), key=lambda P: abs(np.median(P[:, 0]) - zj), default=None)
            if best is None: continue
            A, r = line_sA(best, umb, s); a0 = A.min(); A = A + th_ref + ((a0 - th_ref) % TWO_PI) - a0
            tgt = th_ref + ((s * phi - th_ref) % TWO_PI)
            for i in range(len(A) - 1):
                lo, hi = sorted((A[i], A[i + 1]))
                if lo <= tgt < hi and hi - lo < 0.2:
                    cr[k] = r[i] + (tgt - A[i]) / (A[i + 1] - A[i]) * (r[i + 1] - r[i]); break
        ks = sorted(cr)
        one = [cr[b] - cr[a] for a, b in zip(ks, ks[1:]) if b - a == 1]
        if len(one) < 5: continue
        p = np.median(one)
        for a, b in zip(ks, ks[1:]):
            ratios.setdefault(f"{a}->{b}", []).append((cr[b] - cr[a]) / p)
out = {k: dict(n=len(v), median_ratio=round(float(np.median(v)), 3), q25=round(float(np.quantile(v, .25)), 3), q75=round(float(np.quantile(v, .75)), 3))
       for k, v in sorted(ratios.items(), key=lambda kv: int(kv[0].split("->")[0])) if len(v) >= 20}
json.dump(out, open(os.path.join(ROOT, "work/spiral_accuracy/out/diag_label_gaps.json"), "w"), indent=1)
for k, v in out.items(): print(k, v)
