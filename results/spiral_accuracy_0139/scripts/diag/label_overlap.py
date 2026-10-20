"""Diagnostic (labels only): for each pair of used wraps (k, k+1..k+3) the share of k's band vertices within 4 vx (normal distance,
fit_compare.surface_match) of the other wrap: two labels on one sheet would show a large share. -> out/diag_label_overlap.json"""
import importlib.util, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from acc_lib import D, ROOT, SCORE_Z, load_mesh, truth_meshes
spec = importlib.util.spec_from_file_location("fc", os.path.join(ROOT, "lanes/A_geometry/tools/fit_compare.py"))
fc = importlib.util.module_from_spec(spec); spec.loader.exec_module(fc)
tm = truth_meshes("9um"); S = {}
for k, d in tm.items():
    C = load_mesh(d)["C"]; n = np.cross(np.gradient(C, axis=0), np.gradient(C, axis=1)); n /= np.linalg.norm(n, axis=-1, keepdims=True) + 1e-12
    ok = np.isfinite(C).all(-1) & np.isfinite(n).all(-1) & (C[..., 0] >= SCORE_Z[0]) & (C[..., 0] < SCORE_Z[1])
    S[k] = dict(P=C[ok], N=n[ok])
out = {}
for k in sorted(S):
    for k2 in sorted(S):
        if 0 < k2 - k <= 3:
            ds, _ = fc.surface_match(S[k2], S[k]["P"])
            out[f"{k}->{k2}"] = dict(share_within4=round(float((ds <= 4).mean()), 4), median=round(float(np.median(ds)), 2))
json.dump(out, open(os.path.join(ROOT, "work/spiral_accuracy/out/diag_label_overlap.json"), "w"), indent=1)
for kk, v in out.items():
    if v["share_within4"] > 0.05: print(kk, v)
