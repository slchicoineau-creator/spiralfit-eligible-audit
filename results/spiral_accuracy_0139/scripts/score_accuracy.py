"""SPIRAL-ACCURACY-0139 scorer (prereg/2026-10-13_partner_spiral_accuracy.md section 5). CPU only.

For each fit run directory: truth vertices of the 34 used wraps in the scoring range -> nearest fit surface with the package
instrument's `surface_match` (lanes/A_geometry/tools/fit_compare.py, unchanged) -> matched if <= 6 vx -> Delta = tau_fit - tau_truth,
c* = mode(round(Delta)), e = round(Delta) - c*. E = share of matched vertices with e != 0 (primary, lower is better).
Angles use the UNSHIFTED mapped umbilicus for every fit (the truth frame).
Usage (repo root): .venv/bin/python -I work/spiral_accuracy/scripts/score_accuracy.py OUT.json NAME=RUN_DIR [NAME=RUN_DIR ...]"""
import importlib.util, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from acc_lib import (D, ROOT, BAND, SCORE_Z, TWO_PI, USED, load_umbilicus, truth_meshes, fit_meshes, spiral_sense, seam_angles,
                     theta_ref_from, mesh_tau)

spec = importlib.util.spec_from_file_location("fit_compare", os.path.join(ROOT, "lanes/A_geometry/tools/fit_compare.py"))
fc = importlib.util.module_from_spec(spec); spec.loader.exec_module(fc)
UMB = os.path.join(D, "PHerc0139_20250728140407", "umbilicus.json")
TOL = 6.0
THIRDS = {"inner": (23, 34), "middle": (35, 46), "outer": (47, 59)}
DROP_CLEAN = (34, 35, 41, 42)     # amendment 1 (label-only consistency rule)


def breaks(modal):
    ks = sorted(modal); st = [f"{a}->{b}:{modal[b] - modal[a]:+d}" for a, b in zip(ks, ks[1:]) if modal[a] != modal[b]]
    return dict(n=len(st), steps=st)


def mode_int(a):
    u, c = np.unique(a, return_counts=True); return int(u[np.argmax(c)])


def score(run_dir, umb, s, truth, truth_seams):
    fm = fit_meshes(run_dir)
    fseams = seam_angles(fm, umb, s, zr=BAND)
    theta_ref, gap = theta_ref_from(np.concatenate([truth_seams, fseams]))
    # fit surface (all windings, band z only)
    FP, FT, FN = [], [], []
    for W, d in fm.items():
        P, T, N = mesh_tau(d, W, umb, s, theta_ref, zr=(BAND[0] - 40, BAND[1] + 40))
        FP.append(P); FT.append(T); FN.append(N)
    B = dict(P=np.concatenate(FP), N=np.concatenate(FN)); FT = np.concatenate(FT)
    TP, TT, TK = [], [], []
    for k, d in truth.items():
        P, T, _ = mesh_tau(d, k, umb, s, theta_ref, zr=SCORE_Z)
        TP.append(P); TT.append(T); TK.append(np.full(len(T), k))
    TP, TT, TK = np.concatenate(TP), np.concatenate(TT), np.concatenate(TK)
    ds, j = fc.surface_match(B, TP)
    on = {t: ds <= t for t in (4.0, 6.0, 8.0)}
    m = on[TOL]
    delta = FT[j[m]] - TT[m]
    rd = np.round(delta).astype(int)
    cstar = mode_int(rd); e = rd - cstar
    E = float((e != 0).mean()) if m.any() else None
    # per cell: 8 sectors (angle around the umbilicus, raw atan2) x label thirds
    from acc_lib import raw_angles
    th, _ = raw_angles(TP[m], umb); ang = np.degrees(th) % 360; km = TK[m]
    cells = {}
    for name, (lo, hi) in THIRDS.items():
        cells[name] = {}
        for sct in range(8):
            mm = (km >= lo) & (km <= hi) & (ang >= 45 * sct) & (ang < 45 * (sct + 1))
            if mm.sum() >= 500:
                cells[name][f"{45 * sct:03d}"] = dict(modal_e=mode_int(e[mm]), wrong_share=round(float((e[mm] != 0).mean()), 4),
                                                      n=int(mm.sum()))
    om = [v["modal_e"] for v in cells["outer"].values()]
    per_wrap = {int(k): round(float((e[km == k] != 0).mean()), 4) for k in np.unique(km)}
    modal_w = {int(k): mode_int(e[km == k]) for k in np.unique(km)}
    cm = ~np.isin(km, DROP_CLEAN)
    if cm.any():
        cc = mode_int(rd[cm]); ec = rd[cm] - cc
        clean = dict(E_clean=round(float((ec != 0).mean()), 5), c_star_clean=cc, n_matched_clean=int(cm.sum()),
                     M2_clean_6vx=round(float(cm.sum() / (~np.isin(TK, DROP_CLEAN)).sum()), 4),
                     breaks_clean=breaks({k: v - (cc - cstar) for k, v in modal_w.items() if k not in DROP_CLEAN}))
    else:
        clean = dict(E_clean=None)
    u, c = np.unique(e, return_counts=True)
    return dict(
        E_wrong_turn_share=round(E, 5) if E is not None else None,
        n_truth=int(len(TT)), n_matched=int(m.sum()),
        M2_on_sheet={f"{int(t)}vx": round(float(on[t].mean()), 4) for t in on},
        M1_right_turn_share=round(float(((e == 0).sum()) / len(TT)), 4),
        c_star=cstar, e_distribution={int(a): int(b) for a, b in zip(u, c)},
        frac_residual_median=round(float(np.median(np.abs(delta - np.round(delta)))), 4),
        theta_ref_deg=round(float(np.degrees(theta_ref)), 1), seam_gap_deg=round(float(np.degrees(gap)), 1),
        fit_windings=[int(min(fm)), int(max(fm))], n_fit_vertices=int(len(FT)),
        cells=cells, outer_sector_range=(max(om) - min(om)) if om else None, wrong_share_per_wrap=per_wrap,
        modal_e_per_wrap=modal_w, breaks=breaks(modal_w), **clean,
        point_to_surface_matched_median=round(float(np.median(ds[m])), 3) if m.any() else None)


if __name__ == "__main__":
    out = sys.argv[1]
    umb = load_umbilicus(UMB)
    truth = truth_meshes("9um")
    s, sv = spiral_sense(truth, umb)
    tseams = seam_angles(truth, umb, s, zr=BAND)
    res = json.load(open(out)) if os.path.exists(out) else {}
    res["_meta"] = dict(spiral_sense=s, sense_detail=sv, tol_vx=TOL, score_z=SCORE_Z, band=BAND, truth_wraps=sorted(truth),
                        instrument="lanes/A_geometry/tools/fit_compare.py surface_match")
    for arg in sys.argv[2:]:
        name, rd = arg.split("=", 1)
        r = score(rd, umb, s, truth, tseams); r["run_dir"] = os.path.relpath(rd, ROOT)
        res[name] = r
        print(f"{name}: E {r['E_wrong_turn_share']} | E_clean {r.get('E_clean')} | breaks {r['breaks']['n']} | M2(6vx) {r['M2_on_sheet']['6vx']} | M1 {r['M1_right_turn_share']} | c* {r['c_star']} "
              f"| frac {r['frac_residual_median']} | seam gap {r['seam_gap_deg']} | outer range {r['outer_sector_range']}", flush=True)
        json.dump(res, open(out, "w"), indent=1)
