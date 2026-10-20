"""SPIRAL-ACCURACY-0139, judge check (read-only, CPU): per-truth-vertex concordance of the turn error e between fits,
inside the labelled footprint. Backs judge finding 7 of reviews/2026-10-14_partner_judge_spiral_accuracy.md (ba22d215):
"most of the extra renumbering goes from one wrong number to another, which E does not count".

Re-implements score_accuracy.score() (same acc_lib calls, same fit_compare.surface_match, same 6 vx tolerance, own c*) but
keeps e per truth vertex, then compares e between two fits on the truth vertices matched by both. Committed copy of the
judge's run script (scratchpad judge_concord.py, 2026-10-14 06:59), with repo-relative paths and a summary block added.

Needs the local fits (work/spiral_accuracy/fits/<fit>/meshes/fitted) and data (work/spiral_accuracy/data); both git-ignored.
Usage (repo root): .venv/bin/python -I work/spiral_accuracy/judge_concordance/concordance.py
Writes work/spiral_accuracy/judge_concordance/concordance.json."""
import importlib.util, itertools, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(ROOT, "work/spiral_accuracy/scripts"))
from acc_lib import (D, BAND, SCORE_Z, load_umbilicus, truth_meshes, fit_meshes, spiral_sense, seam_angles,  # noqa: E402
                     theta_ref_from, mesh_tau)
spec = importlib.util.spec_from_file_location("fc", os.path.join(ROOT, "lanes/A_geometry/tools/fit_compare.py"))
fc = importlib.util.module_from_spec(spec); spec.loader.exec_module(fc)

FITS = ["base_s1", "base_s2", "base_s3", "y10_s1", "y10_s2", "y25_s1", "y25_s2", "y50_s1", "y50_s2", "x25_s1", "x50_s1"]
UNSET = 10**6

umb = load_umbilicus(os.path.join(D, "PHerc0139_20250728140407", "umbilicus.json"))
truth = truth_meshes("9um"); s, _ = spiral_sense(truth, umb); tseams = seam_angles(truth, umb, s, zr=BAND)
acc = json.load(open(os.path.join(ROOT, "work/spiral_accuracy/out/accuracy.json")))

res, per_fit = {}, {}
for f in FITS:
    fm = fit_meshes(os.path.join(ROOT, "work/spiral_accuracy/fits", f))
    th, _ = theta_ref_from(np.concatenate([tseams, seam_angles(fm, umb, s, zr=BAND)]))
    FP, FT, FN = [], [], []
    for W, d in fm.items():
        P, T, N = mesh_tau(d, W, umb, s, th, zr=(BAND[0] - 40, BAND[1] + 40)); FP.append(P); FT.append(T); FN.append(N)
    B = dict(P=np.concatenate(FP), N=np.concatenate(FN)); FT = np.concatenate(FT)
    TP, TT, TK = [], [], []
    for k, d in truth.items():
        P, T, _ = mesh_tau(d, k, umb, s, th, zr=SCORE_Z); TP.append(P); TT.append(T); TK.append(np.full(len(T), k))
    TP, TT, TK = np.concatenate(TP), np.concatenate(TT), np.concatenate(TK)
    ds, j = fc.surface_match(B, TP)
    m = ds <= 6.0
    rd = np.full(len(TT), UNSET, int); rd[m] = np.round(FT[j[m]] - TT[m]).astype(int)
    u, c = np.unique(rd[m], return_counts=True); cstar = int(u[np.argmax(c)])
    e = np.where(m, rd - cstar, UNSET)
    E = round(float((e[m] != 0).mean()), 5)
    E_ref = acc[f]["E_wrong_turn_share"]
    per_fit[f] = dict(E=E, E_accuracy_json=E_ref, matches=bool(abs(E - E_ref) < 5e-6), n_matched=int(m.sum()), c_star=cstar)
    res[f] = dict(e=e, m=m, TK=TK)
    print(f, per_fit[f], flush=True)

pairs = {}
for a, b in itertools.combinations(FITS, 2):
    A, Bb = res[a], res[b]
    both = A["m"] & Bb["m"]
    ea, eb, tk = A["e"][both], Bb["e"][both], A["TK"][both]
    diff = ea != eb
    modal = sum(1 for k in np.unique(tk)
                if np.bincount(ea[tk == k] + 50).argmax() != np.bincount(eb[tk == k] + 50).argmax())
    pairs[f"{a}|{b}"] = dict(n_both=int(both.sum()), share_e_differs=round(float(diff.mean()), 4),
                             right_to_wrong=round(float(((ea == 0) & (eb != 0)).mean()), 4),
                             wrong_to_right=round(float(((ea != 0) & (eb == 0)).mean()), 4),
                             wrong_to_other_wrong=round(float(((ea != 0) & (eb != 0) & diff).mean()), 4),
                             wraps_with_modal_change=int(modal))
    print(a, b, pairs[f"{a}|{b}"], flush=True)

# Summary by pair type, as in the review's finding-7 table: base vs base, and each shifted arm vs the 3 bases.
# right->wrong / wrong->right are pooled over both orientations (a pair "base|arm" counts base->arm).
BASES = ["base_s1", "base_s2", "base_s3"]
groups = {"base vs base": [f"{a}|{b}" for a, b in itertools.combinations(BASES, 2)]}
for arm in ("y10", "y25", "y50", "x25", "x50"):
    groups[f"{arm} vs bases"] = [f"{b}|{f}" for b in BASES for f in FITS if f.startswith(arm + "_")]
summary = {}
for g, keys in groups.items():
    rows = [pairs[k] for k in keys]
    # 2 decimals in percent = exact from the 4-decimal pair values (the review's table rounds these to 1 decimal).
    rng = lambda col: [round(min(r[col] for r in rows) * 100, 2), round(max(r[col] for r in rows) * 100, 2)]
    summary[g] = dict(n_pairs=len(rows), e_differs_pct=rng("share_e_differs"), right_to_wrong_pct=rng("right_to_wrong"),
                      wrong_to_right_pct=rng("wrong_to_right"), wrong_to_other_wrong_pct=rng("wrong_to_other_wrong"),
                      wraps_modal_change=[min(r["wraps_with_modal_change"] for r in rows),
                                          max(r["wraps_with_modal_change"] for r in rows)])
    print(g, summary[g], flush=True)

json.dump(dict(_meta=dict(what="per-truth-vertex concordance of turn error e between fits (judge finding 7)",
                          tol_vx=6.0, score_z=list(SCORE_Z), band=list(BAND), fits=FITS,
                          instrument="lanes/A_geometry/tools/fit_compare.py surface_match (= package scripts/fit_compare.py)",
                          all_E_match_accuracy_json=all(v["matches"] for v in per_fit.values())),
               per_fit=per_fit, summary_by_pair_type=summary, pairs=pairs),
          open(os.path.join(HERE, "concordance.json"), "w"), indent=1)
print("wrote", os.path.join(HERE, "concordance.json"))
