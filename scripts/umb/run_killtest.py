"""ChargeCore kill test, evaluated once on the frozen run (pre-registration: 04_Work_Log.md 20:02 + 20:06 CDT).
Steps: verify chargecore.py sha256; SI 19:40 version (all mod-50 samples) and C2 version (excluded-window samples
dropped), both smoothed the same way; ChargeCore on the C2 z list; C5 coverage/gap rule; C6 N+; C7 outcome table.
Usage: run_killtest.py [sigma]   (sigma only for the pre-registered MARGINAL LOSO second look)."""
import hashlib, json, os, sys, numpy as np
from concurrent.futures import ProcessPoolExecutor
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import score, gridio, chargecore as cc
SIGMA = float(sys.argv[1]) if len(sys.argv) > 1 else cc.SIGMA
TAG = "" if SIGMA == cc.SIGMA else f"_s{int(SIGMA)}"
assert hashlib.sha256(open(os.path.join(H, "chargecore.py"), "rb").read()).hexdigest() == open(os.path.join(H, "chargecore.sha256")).read().split()[0], "chargecore.py changed after freeze"
MM = score.MM


def in_excl(sc, z):
    return any(a <= z <= b for a, b in score.EXCL.get(sc, []))


def smooth(Z, XY, K=5):
    """running median in x and y separately over +-K positions of the z list; missing (NaN) skipped; output only
    where the position has its own estimate."""
    out = []
    for i in range(len(Z)):
        if np.isnan(XY[i, 0]):
            continue
        w = XY[max(0, i - K):i + K + 1]; w = w[~np.isnan(w[:, 0])]
        out.append((Z[i], np.median(w[:, 0]), np.median(w[:, 1])))
    return np.array(out)


def score_c5(sc, zxy):
    """19:40 scoring plus C5: a truth point with no estimate within +-500 z scores 10 mm."""
    P = score.scored_points(sc); zxy = zxy[np.argsort(zxy[:, 0])]
    ex = np.interp(P[:, 2], zxy[:, 0], zxy[:, 1]); ey = np.interp(P[:, 2], zxy[:, 0], zxy[:, 2])
    e = np.hypot(ex - P[:, 0], ey - P[:, 1]) * MM
    near = np.array([np.min(np.abs(zxy[:, 0] - z)) <= 500 for z in P[:, 2]])
    e[~near] = 10.0
    return float(np.median(e)), int((~near).sum()), e


def run_slice(args):
    sc, z = args
    f = os.path.join(H, "ng", sc, "xy", f"{z:06d}.grid")
    b = gridio.bounds(f); o = cc.detect(gridio.read_grid(f), b["w"], b["h"], sigma=SIGMA)
    return z, o["est"], o["n_plus"], o["n_cands"]


def main():
    SI, SIc, CCs, rep = {}, {}, {}, {}
    env = json.load(open(os.path.join(H, "env_scores.json")))
    for sc in score.SCROLLS:
        cp = json.load(open(os.path.join(H, f"strong_{sc}.json")))["control_points"]
        A = np.array(sorted((p["z"], p["x"], p["y"]) for p in cp), float)
        SI[sc] = score.score(sc, smooth(A[:, 0], A[:, 1:]))["median_mm"]                     # 19:40 version
        Ac = A[[not in_excl(sc, z) for z in A[:, 0]]]
        SIc[sc] = score.score(sc, smooth(Ac[:, 0], Ac[:, 1:]))["median_mm"]                  # C2 version (gated)
        zl = [int(z) for z in Ac[:, 0]]
        with ProcessPoolExecutor(12) as ex:
            res = list(ex.map(run_slice, [(sc, z) for z in zl]))
        XY = np.array([[e[0], e[1]] if e else [np.nan, np.nan] for _, e, _, _ in res], float)
        raw = np.array([(z, *xy) for z, xy in zip(zl, XY)], float)
        np.save(os.path.join(H, f"chargecore_raw_{sc}{TAG}.npy"), raw)
        S = smooth(np.array(zl, float), XY); np.save(os.path.join(H, f"chargecore_{sc}{TAG}.npy"), S)
        cov = float(np.mean(~np.isnan(XY[:, 0]))); nplus = float(np.median([r[2] for r in res]))
        D, ngap, e = score_c5(sc, S) if len(S) else (10.0, -1, None)
        B = min(SIc[sc], env[sc]["median_mm"])
        ks = [r[1][2] for r in res if r[1]]
        CCs[sc] = D
        rep[sc] = dict(D=D, SI_c2=SIc[sc], SI_1940=SI[sc], ENV=env[sc]["median_mm"], B=B, beats=bool(D <= 0.9 * B),
                       ratio_D_SI=D / SIc[sc], coverage=cov, N_plus=nplus, n_slices=len(zl), n_truth_gap10mm=ngap,
                       share_le_094=float(np.mean(e <= 0.94)) if e is not None else 0.0,
                       D_smooth_truth=score.score(sc, S)["median_mm_smooth_truth"] if len(S) else None,
                       n_core=ks.count("core"), n_pair=ks.count("pair"))
        print(sc, json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in rep[sc].items()}), flush=True)
    s = score.SCROLLS
    gate = sum(rep[x]["D"] <= 0.94 for x in s) >= 3 and sum(rep[x]["D"] < rep[x]["SI_c2"] for x in s) >= 3
    kill_legs = {
        "D>B on >=2": sum(rep[x]["D"] > rep[x]["B"] for x in s) >= 2,
        "D>1.5xSI on 0125/0211/0139": any(rep[x]["D"] > 1.5 * rep[x]["SI_c2"] for x in ["PHerc0125", "PHerc0211", "PHerc0139"]),
        "coverage<70% on >=2": sum(rep[x]["coverage"] < 0.7 for x in s) >= 2,
        "N+>6 on >=2": sum(rep[x]["N_plus"] > 6 for x in s) >= 2}
    cont = (sum(rep[x]["beats"] for x in s) >= 3 and all(rep[x]["coverage"] >= 0.7 for x in s)
            and sum(rep[x]["N_plus"] <= 6 for x in s) >= 3)
    outcome = "KILL" if any(kill_legs.values()) else ("GATE PASS" if gate else ("CONTINUE" if cont else "MARGINAL"))
    if gate and any(kill_legs.values()):
        outcome = "GATE PASS (flagged: KILL leg fired; to the user before any build)"
    ratios = {x: rep[x]["ratio_D_SI"] for x in s}
    summary = dict(sigma=SIGMA, outcome=outcome, gate_1940=bool(gate), kill_legs=kill_legs, continue_rule=bool(cont),
                   mechanistic_prediction_0826_lowest_ratio=bool(min(ratios, key=ratios.get) == "PHerc0826"), per_scroll=rep)
    json.dump(summary, open(os.path.join(H, f"killtest_result{TAG}.json"), "w"), indent=1)
    print("OUTCOME:", outcome, "| kill legs:", kill_legs, "| 0826 lowest D/SI:", summary["mechanistic_prediction_0826_lowest_ratio"])


if __name__ == "__main__":
    main()
