"""EXPLORATORY, post-kill (not gated, not pre-registered): (1) x/y-swap sanity check; (2) is the charge signal present
but badly selected? distance from the truth to the nearest surviving +1/2 defect and to the nearest pair midpoint,
on up to 40 evenly spaced scored slices per scroll."""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score, truth, gridio, chargecore as cc
for sc in score.SCROLLS:
    R = np.load(f"chargecore_raw_{sc}.npy"); R = R[~np.isnan(R[:, 1])]
    T = truth.at(sc, R[:, 0]); ok = ~np.isnan(T[:, 0])
    e = np.median(np.hypot(*(R[ok, 1:] - T[ok]).T)); es = np.median(np.hypot(*(R[ok, 1:][:, ::-1] - T[ok]).T))
    zs = R[ok, 0][np.linspace(0, ok.sum() - 1, 40).astype(int)]; dp, dm, n = [], [], []
    for z in zs:
        f = f"ng/{sc}/xy/{int(z):06d}.grid"; b = gridio.bounds(f); o = cc.detect(gridio.read_grid(f), b["w"], b["h"])
        t = truth.at(sc, [z])[0]; P = np.array(o["plus_surv"]).reshape(-1, 2)
        if len(P):
            dp.append(np.min(np.hypot(*(P - t).T)))
            M = np.array([(P[i] + P[j]) / 2 for i in range(len(P)) for j in range(i + 1, len(P))]).reshape(-1, 2)
            if len(M): dm.append(np.min(np.hypot(*(M - t).T)))
        n.append(len(P))
    mm = score.MM
    print(f"{sc}: raw per-slice median err {e*mm:.2f} mm (x/y swapped: {es*mm:.2f}); nearest +1/2 to truth median "
          f"{np.median(dp)*mm:.2f} mm; best possible pair midpoint median {np.median(dm)*mm:.2f} mm; +1/2 per slice median {np.median(n):.0f}")
