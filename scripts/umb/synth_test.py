"""Code-correctness test for chargecore.py (clarification 7; not evidence for the method).
(a) round Archimedean spiral -> one +1 core at the centre; (b) the same squeezed 5x in y -> +1/2 pair at the slit tips,
estimate at the centre; (c) parallel layers with a nested-ellipse buckle eye -> net charge 0 over the eye."""
import numpy as np
import chargecore as cc
rng = np.random.default_rng(0); C = 1500.0


def pieces(pts):
    """cut a dense curve into normal-grid-like polylines: 20 vx steps, 3 points per piece, 10% of pieces dropped."""
    s = np.r_[0, np.cumsum(np.hypot(*np.diff(pts, axis=0).T))]; t = np.arange(0, s[-1], 20.0)
    q = np.c_[np.interp(t, s, pts[:, 0]), np.interp(t, s, pts[:, 1])]
    out = [np.rint(q[i:i + 3]).astype(np.int32) for i in range(0, len(q) - 2, 3) if rng.random() > 0.1]
    return out


def spiral(r0, pitch, rmax, squeeze=1.0):
    th = np.linspace(0, (rmax - r0) / pitch * 2 * np.pi, 200000); r = r0 + pitch * th / (2 * np.pi)
    return pieces(np.c_[C + r * np.cos(th), C + squeeze * r * np.sin(th)])


def report(name, paths, expect):
    o = cc.detect(paths, 3000, 3000)
    e = o["est"]; err = None if e is None else float(np.hypot(e[0] - C, e[1] - C))
    print(f"{name}: est={e} err={err if err is None else round(err, 1)} vx  n_plus={o['n_plus']} n_cands={o['n_cands']} "
          f"cores={o['n_cores']} surv+={np.round(o['plus_surv'], 0).tolist()[:6]}  expect: {expect}")
    return o


report("(a) round", spiral(20, 20, 1300), "one core within ~20 vx of centre")
report("(b) squeezed 5x", spiral(300, 70, 1400, 0.2), "+1/2 pair near x=1500+-300, estimate near centre")
# (c) buckle eye: nested ellipses in a band of parallel layers; outer layers bulge around it
paths = []
for j in range(6):
    a, b = 420 - 45 * j, 84 - 14 * j
    t = np.linspace(0, 2 * np.pi, 4000); paths += pieces(np.c_[C + a * np.cos(t), C + b * np.sin(t)])
x = np.linspace(100, 2900, 6000)
for k in list(range(-60, -6)) + list(range(7, 61)):
    bump = 14 * np.clip(1 - ((x - C) / 520) ** 2, 0, None) * np.sign(k) * (7 / abs(k)) ** 0.5 * 1.0
    paths += pieces(np.c_[x, C + 14 * k + bump])
o = cc.detect(paths, 3000, 3000, return_field=True)
th, ok = o["theta"], o["valid"] & o["body"]
a_, b_, c_, d_ = th[:-1, :-1], th[:-1, 1:], th[1:, 1:], th[1:, :-1]; w = lambda v: (v + np.pi) % (2 * np.pi) - np.pi
q = (w(b_ - a_) + w(c_ - b_) + w(d_ - c_) + w(a_ - d_)) / (4 * np.pi)
okp = ok[:-1, :-1] & ok[:-1, 1:] & ok[1:, 1:] & ok[1:, :-1]
i0, i1, j0, j1 = int((C - 300) / 8), int((C + 300) / 8), int((C - 800) / 8), int((C + 800) / 8)
print(f"(c) buckle eye: net charge in box +-800 x +-300 = {np.rint(2 * q[i0:i1, j0:j1] * okp[i0:i1, j0:j1]).sum() / 2}  "
      f"(expect 0); est={o['est']} n_plus={o['n_plus']} n_cands={o['n_cands']} (spec property, reported only)")
