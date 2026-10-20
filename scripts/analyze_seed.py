"""Seed-replicated pilot analysis (pre-registered): per scroll, effect = mean over seeds of C minus mean over seeds of A
(per-seed means; repeats averaged); noise = range of A's 3 seed means; 'real' iff effect > 0 and effect > noise."""
import json, os, sys, numpy as np
HERE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
K = json.load(open(os.path.join(HERE, "seed_rating_key.json"))); R = json.load(open(sys.argv[1]))
assert R.get("complete") and [a["id"] for a in R["answers"]] == K["order"], "incomplete or mismatched"
by = {}
for a in R["answers"]:
    by.setdefault(a["id"], []).append(a["rating"])
m = {}
for iid, v in K["key"].items():
    m.setdefault((v["scroll"], v["arm"], v["seed"]), []).append(float(np.mean(by[iid])))
real = 0
for sc in sorted({k[0] for k in m}):
    A = [np.mean(m[(sc, "A", s)]) for s in ("s1", "s2", "s3")]; C = [np.mean(m[(sc, "C", s)]) for s in ("s1", "s2", "s3")]
    eff, noise = np.mean(C) - np.mean(A), max(A) - min(A); ok = eff > 0 and eff > noise; real += ok
    print(f"{sc}: A seeds {np.round(A, 2)} C seeds {np.round(C, 2)} effect {eff:+.2f} A-noise range {noise:.2f} -> {'REAL' if ok else 'not shown'}")
print("VERDICT:", {2: "supported (2 scrolls)", 1: "mixed", 0: "not supported"}[real])
