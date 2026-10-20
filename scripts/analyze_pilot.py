"""TrackFix pilot analysis, exactly as pre-registered (work log, pilot pre-registration): per scroll, mean blind rating
per arm (a repeated picture's two ratings averaged); CONTINUE to the full run iff mean(C)-mean(A) >= +0.25 on >= 1
scroll AND >= -0.25 on the other; otherwise STOP. B reported the same way (if B >= C on both scrolls the full run keeps
B as the key comparator). Usage: python analyze_pilot.py trackfix_pilot_ratings.json"""
import json, os, sys, numpy as np
HERE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
K = json.load(open(os.path.join(HERE, "pilot_rating_key.json"))); R = json.load(open(sys.argv[1]))
ans = R["answers"]
assert R.get("complete") and [a["id"] for a in ans] == K["order"], "ratings incomplete or do not match the key order"
by = {}
for a in ans:
    by.setdefault(a["id"], []).append(a["rating"])
reps = [v for v in by.values() if len(v) == 2]
print(f"self-consistency: {len(reps)} repeats, exact {np.mean([a == b for a, b in reps]):.2f}, within 1 {np.mean([abs(a - b) <= 1 for a, b in reps]):.2f}")
m = {}
for iid, v in K["key"].items():
    m.setdefault((v["scroll"], v["arm"]), []).append(float(np.mean(by[iid])))
res = {}
for sc in sorted({s for s, _ in m}):
    A, B, C = (float(np.mean(m[(sc, x)])) for x in "ABC")
    res[sc] = dict(A=A, B=B, C=C, C_minus_A=C - A, B_minus_A=B - A, C_minus_B=C - B, n={x: len(m[(sc, x)]) for x in "ABC"})
    print(f"{sc}: A {A:.2f}  B {B:.2f}  C {C:.2f}   C-A {C - A:+.2f}  B-A {B - A:+.2f}  C-B {C - B:+.2f}")
d = [r["C_minus_A"] for r in res.values()]
cont = any(x >= 0.25 for x in d) and all(x >= -0.25 for x in d)
bge = all(r["B"] >= r["C"] for r in res.values())
print("PILOT DECISION:", "CONTINUE to full run" if cont else "STOP", "| B >= C on both scrolls:" , bge)
json.dump(dict(per_scroll=res, continue_=cont, B_ge_C_both=bge), open(os.path.join(HERE, "pilot_results.json"), "w"), indent=1)
