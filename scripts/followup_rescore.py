"""Re-score the follow-up tolerance tests (README §4.3) from their shipped comparison files.

Each results/followups/<test>/compares.json holds fit_compare.py outputs for pairs "A|B" (A = an unshifted baseline fit,
B = a shifted fit, or two fits of the same shift). The frozen rule (README §4.2): a pair is NOT tolerated iff
  turn_share_at_mode < floor   (share leg)   or   outer sector range > 1   (range leg),
with floor = lowest unshifted seed-pair share - 0.10 (0.538 on PHerc0826 [8136, 9136), 0.791 on PHerc0211 [10496, 11496)).

Prints, per test and arm: cross pairs not tolerated (and by which leg), within-shift pairs not tolerated, the count at
share margins 0.04 / 0.10 / 0.15, and the floor-free ordering (max cross share vs min within-shift share).
Usage: python scripts/followup_rescore.py [results/followups]"""
import json, os, sys

B26, B11 = ["seed1", "seed2", "seed3", "seed4"], ["s1", "s2", "s3"]
s = lambda p, n: [f"{p}_s{i}" for i in range(1, n + 1)]
TESTS = {   # test: (floor, baselines, {arm: shifted fits})
    "seed5": (0.538, B26, {"+x 0.47 mm": ["x50"] + s("x50", 5)[1:], "+y 0.47 mm": ["y50"] + s("y50", 5)[1:]}),
    "tol25": (0.538, B26, {"+y 0.23 mm": ["y25", "y25_s2", "y25_s3"]}),
    "drift": (0.538, B26, {"tilt": s("tilt82", 3)}),
    "baseswap": (0.538, ["seed2", "seed3", "seed4"], {"+x 0.47 mm": ["x50"], "+y 0.47 mm": ["y50"]}),
    "tol0211": (0.791, B11, {"+y 0.47 mm": s("y50", 3), "+x 0.47 mm": s("x50", 3)}),
    "y25_0211": (0.791, B11, {"+y 0.23 mm": s("y25", 3)}),
    "y10_0211": (0.791, B11, {"+y 0.09 mm": s("y10", 3)}),
    "realaxis_0211": (0.791, B11, {"real axis": s("realax", 3)}),
}


def rng(c):
    v = list(c["sector_modal_offset"]["outer"].values())
    return max(v) - min(v)


def main(root):
    for t, (floor, bases, arms) in TESTS.items():
        c = json.load(open(os.path.join(root, t, "compares.json"), encoding="utf-8"))
        get = lambda a, b: c.get(f"{a}|{b}") or c.get(f"{b}|{a}")
        for arm, fits in arms.items():
            cross = [get(b, f) for b in bases for f in fits if get(b, f)]
            within = [get(a, b) for i, a in enumerate(fits) for b in fits[i + 1:] if get(a, b)]
            bad = lambda ps, fl: sum(p["turn_share_at_mode"] < fl or rng(p) > 1 for p in ps)
            share_leg = sum(p["turn_share_at_mode"] < floor for p in cross)
            range_leg = sum(rng(p) > 1 for p in cross)
            sweep = {m: bad(cross, floor + 0.10 - m) for m in (0.04, 0.10, 0.15)}
            mx = max(p["turn_share_at_mode"] for p in cross)
            mn = min((p["turn_share_at_mode"] for p in within), default=None)
            print(f"{t:14} {arm:12} cross not tolerated {bad(cross, floor)}/{len(cross)} "
                  f"(share leg {share_leg}, range leg {range_leg}); within not tolerated {bad(within, floor)}/{len(within)}; "
                  f"margins 0.04/0.10/0.15 -> {sweep[0.04]}/{sweep[0.10]}/{sweep[0.15]}; "
                  f"max cross {mx:.3f}" + (f" vs min within {mn:.3f}" if mn is not None else ""))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "followups"))
