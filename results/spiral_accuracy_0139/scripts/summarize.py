"""Apply the frozen decision rules (prereg section 6) to out/accuracy.json (+ out/h2_fit_compare.json if present).
Writes out/verdict.json and prints a markdown table. Run from the repo root:
  .venv/bin/python -I work/spiral_accuracy/scripts/summarize.py"""
import json, os, re, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
O = os.path.join(ROOT, "work/spiral_accuracy/out")
acc = json.load(open(os.path.join(O, "accuracy.json")))
fits = {k: v for k, v in acc.items() if not k.startswith("_") and not k.startswith("smoke")}
base = {k: v for k, v in fits.items() if k.startswith("base_")}
scorable = {k: v for k, v in fits.items() if v["M2_on_sheet"]["6vx"] >= 0.2 and v["E_wrong_turn_share"] is not None}
bs = {k: v for k, v in base.items() if k in scorable}
out = dict(rules="prereg/2026-10-13_partner_spiral_accuracy.md section 6", base_fits=sorted(bs), unscorable=sorted(set(fits) - set(scorable)))
void = len(bs) < 2 or np.median([v["M2_on_sheet"]["6vx"] for v in base.values()]) < 0.5
out["VOID"] = bool(void)


def rules(key):
    Eb = [v[key] for v in bs.values()]
    Emax, Emin = max(Eb), min(Eb)
    noisy = Emax - Emin > 0.10
    arms = {}
    for k, v in scorable.items():
        mm = re.match(r"([xy])(\d+)_s(\d+)", k)
        if mm and v.get(key) is not None:
            arms.setdefault((mm.group(1), int(mm.group(2))), {})[k] = v[key]
    res = {}
    for (ax, sh), d in sorted(arms.items()):
        E = list(d.values()); r = dict(fits=d, n=len(E), mm_shift=round(sh * 9.362 / 1000, 3),
                                       headline_pp=round(100 * (float(np.median(E)) - float(np.median(Eb))), 2))
        for m in (0.005, 0.01, 0.02):
            if len(E) < 2:
                c = "NOT RUN (arm incomplete: 1 of 2 fits)"
            elif noisy:
                c = "UNDECIDED"
            elif all(e < Emin - m for e in E):
                c = "GAIN"
            elif all(e > Emax + m for e in E):
                c = "ACCURACY COST"
            elif all(e <= Emax + m for e in E):
                c = "TOLERATED"
            else:
                c = "UNDECIDED"
            r[f"verdict_m{m}"] = c
        r["verdict"] = r["verdict_m0.01"]
        res[f"{ax}{sh}"] = r
    return dict(sorted((k, v[key]) for k, v in bs.items())), Emax, Emin, noisy, res


if not void:
    Eb_, Emax, Emin, noisy, res = rules("E_wrong_turn_share")
    out.update(E_base=Eb_, E_max=Emax, E_min=Emin, seed_noise_too_large=bool(noisy), arms=res)
    Ebc, Emaxc, Eminc, noisyc, resc = rules("E_clean")
    out.update(E_base_clean=Ebc, E_max_clean=Emaxc, E_min_clean=Eminc, seed_noise_too_large_clean=bool(noisyc), arms_clean=resc)
# H2 (descriptive): the package's fit-vs-fit rule on the same fits
h2p = os.path.join(O, "h2_fit_compare.json")
if os.path.exists(h2p):
    h2 = json.load(open(h2p))
    bp = {k: v for k, v in h2.items() if k.startswith("base_") and "|base_" in k}
    floor = min(v["turn_share_at_mode"] for v in bp.values()) if bp else None
    tab = {}
    for k, v in h2.items():
        if k == "grid_floor" or k in bp:
            continue
        other = k.split("|")[1]
        mm = re.match(r"([xy]\d+)_s\d+", other)
        if not mm:
            continue
        nt = (floor is not None and v["turn_share_at_mode"] < floor) or (v["outer_sector_range"] is not None and v["outer_sector_range"] > 1)
        tab.setdefault(mm.group(1), []).append(bool(nt))
    out["H2"] = dict(floor=floor, base_pairs={k: v["turn_share_at_mode"] for k, v in bp.items()},
                     not_tolerated_pairs={a: f"{sum(x)}/{len(x)}" for a, x in tab.items()})
json.dump(out, open(os.path.join(O, "verdict.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
