"""Umbilicus tolerance curve (pre-registered 2026-09-29 19:36 CDT, 04_Work_Log.md): PHerc0826 [8136, 9136), arm A.
Each fit (published umbilicus shifted by a fixed vector, seed 1; or unshifted with another seed) is compared with the
baseline fit (published umbilicus, seed 1) by fit_compare.py.
Rule: a shift is TOLERATED iff outer sector range <= 1 AND share at the modal turn offset >= (lowest seed pair) - 0.10.
Tolerance per direction = the largest shift on the ladder tolerated with every smaller shift in that direction also
tolerated. Headline = min(+x, +y). If a seed pair itself has outer sector range > 1, the range limit becomes the largest
seed range. Writes results/tolerance_curve.json.
Usage: python tol_analyze.py FITS_ROOT UMBILICUS.json   (FITS_ROOT holds full/, sens/, tol/ as written by the scripts)"""
import glob, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fit_compare import compare
F, U = sys.argv[1], sys.argv[2]
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "tolerance_curve.json")
d = lambda p: (glob.glob(os.path.join(F, p, "*/")) or [None])[0]
cp = sorted(json.load(open(U))["control_points"], key=lambda c: c["z"])
umb = tuple(np.array([c[k] for c in cp], float) for k in ("z", "y", "x"))
base = d("full/PHerc0826_8136_A")
fits = {"seed2": "sens/PHerc0826_8136_R", "seed3": "tol/PHerc0826_8136_seed3", "seed4": "tol/PHerc0826_8136_seed4",
        "x50": "tol/PHerc0826_8136_x50", "x100": "tol/PHerc0826_8136_x100", "x250": "sens/PHerc0826_8136_S",
        "x500": "tol/PHerc0826_8136_x500", "y50": "tol/PHerc0826_8136_y50", "y100": "tol/PHerc0826_8136_y100",
        "y250": "tol/PHerc0826_8136_y250", "y500": "tol/PHerc0826_8136_y500", "xm250": "tol/PHerc0826_8136_xm250"}
res = json.load(open(OUT))["fits"] if os.path.exists(OUT) else {}
for name, p in fits.items():
    run = d(p)
    if run is None or not glob.glob(os.path.join(run, "meshes", "fitted", "w*")):
        print(f"{name}: missing (fit failed or not run)"); continue
    if name not in res:
        res[name] = compare(base, run, umb)
        sat = json.load(open(os.path.join(run, "satisfaction_metrics_fitted.json")))["summary"]
        res[name]["satisfied_tracks_fraction"] = round(sat["satisfied_tracks_fraction"], 4)
        res[name]["satisfied_track_points_fraction"] = round(sat["satisfied_track_points_fraction"], 4)
    r = res[name]
    print(f"{name}: p2s median {r['point_to_surface']['median']} within4 {r['point_to_surface']['within4']} | "
          f"turn {r['turn_modal_offset']:+d} on {r['turn_share_at_mode']} | outer sectors {r['sector_modal_offset']['outer']} "
          f"range {r['outer_sector_range']}", flush=True)
    json.dump({"fits": res}, open(OUT, "w"), indent=1)

seeds = [res[s] for s in ("seed2", "seed3", "seed4") if s in res]
share_floor = min(s["turn_share_at_mode"] for s in seeds) - 0.10
range_limit = max(1, max(s["outer_sector_range"] for s in seeds))
tol = {}
for direction, ladder in (("+x", [50, 100, 250, 500]), ("+y", [50, 100, 250, 500])):
    t = 0; verdicts = {}
    for mag in ladder:
        key = f"{direction[1]}{mag}"
        if key not in res:
            verdicts[mag] = "missing"; break
        ok = res[key]["outer_sector_range"] <= range_limit and res[key]["turn_share_at_mode"] >= share_floor
        verdicts[mag] = "tolerated" if ok else "NOT tolerated"
        if not ok:
            break
        t = mag
    complete = "missing" not in verdicts.values()
    tol[direction] = dict(tolerance_vx=t if complete else None, tolerance_mm=round(t * 9.362e-3, 2) if complete else None,
                          all_tolerated=complete and t == ladder[-1], verdicts=verdicts)
if any(v["tolerance_vx"] is None for v in tol.values()):
    headline = "incomplete (a required fit is missing)"
else:
    head = min(v["tolerance_vx"] for v in tol.values())
    headline = ("< 50 vx (0.47 mm)" if head == 0 else
                f">= 500 vx (4.68 mm): every shift tolerated" if all(v["all_tolerated"] for v in tol.values()) else
                f"{head} vx ({head * 9.362e-3:.2f} mm)")
summary = dict(share_floor=round(share_floor, 3), range_limit=range_limit, per_direction=tol, headline=headline)
json.dump({"rule": __doc__.split("\n")[2:6], "summary": summary, "fits": res}, open(OUT, "w"), indent=1)
print(json.dumps(summary, indent=1))
