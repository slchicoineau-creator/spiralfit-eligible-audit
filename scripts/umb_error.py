"""Errors of the two umbilicus estimators against the published umbilici (xy distance to the published polyline,
linearly interpolated in z), in voxels and mm. Inputs: results/umb_est_*.json (lasagna-normal convergence,
scripts/umb_estimate.py) and results/visual_umbilicus_readings.json (visual core marking). Published files are fetched
from the open-data bucket (or read from TRUTH_DIR if given). Writes results/umbilicus_estimate_errors.json.
Usage: python umb_error.py [TRUTH_DIR]"""
import glob, json, os, re, sys, urllib.request, numpy as np
R = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
S3 = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"


def published(sc, vid):
    if len(sys.argv) > 1:
        return json.load(open(glob.glob(os.path.join(sys.argv[1], f"{sc}_{vid}-umbilicus-*.json"))[-1]))["control_points"]
    keys = re.findall(r"<Key>([^<]*)</Key>", urllib.request.urlopen(f"{S3}/?list-type=2&prefix={sc}/representations/umbilicus/", timeout=60).read().decode())
    return json.loads(urllib.request.urlopen(f"{S3}/{sorted(k for k in keys if vid in k)[-1]}", timeout=60).read())["control_points"]


def errors(sc, vid, pts):
    cp = sorted(published(sc, vid), key=lambda c: c["z"]); z = [c["z"] for c in cp]
    return [float(np.hypot(x - np.interp(q, z, [c["x"] for c in cp]), y - np.interp(q, z, [c["y"] for c in cp]))) for q, x, y in pts]


out = {"lasagna_normal_convergence": {}, "visual_core_marking": {}}
for f in sorted(glob.glob(os.path.join(R, "umb_est_*.json"))):
    d = json.load(open(f)); sc = os.path.basename(f)[8:-5]; vid = d["metadata"]["volume"]
    e = errors(sc, vid, [(c["z"], c["x"], c["y"]) for c in d["control_points"]])
    out["lasagna_normal_convergence"][sc] = dict(n=len(e), median_vx=round(float(np.median(e))), median_mm=round(float(np.median(e)) * 9.362e-3, 2))
V = json.load(open(os.path.join(R, "visual_umbilicus_readings.json"))); pooled = []
for key, v in V.items():
    if key.startswith("_") or not v["points"] or "truth" in v:
        continue
    sc, vid = key.split("_"); e = [x * v["voxel_um"] * 1e-3 for x in errors(sc, vid, v["points"])]; pooled += e
    out["visual_core_marking"][key] = dict(n=len(e), median_mm=round(float(np.median(e)), 2), per_point_mm=[round(x, 2) for x in e])
out["visual_core_marking"]["pooled_2.4um"] = dict(n=len(pooled), median_mm=round(float(np.median(pooled)), 2))
out["visual_core_marking"]["PHerc0139_20250728140407 (9.362 um)"] = "not recomputed here; logged median 2.34 mm (see readings file)"
json.dump(out, open(os.path.join(R, "umbilicus_estimate_errors.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
