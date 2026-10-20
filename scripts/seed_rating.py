"""TrackFix pilot blind rating page (pre-registered). Collects renders/<SCROLL>_<ARM>/crop*.png (12 per arm per scroll),
builds pilot_rating.html from the GrowGate rating template (opaque ids, random order, ~10% repeats, no two consecutive
from the same scroll+arm, repeats >= 20 apart) and pilot_rating_key.json (kept apart for blinding).
Usage: python pilot_rating.py [--allow-partial]"""
import base64, glob, hashlib, io, json, os, sys, numpy as np
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__)); GG = HERE
RNG = np.random.default_rng(2029); MIN_GAP = 20
items, key = [], {}
for d in sorted(glob.glob(os.path.join(HERE, "renders_seed", "*_*_s*"))):
    sc, arm, seed = os.path.basename(d).rsplit("_", 2)
    for f in sorted(glob.glob(os.path.join(d, "crop*.png"))):
        iid = hashlib.sha256(f"{sc}|{arm}|{seed}|{os.path.basename(f)}|trackfix-seedrep".encode()).hexdigest()[:10]
        b = io.BytesIO(); Image.open(f).save(b, "JPEG", quality=90)
        key[iid] = dict(scroll=sc, arm=arm, seed=seed, crop=os.path.basename(f))
        items.append(dict(id=iid, unit=f"{sc}|{arm}|{seed}", src="data:image/jpeg;base64," + base64.b64encode(b.getvalue()).decode()))
arms = {(v["scroll"], v["arm"], v["seed"]) for v in key.values()}
if len(arms) < 12 and "--allow-partial" not in sys.argv:
    sys.exit(f"only {len(arms)} scroll-arm sets rendered; need 12")
dups = [items[i] for i in RNG.choice(len(items), int(round(0.1 * len(items))), replace=False)]
def build_order():
    """Greedy randomised ordering: no two consecutive items from the same scroll+arm (pick among allowed units with
    the most items left, ties random), then insert each repeat >= MIN_GAP after its first showing, not adjacent to its
    own unit."""
    pools = {}
    for x in items:
        pools.setdefault(x["unit"], []).append(x)
    for u in pools:
        RNG.shuffle(pools[u])
    seq, prev = [], None
    while any(pools.values()):
        cand = [u for u, v in pools.items() if v and u != prev] or [u for u, v in pools.items() if v]
        m = max(len(pools[u]) for u in cand); top = [u for u in cand if len(pools[u]) == m]
        u = top[int(RNG.integers(len(top)))]; seq.append(pools[u].pop()); prev = u
    for dd in dups:
        first = next(i for i, x in enumerate(seq) if x["id"] == dd["id"])
        slots = [p for p in range(first + MIN_GAP, len(seq) + 1) if (p == 0 or seq[p - 1]["unit"] != dd["unit"]) and (p == len(seq) or seq[p]["unit"] != dd["unit"])]
        if not slots:
            return None
        seq.insert(int(RNG.choice(slots)), dd)
    return seq


seq = None
while seq is None:
    seq = build_order()
assert not any(a["unit"] == b["unit"] for a, b in zip(seq, seq[1:])), "adjacent same unit"
small = json.loads(open(os.path.join(GG, "reference_examples.html"), encoding="utf-8").read().split("const D = ")[1].split(";\nconst SCALE")[0])["lesson"]
html = open(os.path.join(GG, "rating_template.html"), encoding="utf-8").read()
html = html.replace("/*DATA*/null", json.dumps(dict(items=[dict(id=x["id"], src=x["src"]) for x in seq], ref={g: [dict(src=e["src"]) for e in v] for g, v in small.items()})))
html = html.replace("Sheet or Swirl Rating", "Sheet or Swirl Rating (seed check)").replace("sheet_or_swirl_ratings.json", "trackfix_seed_ratings.json")
tag = "_TEST_partial" if "--allow-partial" in sys.argv else ""
open(os.path.join(HERE, f"seed_rating{tag}.html"), "w", encoding="utf-8").write(html)
json.dump(dict(key=key, order=[x["id"] for x in seq], repeats=[d["id"] for d in dups]), open(os.path.join(HERE, f"seed_rating_key{tag}.json"), "w"), indent=1)
print(f"{len(items)} crops + {len(dups)} repeats = {len(seq)} pictures from {len(arms)} scroll-arm sets")
