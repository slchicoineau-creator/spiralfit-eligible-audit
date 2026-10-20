"""Diagnostic: WHERE is the giant (size-capped) component on each census ribbon? For each scroll's middle census ribbon,
report the share of giant-component pixels that lie OUTSIDE the scroll body (the filled hull of the kept components,
dilated 3 px) and save a small image (grey = kept, red = giant). Writes results/giant_where.json."""
import json, os, re, sys, urllib.request, numpy as np, cc3d
from scipy.ndimage import binary_fill_holes, binary_dilation, binary_closing
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib")); os.environ.setdefault("SB_CACHE", os.path.join(REPO, "cache", "mask"))
from zio import RemoteArray  # noqa
VOLID = {"PHerc0125": "20250821151825", "PHerc0191": "20250821151635", "PHerc0211": "20250821151803", "PHerc0257": "20250821151750",
         "PHerc0268": "20251110183117", "PHerc0358": "20250821151737", "PHerc0813": "20250821151723", "PHerc0826": "20250821151701", "PHerc0139": "20250728140407"}
os.makedirs(os.path.join(REPO, "figs", "giant_viz"), exist_ok=True); RES = {}
for SC, VID in VOLID.items():
    c = json.load(open(os.path.join(REPO, "results", f"census_{SC}.json"))); Z = c["ribbons"][1]["z"]; p = c["params"]; ds = p["downsample_factor"]
    base = f"https://dl.ash2txt.org/datasets/spiral_datasets/{SC}/{VID}/tracks/"
    ej = re.findall(r'href="([^"]*extract\.json)"', urllib.request.urlopen(base, timeout=60).read().decode())[0]
    meta = json.loads(urllib.request.urlopen(base + ej, timeout=60).read()); (y0, x0), (y1, x1) = meta["yx_bounds"]
    url = meta["array_url"].replace("s3://vesuvius-challenge-open-data", "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com")
    arr = RemoteArray(url.rsplit("/", 1)[0], url.rsplit("/", 1)[1])
    pred = (arr.read(Z, Z + p["z_chunk_depth_h"], y0, y1 + 1, x0, x1 + 1, workers=32) > 0).astype(np.uint8)
    cc = cc3d.connected_components(pred, connectivity=6); cs = tuple(s - s % ds for s in cc.shape)
    lab = cc[:cs[0], :cs[1], :cs[2]].reshape(cs[0] // ds, ds, cs[1] // ds, ds, cs[2] // ds, ds).max((1, 3, 5))[0]
    ids, cnt = np.unique(lab[lab > 0], return_counts=True); giant = np.isin(lab, ids[cnt > p["max_area_threshold"] // ds ** 3]); kept = (lab > 0) & ~giant
    body = binary_dilation(binary_fill_holes(binary_closing(kept, iterations=8)), iterations=3)
    out = float((giant & ~body).sum() / max(giant.sum(), 1))
    img = np.zeros(lab.shape + (3,), np.uint8); img[kept] = (150, 150, 150); img[giant] = (220, 50, 40)
    im = Image.fromarray(img); im.thumbnail((500, 500)); im.save(os.path.join(REPO, "figs", "giant_viz", f"{SC}_z{Z}.png"))
    RES[SC] = dict(z=int(Z), giant_share_of_fg=round(float(giant.sum() / max((lab > 0).sum(), 1)), 3), giant_share_outside_body=round(out, 3))
    print(f"{SC} z={Z}: giant share of fg {giant.sum() / max((lab > 0).sum(), 1):.3f}; share of giant OUTSIDE scroll body {out:.3f}", flush=True)
json.dump(RES, open(os.path.join(REPO, "results", "giant_where.json"), "w"), indent=1)
