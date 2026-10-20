"""W2 descriptive census of VERTICAL slabs (advisor check 2; not a gate): dropped-foreground share on 3 vy slabs of one
scroll, replicating extract_vertical + prepare_cc_labels with the scroll's published extract.json parameters.
Usage: python census_v.py SCROLL VOLID"""
import json, os, sys, urllib.request, re, numpy as np, cc3d
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib")); os.environ.setdefault("SB_CACHE", os.path.join(REPO, "cache", "mask"))
from zio import RemoteArray  # noqa


def main():
    SC, VID = sys.argv[1], sys.argv[2]
    base = f"https://dl.ash2txt.org/datasets/spiral_datasets/{SC}/{VID}/tracks/"
    ej = re.findall(r'href="([^"]*extract\.json)"', urllib.request.urlopen(base, timeout=60).read().decode())[0]
    meta = json.loads(urllib.request.urlopen(base + ej, timeout=60).read()); p = meta["parameters"]; ds = p["downsample_factor"]
    url = meta["array_url"].replace("s3://vesuvius-challenge-open-data", "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com")
    arr = RemoteArray(url.rsplit("/", 1)[0], url.rsplit("/", 1)[1])
    zlo, zhi = max(p["z_min"], 0), min(p["z_max"], arr.shape[0]); (y0, x0), (y1, x1) = meta["yx_bounds"]
    out = []
    for f in (0.35, 0.5, 0.65):
        w = y0 + int(round(f * (y1 - y0) / p["yx_stride_v"])) * p["yx_stride_v"]
        pred = (arr.read(zlo, zhi, w, w + p["yx_chunk_thickness_v"], 0, arr.shape[2], workers=32) > 0).astype(np.uint8)
        cc = cc3d.connected_components(pred, connectivity=6)
        c = tuple(s - s % ds for s in cc.shape)
        lab = cc[:c[0], :c[1], :c[2]].reshape(c[0] // ds, ds, c[1] // ds, ds, c[2] // ds, ds).max((1, 3, 5))
        ids, cnt = np.unique(lab[lab > 0], return_counts=True)
        big = cnt > p["max_area_threshold"] // ds ** 3
        d = float(cnt[big].sum() / cnt.sum()) if cnt.sum() else 0.0
        out.append(dict(w=w, fg=int(cnt.sum()), n_big=int(big.sum()), dropped=d, largest=int(cnt.max()) if len(cnt) else 0))
        print(f"{SC} vy:{w}: fg {int(cnt.sum())} big comps {int(big.sum())} dropped {d:.3f} largest {out[-1]['largest']}", flush=True)
    json.dump(out, open(os.path.join(REPO, "results", f"census_v_{SC}.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
