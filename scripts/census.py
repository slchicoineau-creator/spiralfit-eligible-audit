"""W2 TrackFix day-2 census (operationalisation pre-registered in the work log before running). For each scroll, re-runs
the published extractor's horizontal-ribbon steps on 3 ribbons and measures the dropped-foreground share and the
junction-adjacent-or-merged share of track vertices. Usage: python census.py SCROLL  (writes census_<SCROLL>.json)"""
import json, os, sys, urllib.request, numpy as np, cc3d, kimimaro, networkx as nx
from scipy.ndimage import distance_transform_edt
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib")); os.environ.setdefault("SB_CACHE", os.path.join(REPO, "cache", "mask"))
from zio import RemoteArray  # noqa
def main():
    SC = sys.argv[1]
    VOLID = {"PHerc0125": "20250821151825", "PHerc0191": "20250821151635", "PHerc0211": "20250821151803", "PHerc0257": "20250821151750",
             "PHerc0268": "20251110183117", "PHerc0358": "20250821151737", "PHerc0813": "20250821151723", "PHerc0826": "20250821151701",
             "PHerc1203": "20250820131727", "PHerc1545": "20250821151648", "PHerc0139": "20250728140407"}
    base = f"https://dl.ash2txt.org/datasets/spiral_datasets/{SC}/{VOLID[SC]}/tracks/"
    listing = urllib.request.urlopen(base, timeout=60).read().decode()
    ej = [h for h in __import__("re").findall(r'href="([^"]*extract\.json)"', listing)]
    if not ej:
        json.dump(dict(scroll=SC, error="no track store"), open(os.path.join(REPO, "results", f"census_{SC}.json"), "w")); return print(f"{SC}: no track store")
    meta = json.loads(urllib.request.urlopen(base + ej[0], timeout=60).read()); p = meta["parameters"]
    url = meta["array_url"].replace("s3://vesuvius-challenge-open-data", "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com")
    arr = RemoteArray(url.rsplit("/", 1)[0], url.rsplit("/", 1)[1])
    ds = p["downsample_factor"]; scale = ds ** 3; zlo, zhi = meta["effective_z_bounds"]; (y0, x0), (y1, x1) = meta["yx_bounds"]


    def pool(x):
        c = tuple(s - s % ds for s in x.shape); x = x[:c[0], :c[1], :c[2]]
        return x.reshape(c[0] // ds, ds, c[1] // ds, ds, c[2] // ds, ds).max((1, 3, 5))


    rows = []
    for f in (0.25, 0.5, 0.75):
        z = zlo + int(round(f * (zhi - zlo) / 8)) * 8
        pred = (arr.read(z, z + p["z_chunk_depth_h"], y0, y1 + 1, x0, x1 + 1, workers=32) > 0).astype(np.uint8)
        cc = cc3d.connected_components(pred, connectivity=6)
        lab = pool(cc); fg = lab > 0
        ids, cnt = np.unique(lab[fg], return_counts=True)
        big = set(ids[cnt > p["max_area_threshold"] // scale]); small = set(ids[cnt < max(1, p["dust_threshold"] // scale)])
        dropped = float(cnt[np.isin(ids, list(big))].sum() / cnt.sum()) if cnt.sum() else 0.0
        keep = lab.copy(); keep[np.isin(keep, list(big | small))] = 0
        sk = kimimaro.skeletonize(keep, teasar_params={"scale": 1.0, "const": 2.0}, anisotropy=(ds, ds, ds), dust_threshold=p["dust_threshold"],
                                  fix_branching=True, fix_borders=True, fill_holes=False, parallel=12, parallel_chunk_size=250, in_place=False)
        diam = 2 * distance_transform_edt(keep[0] > 0)                       # in-slice diameter (pooled px); ribbon pools to 1 slice
        verts, jadj, dvals = 0, 0, []
        per_track = []
        for s in sk.values():
            g = nx.Graph(); g.add_edges_from(s.edges); deg = dict(g.degree())
            crit = [n for n in g.nodes() if deg[n] != 2]; seen = set()
            for n0 in crit:
                for nb in g.neighbors(n0):
                    e = tuple(sorted((n0, nb)))
                    if e in seen:
                        continue
                    path = [n0, nb]; seen.add(e); cur = nb
                    while deg[cur] == 2:
                        nx_ = [m for m in g.neighbors(cur) if m != path[-2]]
                        if not nx_:
                            break
                        e = tuple(sorted((cur, nx_[0])))
                        if e in seen:
                            break
                        path.append(nx_[0]); seen.add(e); cur = nx_[0]
                    if len(path) < 10:
                        continue
                    v = (s.vertices[path] / ds).astype(int); v[:, 1] = np.clip(v[:, 1], 0, diam.shape[0] - 1); v[:, 2] = np.clip(v[:, 2], 0, diam.shape[1] - 1)
                    d = diam[v[:, 1], v[:, 2]]
                    ja = deg[path[0]] >= 3 or deg[path[-1]] >= 3
                    per_track.append((len(path), ja, d)); dvals.append(d)
        rows.append(dict(z=z, fg=int(cnt.sum()), n_big=len(big), dropped=dropped, per_track=per_track))
        print(f"{SC} z={z}: fg {int(cnt.sum())} big comps {len(big)} dropped {dropped:.3f} tracks {len(per_track)}", flush=True)
    allv = np.concatenate([d for r in rows for (_, _, d) in r["per_track"]]) if any(r["per_track"] for r in rows) else np.array([1.0])
    thr = 1.5 * float(np.median(allv))
    tot = un = ja_v = mg_v = 0
    for r in rows:
        for n, ja, d in r["per_track"]:
            m = d > thr; tot += n; ja_v += n if ja else 0; mg_v += int(m.sum()); un += n if ja else int(m.sum())
    fgs = sum(r["fg"] for r in rows); drop = sum(r["dropped"] * r["fg"] for r in rows) / fgs if fgs else 0.0
    out = dict(scroll=SC, params=p, median_diam=float(np.median(allv)), merged_thr=thr, track_vertices=tot, junction_share=ja_v / tot if tot else 0,
               merged_share=mg_v / tot if tot else 0, union_share=un / tot if tot else 0, dropped_share=drop,
               ribbons=[dict(z=r["z"], fg=r["fg"], n_big=r["n_big"], dropped=r["dropped"], n_tracks=len(r["per_track"])) for r in rows])
    json.dump(out, open(os.path.join(REPO, "results", f"census_{SC}.json"), "w"), indent=1)
    print(f"{SC}: union {out['union_share']:.3f} (junction {out['junction_share']:.3f}, merged {out['merged_share']:.3f}), dropped {drop:.3f}")


if __name__ == "__main__":
    main()
