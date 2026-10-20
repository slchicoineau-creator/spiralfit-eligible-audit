"""Resample PHerc0139's lasagna normals (2.399 um volume 20260102150214, group 4 = scale 16) onto the fit frame (9.362 um
volume 20250728140407, group 2 = lasagna_scale 4), for the band +-64 slices (prereg section 3.3), then run gate S0-b.
CPU only. Run from the repo root:  .venv/bin/python -I work/spiral_accuracy/scripts/map_normals.py [fetch|map|check|all]"""
import json, os, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import numcodecs
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from acc_lib import D, ROOT, BAND, SCORE_Z, load_mesh, lines

S3 = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
RUN = "PHerc0139/representations/predictions/lasagna/20260102150214-lasagna-20260419180421-L2/PHerc0139-20260102150214-lasagna-20260724"
SRC_G, SRC_S = "4", 16
TGT_S, TGT_SHAPE_FULL = 4, (20974, 6621, 6621)
TGT_SHAPE = tuple(int(np.ceil(s / TGT_S)) for s in TGT_SHAPE_FULL)
MARGIN = 64
CH = 32
SRC_CACHE = os.path.join(D, "lasagna_src_g4")
DS = os.path.join(D, "PHerc0139_20250728140407")
OUTJ = os.path.join(ROOT, "work/spiral_accuracy/out/s0_normals.json")
cat = json.load(open(os.path.join(D, "catalog_0139.json")))
M = np.array(cat["properties"]["20260102150214"]["transforms"][0]["transformation_matrix"], float)   # 2.399 xyz -> 9.362 xyz
L, t = M[:, :3], M[:, 3]
Linv = np.linalg.inv(L)
NT = np.linalg.inv(L).T            # normals map with the inverse transpose
blosc = numcodecs.Blosc(cname="lz4", clevel=5, shuffle=1)


def get(u):
    for _ in range(8):
        try:
            return urllib.request.urlopen(u, timeout=120).read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return b""
        except Exception:
            pass
    raise RuntimeError(u)


def tgt_iz_range():
    return (BAND[0] - MARGIN) // TGT_S, -(-(BAND[1] + MARGIN) // TGT_S)


def src_z_range():
    """source group-4 z index range touched by the target planes (all xy corners)."""
    iz0, iz1 = tgt_iz_range()
    zs = []
    for iz in (iz0, iz1):
        for y in (0, TGT_SHAPE_FULL[1]):
            for x in (0, TGT_SHAPE_FULL[2]):
                p2 = Linv @ (np.array([x, y, iz * TGT_S], float) - t)
                zs.append(p2[2] / SRC_S)
    return int(np.floor(min(zs))) - 1, int(np.ceil(max(zs))) + 2


def src_shape():
    za = json.loads(get(f"{S3}/{RUN}_nx.ome.zarr/{SRC_G}/.zarray"))
    return tuple(za["shape"]), za


def fetch():
    shp, za = src_shape()
    assert za["chunks"] == [CH, CH, CH] and za["dtype"] == "|u1", za
    z0, z1 = src_z_range()
    jobs = []
    for ch in ("nx", "ny"):
        for cz in range(z0 // CH, -(-z1 // CH)):
            for cy in range(-(-shp[1] // CH)):
                for cx in range(-(-shp[2] // CH)):
                    jobs.append((ch, cz, cy, cx))

    def one(j):
        ch, cz, cy, cx = j
        p = os.path.join(SRC_CACHE, ch, f"{cz}.{cy}.{cx}")
        if os.path.exists(p):
            return
        b = get(f"{S3}/{RUN}_{ch}.ome.zarr/{SRC_G}/{cz}/{cy}/{cx}")
        os.makedirs(os.path.dirname(p), exist_ok=True); open(p + ".part", "wb").write(b); os.replace(p + ".part", p)
    with ThreadPoolExecutor(48) as ex:
        list(ex.map(one, jobs))
    json.dump(dict(shape=shp, src_z=[z0, z1], n_chunks=len(jobs)), open(os.path.join(SRC_CACHE, "fetch.json"), "w"))
    print("fetched", len(jobs), "chunk files; src z", z0, z1)


def load_src():
    info = json.load(open(os.path.join(SRC_CACHE, "fetch.json")))
    shp = info["shape"]; z0, z1 = info["src_z"]
    cz0, cz1 = z0 // CH, -(-z1 // CH)
    out = {}
    for ch in ("nx", "ny"):
        A = np.zeros(((cz1 - cz0) * CH, shp[1], shp[2]), np.uint8)
        for cz in range(cz0, cz1):
            for cy in range(-(-shp[1] // CH)):
                for cx in range(-(-shp[2] // CH)):
                    b = open(os.path.join(SRC_CACHE, ch, f"{cz}.{cy}.{cx}"), "rb").read()
                    if not b:
                        continue
                    blk = np.frombuffer(blosc.decode(b), np.uint8).reshape(CH, CH, CH)
                    ys, xs = cy * CH, cx * CH
                    h, w = min(CH, shp[1] - ys), min(CH, shp[2] - xs)
                    A[(cz - cz0) * CH:(cz - cz0 + 1) * CH, ys:ys + h, xs:xs + w] = blk[:, :h, :w]
        out[ch] = A
    return out, cz0 * CH


def decode(nxu, nyu):
    nx = (nxu.astype(np.float32) - 128.0) / 127.0; ny = (nyu.astype(np.float32) - 128.0) / 127.0
    nz = np.sqrt(np.clip(1.0 - nx * nx - ny * ny, 0.0, 1.0))
    return nx, ny, nz


def encode(nx, ny, nz):
    flip = np.where(nz < 0.0, -1.0, 1.0)
    e = lambda v: np.clip(np.round(v * flip * 127.0 + 128.0), 0, 255).astype(np.uint8)
    return e(nx), e(ny)


def map_store():
    src, zoff = load_src()
    iz0, iz1 = tgt_iz_range()
    ny_t, nx_t = TGT_SHAPE[1], TGT_SHAPE[2]
    tgt = {c: np.zeros((iz1 - iz0, ny_t, nx_t), np.uint8) for c in ("nx", "ny")}
    gy, gx = np.meshgrid(np.arange(ny_t) * TGT_S, np.arange(nx_t) * TGT_S, indexing="ij")
    S = src["nx"].shape
    for iz in range(iz0, iz1):
        P = np.stack([gx.ravel(), gy.ravel(), np.full(gx.size, iz * TGT_S)], 0).astype(np.float64)     # xyz
        Q = Linv @ (P - t[:, None])
        ix, iy, izs = (np.round(Q / SRC_S).astype(np.int64))
        izs = izs - zoff
        ok = (ix >= 0) & (ix < S[2]) & (iy >= 0) & (iy < S[1]) & (izs >= 0) & (izs < S[0])
        a = np.zeros(gx.size, np.uint8); b = np.zeros(gx.size, np.uint8)
        a[ok] = src["nx"][izs[ok], iy[ok], ix[ok]]; b[ok] = src["ny"][izs[ok], iy[ok], ix[ok]]
        valid = ok & ~((a == 0) & (b == 0))
        nx, ny, nz = decode(a[valid], b[valid])
        v = NT @ np.stack([nx, ny, nz], 0).astype(np.float64)
        v /= np.linalg.norm(v, axis=0, keepdims=True) + 1e-12
        ea, eb = encode(v[0], v[1], v[2])
        oa = np.zeros(gx.size, np.uint8); ob = np.zeros(gx.size, np.uint8)
        oa[valid] = ea; ob[valid] = eb
        oa[valid & (oa == 0) & (ob == 0)] = 1          # never write a valid normal as the no-data code
        tgt["nx"][iz - iz0] = oa.reshape(ny_t, nx_t); tgt["ny"][iz - iz0] = ob.reshape(ny_t, nx_t)
    for ch in ("nx", "ny"):
        root = os.path.join(DS, "lasagna_inputs", f"las_008_{ch}.ome.zarr")
        g = os.path.join(root, "2"); os.makedirs(g, exist_ok=True)
        json.dump({"zarr_format": 2}, open(os.path.join(root, ".zgroup"), "w"))
        json.dump({"multiscales": [{"version": "0.4", "name": ch, "axes": [{"name": a, "type": "space", "unit": "pixel"} for a in "zyx"],
                                    "datasets": [{"path": "2", "coordinateTransformations": [{"type": "scale", "scale": [4.0, 4.0, 4.0]}]}]}],
                   "spiral_accuracy_provenance": dict(source=f"s3://vesuvius-challenge-open-data/{RUN}_{ch}.ome.zarr/{SRC_G}",
                                                      transform="catalog 20260102150214 transforms[0]", method="nearest, inverse-transpose normal map",
                                                      z_planes=[iz0, iz1])}, open(os.path.join(root, ".zattrs"), "w"), indent=1)
        json.dump({"shape": list(TGT_SHAPE), "chunks": [CH, CH, CH], "dtype": "|u1", "fill_value": 0, "order": "C", "filters": None,
                   "dimension_separator": "/", "compressor": {"id": "blosc", "cname": "lz4", "clevel": 5, "shuffle": 1, "blocksize": 0},
                   "zarr_format": 2}, open(os.path.join(g, ".zarray"), "w"), indent=1)
        A = tgt[ch]; nw = 0
        for cz in range(iz0 // CH, -(-iz1 // CH)):
            blkz = np.zeros((CH, ny_t, nx_t), np.uint8)
            lo, hi = max(cz * CH, iz0), min((cz + 1) * CH, iz1)
            blkz[lo - cz * CH:hi - cz * CH] = A[lo - iz0:hi - iz0]
            for cy in range(-(-ny_t // CH)):
                for cx in range(-(-nx_t // CH)):
                    blk = np.zeros((CH, CH, CH), np.uint8)
                    sub = blkz[:, cy * CH:(cy + 1) * CH, cx * CH:(cx + 1) * CH]
                    blk[:, :sub.shape[1], :sub.shape[2]] = sub
                    if not blk.any():
                        continue
                    p = os.path.join(g, str(cz), str(cy)); os.makedirs(p, exist_ok=True)
                    open(os.path.join(p, str(cx)), "wb").write(blosc.encode(np.ascontiguousarray(blk)))
                    nw += 1
        print(ch, "chunks written", nw)
    np.savez_compressed(os.path.join(D, "mapped_g2_band.npz"), nx=tgt["nx"], ny=tgt["ny"], iz0=iz0)


def mesh_normals(d):
    m = load_mesh(d); C = m["C"]
    n = np.cross(np.gradient(C, axis=0), np.gradient(C, axis=1))
    n /= np.linalg.norm(n, axis=-1, keepdims=True) + 1e-12
    ok = np.isfinite(C).all(-1) & np.isfinite(n).all(-1)
    return C[ok], n[ok]                     # zyx


def check():
    import glob
    src, zoff = load_src()
    tg = np.load(os.path.join(D, "mapped_g2_band.npz")); iz0 = int(tg["iz0"])
    res = {}
    for w in ("w025", "w031", "w037", "w043", "w049", "w055"):
        d9 = glob.glob(os.path.join(D, "labels_9um", w + "_*"))[0]; d2 = glob.glob(os.path.join(D, "labels_2um", w + "_*"))[0]
        P9, N9 = mesh_normals(d9)
        sel = (P9[:, 0] >= SCORE_Z[0]) & (P9[:, 0] < SCORE_Z[1]); P9, N9 = P9[sel], N9[sel]
        idx = np.round(P9 / TGT_S).astype(int); idx[:, 0] -= iz0
        a = tg["nx"][idx[:, 0], idx[:, 1], idx[:, 2]]; b = tg["ny"][idx[:, 0], idx[:, 1], idx[:, 2]]
        v = (a > 0) | (b > 0); nx, ny, nz = decode(a[v], b[v])
        cm = np.abs(nx * N9[v, 2] + ny * N9[v, 1] + nz * N9[v, 0])
        P2, N2 = mesh_normals(d2)
        q = (M @ np.concatenate([P2[:, ::-1], np.ones((len(P2), 1))], 1).T).T     # mapped xyz, for the same z selection
        sel2 = (q[:, 2] >= SCORE_Z[0]) & (q[:, 2] < SCORE_Z[1]); P2, N2 = P2[sel2], N2[sel2]
        i2 = np.round(P2 / SRC_S).astype(int); i2[:, 0] -= zoff
        okk = (i2[:, 0] >= 0) & (i2[:, 0] < src["nx"].shape[0])
        a2 = src["nx"][i2[okk, 0], i2[okk, 1], i2[okk, 2]]; b2 = src["ny"][i2[okk, 0], i2[okk, 1], i2[okk, 2]]
        v2 = (a2 > 0) | (b2 > 0); nx2, ny2, nz2 = decode(a2[v2], b2[v2]); N2o = N2[okk][v2]
        cn = np.abs(nx2 * N2o[:, 2] + ny2 * N2o[:, 1] + nz2 * N2o[:, 0])
        res[w] = dict(mapped_n=int(v.sum()), mapped_median=round(float(np.median(cm)), 4), mapped_share_gt09=round(float((cm > 0.9).mean()), 4),
                      native_n=int(v2.sum()), native_median=round(float(np.median(cn)), 4), native_share_gt09=round(float((cn > 0.9).mean()), 4),
                      mapped_nodata_share=round(float(1 - v.mean()), 4))
        res.setdefault("_all_m", []).extend(cm.tolist()); res.setdefault("_all_n", []).extend(cn.tolist())
    am, an = np.array(res.pop("_all_m")), np.array(res.pop("_all_n"))
    res["pooled"] = dict(mapped_median=round(float(np.median(am)), 4), native_median=round(float(np.median(an)), 4),
                         PASS=bool(np.median(am) >= np.median(an) - 0.03))
    json.dump(res, open(OUTJ, "w"), indent=1); print(json.dumps(res, indent=1))


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("fetch", "all"):
        fetch()
    if what in ("map", "all"):
        map_store()
    if what in ("check", "all"):
        check()
