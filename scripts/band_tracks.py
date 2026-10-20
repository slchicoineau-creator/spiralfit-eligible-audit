"""TrackFix band builder (runs in WSL, spiral-fitting venv). For one scroll and z band [Z0, Z1) it writes three track
DBMs with the published key layout (h:<z>, vy:<w>, vx:<w>; values = pickled lists of (N,3) int32 zyx arrays):
  A  published : the team's published tracks lying entirely inside [Z0-M, Z1+M) (read from the public .vctracks by
                 HTTP range requests; nothing re-extracted)
  B  quickfix  : A + tracks cut from every connected component of the band's own ribbons/slabs that NO published
                 track touches (the components the size cap dropped), with the extractor's own interjoint rule
  C  trackfix  : A + the same uncovered components, but tracks continue straight through junctions (branch pairs
                 within 30 degrees of straight) and stretches whose mask diameter exceeds 1.5x the scroll's median
                 published-track diameter at native resolution (fused sheets) are cut out; pieces < 10 vertices dropped
Band-local slabs/ribbons are processed with the scroll's published extract parameters. Only tracks entirely inside
[Z0-M, Z1+M) are kept in every arm (the fitter keeps only tracks wholly inside its ROI).
Usage: python band_tracks.py SCROLL VOLID Z0 Z1 OUTDIR"""
import dbm, json, os, pickle, re, sys, urllib.request, time
import numpy as np, cc3d, kimimaro, networkx as nx
from concurrent.futures import ThreadPoolExecutor
from scipy.ndimage import distance_transform_edt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
os.environ.setdefault("SB_CACHE", os.path.join(os.environ.get("R3", os.path.expanduser("~/r3")), "/maskcache".lstrip("/")))
from zio import RemoteArray  # noqa
TILE = 256; M = 64; STRAIGHT = np.cos(np.deg2rad(30)); FUSED = 1.5; MINV = 10


def get(u, rng=None):
    for a in range(8):
        try:
            req = urllib.request.Request(u, headers={"Range": f"bytes={rng[0]}-{rng[1] - 1}"} if rng else {})
            return urllib.request.urlopen(req, timeout=300).read()
        except Exception:
            time.sleep(2 ** min(a, 5))
    raise RuntimeError(u)


def published(SC, VID, zlo, zhi, cache):
    base = f"https://dl.ash2txt.org/datasets/spiral_datasets/{SC}/{VID}/tracks/"
    name = re.findall(r'href="([^"]*\.dbm\.vctracks)/"', get(base).decode())[0]; V = base + name + "/"
    os.makedirs(cache, exist_ok=True)
    arr = {}
    for f, dt in (("offsets.i64", "<i8"), ("family_codes.i8", "<i1"), ("z_bounds.i32", "<i4"), ("source_ids.u64", "<u8")):
        p = os.path.join(cache, f)
        if not os.path.exists(p):
            open(p, "wb").write(get(V + f))
        arr[f] = np.fromfile(p, dt)
    off, fam, zb = arr["offsets.i64"], arr["family_codes.i8"], arr["z_bounds.i32"].reshape(-1, 2)
    sel = np.nonzero((zb[:, 0] >= zlo) & (zb[:, 1] < zhi) & (fam >= 0))[0]
    # merge byte ranges of selected tracks (coordinates.i32, 12 bytes per point)
    starts, ends = off[sel] * 12, off[sel + 1] * 12
    order = np.argsort(starts); starts, ends, sel = starts[order], ends[order], sel[order]
    blocks = []; bs, be = starts[0], ends[0]
    for s, e in zip(starts[1:], ends[1:]):
        if s - be < 4_000_000 and e - bs < 64_000_000:
            be = max(be, e)
        else:
            blocks.append((bs, be)); bs, be = s, e
    blocks.append((bs, be))
    with ThreadPoolExecutor(16) as ex:
        data = list(ex.map(lambda b: (b, get(V + "coordinates.i32", b)), blocks))
    buf = {b: np.frombuffer(d, "<i4").reshape(-1, 3) for b, d in data}
    bi = 0; tracks = []
    bl = sorted(buf)
    for i, s, e in zip(sel, starts, ends):
        while not (bl[bi][0] <= s and e <= bl[bi][1]):
            bi += 1
        b0 = bl[bi][0]; tracks.append((int(fam[i]), int(arr["source_ids.u64"][i] >> np.uint64(32)), buf[bl[bi]][(s - b0) // 12:(e - b0) // 12].copy()))
    return tracks


def all_keys(p, zlimit, y0, y1, x0, x1):
    """The extractor's DBM key set (sorted order = packed-store key ordinals). Horizontal keys past the volume end are
    never written; vertical keys are written for every stride position in the yx range."""
    ks = [f"h:{z}" for z in range(p["z_min"], p["z_max"], p["z_chunk_stride_h"]) if min(z + p["z_chunk_depth_h"], p["z_max"], zlimit) > z]
    ks += [f"vy:{w}" for w in range(y0, y1, p["yx_stride_v"])] + [f"vx:{w}" for w in range(x0, x1, p["yx_stride_v"])]
    return sorted(k.encode() for k in ks)


def key_of(fam, t, p, y0, x0):
    if fam == 0:
        return f"h:{int(t[:, 0].min()) - (int(t[:, 0].min()) - p['z_min']) % p['z_chunk_stride_h']}"
    ys, xs = t[:, 1], t[:, 2]
    if ys.max() - ys.min() < p["yx_chunk_thickness_v"] + 1:          # constant y -> vy slab
        w = int(ys.min()); return f"vy:{w - (w - y0) % p['yx_stride_v']}"
    w = int(xs.min()); return f"vx:{w - (w - x0) % p['yx_stride_v']}"


def paths_interjoint(g):
    deg = dict(g.degree()); out = []; seen = set()
    for n0 in [n for n in g.nodes() if deg[n] != 2]:
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
            out.append(path)
    return out


def resolve(g, verts):
    """TrackFix: join interjoint branches that continue straight (< 30 deg) through a junction."""
    br = paths_interjoint(g); deg = dict(g.degree())
    ends = {}                                   # junction node -> list of (branch idx, which end)
    for i, b in enumerate(br):
        for which, n in ((0, b[0]), (1, b[-1])):
            if deg[n] >= 3:
                ends.setdefault(n, []).append((i, which))
    def direction(i, which):
        b = br[i]; seq = b if which == 0 else b[::-1]; k = min(8, len(seq) - 1)
        v = verts[seq[k]] - verts[seq[0]]; n = np.linalg.norm(v); return v / n if n else v
    link = {}
    for n, lst in ends.items():
        cand = []
        for a in range(len(lst)):
            for b in range(a + 1, len(lst)):
                c = float(direction(*lst[a]) @ direction(*lst[b]))
                if c < -STRAIGHT:
                    cand.append((c, lst[a], lst[b]))
        used = set()
        for c, a, b in sorted(cand):
            if a in used or b in used:
                continue
            used |= {a, b}; link[a] = b; link[b] = a
    done = set(); tracks = []
    for i in range(len(br)):
        if i in done:
            continue
        # walk to one free end of the chain containing branch i
        cur, which = i, 0
        seen_c = {i}
        while (cur, which) in link:
            j, wj = link[(cur, which)]
            if j in seen_c:
                break
            seen_c.add(j); cur, which = j, 1 - wj
        chain = []; c, w = cur, which; visited = set()
        while True:
            visited.add(c); seg = br[c] if w == 0 else br[c][::-1]
            chain.extend(seg if not chain else seg[1:])
            nxt = link.get((c, 1 - w))
            if nxt is None or nxt[0] in visited:
                break
            c, w = nxt[0], nxt[1]
        done |= visited; tracks.append(chain)
    return tracks


def pool_labels(cc, ds):
    c = tuple(s - s % ds for s in cc.shape)
    return cc[:c[0], :c[1], :c[2]].reshape(c[0] // ds, ds, c[1] // ds, ds, c[2] // ds, ds).max((1, 3, 5))


def native_diam(mask_native, thin):
    return 2 * distance_transform_edt(mask_native.max(axis=thin))


def recover(sub, cov_pts, p, offset, thin, med_diam):
    """sub: native binary slab/ribbon. cov_pts: published track vertices of this key, relative to sub (native zyx).
    Returns (arm B tracks, arm C tracks) from components no published track touches; one skeletonisation."""
    ds = p["downsample_factor"]
    cc = cc3d.connected_components(sub, connectivity=6)
    lab = pool_labels(cc, ds)
    ids = np.unique(lab[lab > 0])
    touched = set()
    ax = [i for i in range(3) if i != thin]
    if len(cov_pts):
        q = cov_pts // ds
        ok = np.all((q >= 0) & (q < np.array(lab.shape)), 1); q = q[ok]
        for dy in (-1, 0, 1):                       # a published vertex touches its pooled cell or an in-plane neighbour
            for dx in (-1, 0, 1):
                qq = q.copy()
                qq[:, ax[0]] = np.clip(qq[:, ax[0]] + dy, 0, lab.shape[ax[0]] - 1)
                qq[:, ax[1]] = np.clip(qq[:, ax[1]] + dx, 0, lab.shape[ax[1]] - 1)
                for t in range(lab.shape[thin]):
                    qq[:, thin] = t; touched.update(np.unique(lab[qq[:, 0], qq[:, 1], qq[:, 2]]).tolist())
    unc = np.setdiff1d(ids, np.array(sorted(touched - {0}), dtype=ids.dtype))
    if not len(unc):
        return [], []
    keep = np.where(np.isin(lab, unc), lab, 0)
    # runtime: cut recovered components into TILE x TILE pooled tiles (1-px zero grid lines) so TEASAR tracing stays
    # tractable on the giant merged components; identical for arms B and C (logged pre-data)
    for a in [i for i in range(3) if i != thin]:
        idx = [slice(None)] * 3; idx[a] = slice(TILE, None, TILE); keep[tuple(idx)] = 0
    keep = cc3d.connected_components(keep > 0, connectivity=6).astype(np.uint32)
    cc3d.dust(keep, threshold=max(1, p["dust_threshold"] // ds ** 3), in_place=True, precomputed_ccl=True)
    if not keep.any():
        return [], []
    dn = native_diam(np.isin(cc, unc), thin)
    sk = kimimaro.skeletonize(keep, teasar_params={"scale": 1.0, "const": 2.0}, anisotropy=(ds, ds, ds),
                              dust_threshold=p["dust_threshold"], fix_branching=True, fix_borders=True, fill_holes=False,
                              parallel=1, in_place=True, progress=False)
    off = np.asarray(offset, np.int64); outB, outC = [], []
    for s in sk.values():
        g = nx.Graph(); g.add_edges_from(s.edges)
        if g.number_of_nodes() == 0:
            continue
        for path in paths_interjoint(g):
            if len(path) >= MINV:
                outB.append((s.vertices[path].astype(np.int64) + off).astype(np.int32))
        for path in resolve(g, s.vertices):
            v = s.vertices[path]; q = v.astype(int)
            a0 = np.clip(q[:, ax[0]], 0, dn.shape[0] - 1); a1 = np.clip(q[:, ax[1]], 0, dn.shape[1] - 1)
            ok = dn[a0, a1] <= FUSED * med_diam
            cur = []
            for vi, o in zip(v, ok):
                if o:
                    cur.append(vi)
                else:
                    if len(cur) >= MINV:
                        outC.append((np.array(cur).astype(np.int64) + off).astype(np.int32))
                    cur = []
            if len(cur) >= MINV:
                outC.append((np.array(cur).astype(np.int64) + off).astype(np.int32))
    return outB, outC


def read_zslices(arr, zs, depth, y0, y1, x0, x1):
    """Read only z-slices [z, z+depth) for z in zs (all inside one chunk row), decompressing each chunk once."""
    cz, cy, cx = arr.chunks; zr = zs[0] // cz
    out = np.zeros((len(zs), depth, y1 - y0, x1 - x0), bool)
    idx = np.concatenate([np.arange(z, z + depth) for z in zs]) - zr * cz
    def one(j):
        iy, ix = j; ch = arr._chunk((zr, iy, ix))[idx]
        ya, yb = max(iy * cy, y0), min((iy + 1) * cy, y1); xa, xb = max(ix * cx, x0), min((ix + 1) * cx, x1)
        if ya < yb and xa < xb:
            out.reshape(-1, y1 - y0, x1 - x0)[:, ya - y0:yb - y0, xa - x0:xb - x0] = ch[:, ya - iy * cy:yb - iy * cy, xa - ix * cx:xb - ix * cx] > 0
    jobs = [(iy, ix) for iy in range(y0 // cy, (y1 - 1) // cy + 1) for ix in range(x0 // cx, (x1 - 1) // cx + 1)]
    with ThreadPoolExecutor(4) as ex:
        list(ex.map(one, jobs))
    return out


G = {}


def task(t):
    """One mask block (h chunk-row, vy chunk-row, vx chunk-col) -> recovered tracks for every ribbon/slab inside it."""
    kind, lo, hi = t; p = G["p"]; arr = G["arr"]; zlo, zhi = G["zlo"], G["zhi"]; (y0, x0), (y1, x1) = G["yx"]
    th, dh = p["yx_chunk_thickness_v"], p["z_chunk_depth_h"]
    if kind == "h":
        zs = [z for z in G["hz"] if lo <= z and z + dh <= hi]
        if not zs:
            return []
        blk = read_zslices(arr, zs, dh, y0, y1 + 1, x0, x1 + 1)
        subs = [(f"h:{z}", (i,), (z, y0, x0), 0) for i, z in enumerate(zs)]
    elif kind == "vy":
        blk = arr.read(zlo, zhi, lo, hi, x0, x1 + 1, workers=4) > 0
        subs = [(f"vy:{w}", (slice(None), slice(w - lo, w - lo + th), slice(None)), (zlo, w, x0), 1) for w in G["vy"] if lo <= w and w + th <= hi]
    else:
        blk = arr.read(zlo, zhi, y0, y1 + 1, lo, hi, workers=4) > 0
        subs = [(f"vx:{w}", (slice(None), slice(None), slice(w - lo, w - lo + th)), (zlo, y0, w), 2) for w in G["vx"] if lo <= w and w + th <= hi]
    out = []
    inside = lambda tr: tr[:, 0].min() >= zlo and tr[:, 0].max() < zhi
    for key, sl, org2, thin in subs:
        sub = np.ascontiguousarray(blk[sl]).astype(np.uint8)
        if not sub.any():
            continue
        pts = (np.concatenate(G["A"][key]).astype(np.int64) - np.asarray(org2)) if G["A"].get(key) else np.zeros((0, 3), np.int64)
        trB, trC = recover(sub, pts, p, org2, thin, G["med"])
        out.append((key, [x for x in trB if inside(x)], [x for x in trC if inside(x)]))
    return out


def main():
    SC, VID, Z0, Z1, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    os.makedirs(OUT, exist_ok=True); t0 = time.time()
    base = f"https://dl.ash2txt.org/datasets/spiral_datasets/{SC}/{VID}/tracks/"
    ej = re.findall(r'href="([^"]*extract\.json)"', get(base).decode())[0]
    meta = json.loads(get(base + ej)); p = meta["parameters"]
    (y0, x0), (y1, x1) = meta["yx_bounds"]; zlo, zhi = Z0 - M, Z1 + M
    pub = published(SC, VID, zlo, zhi, os.path.join(os.environ.get("DATA", os.path.expanduser("~/data")), f"{SC}_{VID}/pubidx"))
    url = meta["array_url"].replace("s3://vesuvius-challenge-open-data", "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com")
    arr = RemoteArray(url.rsplit("/", 1)[0], url.rsplit("/", 1)[1])
    keys = all_keys(p, arr.shape[0], y0, y1, x0, x1); A = {}; bad = 0
    for fam, ko, t in pub:
        k = keys[ko].decode() if ko < len(keys) else key_of(fam, t, p, y0, x0)
        pre, v = k.split(":"); v = int(v)
        axn = {"h": 0, "vy": 1, "vx": 2}[pre]; thk = p["z_chunk_depth_h"] if pre == "h" else p["yx_chunk_thickness_v"]
        if not (t[:, axn].min() >= v and t[:, axn].max() < v + thk):
            bad += 1; k = key_of(fam, t, p, y0, x0)
        A.setdefault(k, []).append(t)
    print(f"key reconstruction: {bad} of {len(pub)} tracks inconsistent with their ordinal key (fell back to geometry)", flush=True)
    print(f"A: {len(pub)} published tracks in [{zlo},{zhi}) over {len(A)} keys ({time.time() - t0:.0f}s)", flush=True)
    hz = list(range(p["z_min"] + ((zlo - p["z_min"]) // p["z_chunk_stride_h"] + 1) * p["z_chunk_stride_h"], zhi - p["z_chunk_depth_h"], p["z_chunk_stride_h"]))
    vy = list(range(y0, y1, p["yx_stride_v"])); vx = list(range(x0, x1, p["yx_stride_v"]))
    cz, cy, cx = arr.chunks
    # median native diameter at published h-track vertices, first full chunk-row of ribbons inside the band
    r0 = zlo - zlo % cz + (cz if zlo % cz else 0)
    blk = arr.read(r0, r0 + cz, y0, y1 + 1, x0, x1 + 1, workers=48) > 0; md = []
    for z in [z for z in hz if r0 <= z and z + p["z_chunk_depth_h"] <= r0 + cz][:12]:
        dn = native_diam(blk[z - r0:z - r0 + p["z_chunk_depth_h"]], 0)
        for t in A.get(f"h:{z}", []):
            q = t[:, 1:].astype(np.int64) - np.array([y0, x0]); ok = np.all((q >= 0) & (q < np.array(dn.shape)), 1)
            md.extend(dn[q[ok, 0], q[ok, 1]].tolist())
    med_diam = float(np.median(md)) if md else 6.0; del blk
    print(f"median native diameter at published h-track vertices {med_diam:.2f} vx (n={len(md)}); C cuts > {FUSED * med_diam:.2f}", flush=True)
    tasks = [("h", max(r, zlo), min(r + cz, zhi)) for r in range(zlo - zlo % cz, zhi, cz)]
    tasks += [("vy", r, r + cy) for r in range(y0 - y0 % cy, y1, cy)] + [("vx", r, r + cx) for r in range(x0 - x0 % cx, x1, cx)]
    G.update(p=p, arr=arr, zlo=zlo, zhi=zhi, yx=((y0, x0), (y1, x1)), hz=hz, vy=vy, vx=vx, A=A, med=med_diam)
    B, C = {k: list(v) for k, v in A.items()}, {k: list(v) for k, v in A.items()}
    import multiprocessing as mp
    nb = nc = 0
    with mp.get_context("fork").Pool(int(os.environ.get("BAND_WORKERS", "10"))) as pool:
        for i, res in enumerate(pool.imap_unordered(task, tasks)):
            for key, trB, trC in res:
                B.setdefault(key, []).extend(trB); C.setdefault(key, []).extend(trC); nb += len(trB); nc += len(trC)
            print(f"  block {i + 1}/{len(tasks)}: +B {nb} +C {nc} ({time.time() - t0:.0f}s)", flush=True)
    # A* (noise-floor arm, pre-registered): A with a random 5% of its tracks removed (rng 7), same keys
    rs = np.random.default_rng(7); Astar = {}
    for k in sorted(A):
        keepm = rs.random(len(A[k])) >= 0.05; Astar[k] = [t for t, kk in zip(A[k], keepm) if kk]
    for name, arm in (("A", A), ("Astar", Astar), ("B", B), ("C", C)):
        d = os.path.join(OUT, name); os.makedirs(d, exist_ok=True)
        path = os.path.join(d, f"{SC}_{VID}_surface_m7_L0_th0.2.dbm")
        with dbm.open(path, "n") as db:
            for k, v in arm.items():
                db[k] = pickle.dumps([np.asarray(t, np.int32) for t in v])
        print(f"arm {name}: {sum(len(v) for v in arm.values())} tracks -> {path}", flush=True)
    json.dump(dict(scroll=SC, volid=VID, z=[Z0, Z1], margin=M, med_native_diam=med_diam, added_B=nb, added_C=nc,
                   n_A=sum(len(v) for v in A.values())), open(os.path.join(OUT, "band_summary.json"), "w"), indent=1)
    print(f"DONE {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
