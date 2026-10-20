"""W2 check 1 (advisor): do the PUBLISHED tracks omit the giant components that the census reproduction says the
extractor drops? PHerc0125, horizontal ribbon at z = 7748 (a census ribbon). Reads the published .vctracks arrays raw
(offsets.i64, coordinates.i32 zyx, family_codes.i8 [0 = h], z_bounds.i32), recomputes the ribbon's pooled CC labels
exactly as the census, and reports track coverage of giant-component vs kept-component footprints. Runs in WSL."""
import json, os, sys, numpy as np, cc3d
from scipy.ndimage import binary_dilation
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
os.environ.setdefault("SB_CACHE", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cache", "mask"))
from zio import RemoteArray  # noqa
D = os.path.join(os.environ.get("DATA", os.path.expanduser("~/data")), f"PHerc0125_20250821151825/tracks")
B = "PHerc0125_20250821151825_surface_m7_L0_th0.2"; V = f"{D}/{B}.dbm.vctracks"
meta = json.load(open(f"{D}/{B}.extract.json")); p = meta["parameters"]; ds = p["downsample_factor"]
Z = int(sys.argv[1]) if len(sys.argv) > 1 else 7748
(y0, x0), (y1, x1) = meta["yx_bounds"]
url = meta["array_url"].replace("s3://vesuvius-challenge-open-data", "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com")
arr = RemoteArray(url.rsplit("/", 1)[0], url.rsplit("/", 1)[1])
pred = (arr.read(Z, Z + p["z_chunk_depth_h"], y0, y1 + 1, x0, x1 + 1, workers=32) > 0).astype(np.uint8)
cc = cc3d.connected_components(pred, connectivity=6)
c = tuple(s - s % ds for s in cc.shape); lab = cc[:c[0], :c[1], :c[2]].reshape(c[0] // ds, ds, c[1] // ds, ds, c[2] // ds, ds).max((1, 3, 5))[0]
ids, cnt = np.unique(lab[lab > 0], return_counts=True)
big = ids[cnt > p["max_area_threshold"] // ds ** 3]; small = ids[cnt < max(1, p["dust_threshold"] // ds ** 3)]
giant = np.isin(lab, big); kept = (lab > 0) & ~giant & ~np.isin(lab, small)
off = np.fromfile(f"{V}/offsets.i64", "<i8"); fam = np.fromfile(f"{V}/family_codes.i8", "<i1")
zb = np.fromfile(f"{V}/z_bounds.i32", "<i4").reshape(-1, 2)
coords = np.memmap(f"{V}/coordinates.i32", "<i4", "r").reshape(-1, 3)
sel = np.nonzero((fam == 0) & (zb[:, 0] >= Z) & (zb[:, 1] < Z + p["z_chunk_depth_h"]))[0]
pts = np.concatenate([np.asarray(coords[off[i]:off[i + 1]]) for i in sel]) if len(sel) else np.zeros((0, 3), int)
print(f"z={Z}: published h-tracks {len(sel)} with {len(pts)} points; giant comps {len(big)}")
py = np.clip((pts[:, 1] - y0) // ds, 0, lab.shape[0] - 1); px = np.clip((pts[:, 2] - x0) // ds, 0, lab.shape[1] - 1)
near = binary_dilation(giant, iterations=1)
in_g = giant[py, px].mean() if len(pts) else 0; in_k = kept[py, px].mean() if len(pts) else 0
cov = np.zeros_like(giant); cov[py, px] = True; cov = binary_dilation(cov, iterations=2)
print(f"share of published track points on giant footprint {in_g:.3f}, on kept components {in_k:.3f}")
print(f"giant footprint within 2 px of a published track: {cov[giant].mean():.3f}; kept footprint: {cov[kept].mean():.3f}")
print(f"giant = {giant.sum()} px ({giant.sum() / max((lab > 0).sum(), 1):.3f} of pooled fg); kept = {kept.sum()} px")
