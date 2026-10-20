"""Figures for the README. (1) census bar chart from results/census_*.json. (2) optional overlay of one published-track
ribbon on the extractor's giant component (needs the downloaded PHerc0125 .vctracks; pass --overlay TRACKS_DIR)."""
import glob, json, os, sys, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(REPO, "results", "census_PHerc*.json")))]
rows = [r for r in rows if "dropped_share" in r]
names = [r["scroll"].replace("PHerc", "") + (" (control)" if r["scroll"] == "PHerc0139" else "") for r in rows]
fig, ax = plt.subplots(figsize=(8, 3.6))
ax.bar(names, [100 * r["dropped_share"] for r in rows], color=["#999" if "control" in n else "#c0392b" for n in names])
ax.set_ylabel("% of m7 foreground in dropped components"); ax.set_ylim(0, 60)
ax.set_title("Published tracks: share of the surface prediction the size cap removes (3 ribbons per scroll)")
for i, r in enumerate(rows):
    ax.text(i, 100 * r["dropped_share"] + 1, f"{100 * r['dropped_share']:.0f}%", ha="center", fontsize=8)
plt.xticks(rotation=30, fontsize=8); plt.tight_layout(); plt.savefig(os.path.join(REPO, "figs", "census_dropped.png"), dpi=130)
print("wrote figs/census_dropped.png")
if "--overlay" in sys.argv:
    sys.path.insert(0, os.path.join(REPO, "lib")); os.environ.setdefault("SB_CACHE", os.path.join(REPO, "cache", "mask"))
    import cc3d
    from zio import RemoteArray
    D = sys.argv[sys.argv.index("--overlay") + 1]; B = "PHerc0125_20250821151825_surface_m7_L0_th0.2"; V = f"{D}/{B}.dbm.vctracks"
    meta = json.load(open(f"{D}/{B}.extract.json")); p = meta["parameters"]; ds = p["downsample_factor"]; Z = 7748
    (y0, x0), (y1, x1) = meta["yx_bounds"]
    url = meta["array_url"].replace("s3://vesuvius-challenge-open-data", "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com")
    arr = RemoteArray(url.rsplit("/", 1)[0], url.rsplit("/", 1)[1])
    pred = (arr.read(Z, Z + p["z_chunk_depth_h"], y0, y1 + 1, x0, x1 + 1, workers=32) > 0).astype(np.uint8)
    cc = cc3d.connected_components(pred, connectivity=6); c = tuple(s - s % ds for s in cc.shape)
    lab = cc[:c[0], :c[1], :c[2]].reshape(c[0] // ds, ds, c[1] // ds, ds, c[2] // ds, ds).max((1, 3, 5))[0]
    ids, cnt = np.unique(lab[lab > 0], return_counts=True); giant = np.isin(lab, ids[cnt > p["max_area_threshold"] // ds ** 3])
    off = np.fromfile(f"{V}/offsets.i64", "<i8"); fam = np.fromfile(f"{V}/family_codes.i8", "<i1"); zb = np.fromfile(f"{V}/z_bounds.i32", "<i4").reshape(-1, 2)
    coords = np.memmap(f"{V}/coordinates.i32", "<i4", "r").reshape(-1, 3)
    sel = np.nonzero((fam == 0) & (zb[:, 0] >= Z) & (zb[:, 1] < Z + p["z_chunk_depth_h"]))[0]
    pts = np.concatenate([np.asarray(coords[off[i]:off[i + 1]]) for i in sel])
    img = np.zeros(lab.shape + (3,), np.uint8); img[lab > 0] = (90, 90, 90); img[giant] = (200, 50, 40)
    py = np.clip((pts[:, 1] - y0) // ds, 0, lab.shape[0] - 1); px = np.clip((pts[:, 2] - x0) // ds, 0, lab.shape[1] - 1); img[py, px] = (60, 220, 255)
    plt.figure(figsize=(7, 7)); plt.imshow(img); plt.axis("off")
    plt.title("PHerc0125 ribbon z=7748: m7 foreground (grey), dropped giant component (red), published tracks (cyan)", fontsize=8)
    plt.tight_layout(); plt.savefig(os.path.join(REPO, "figs", "published_tracks_vs_dropped.png"), dpi=130); print("wrote overlay")
