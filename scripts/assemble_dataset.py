"""Assemble a tracks-only spiral dataset root in WSL, per villa PR #1837's recipe:
spiral-scroll.json, umbilicus.json (published), lasagna_inputs/las_008_{nx,ny}.ome.zarr (only the chunk rows covering
the z band, plus metadata). Tracks are placed separately. Usage: python assemble_dataset.py SCROLL VOLID Z0 Z1 N_WIND"""
import json, os, re, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor
S3 = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
SC, VID, Z0, Z1, NW = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
D = os.path.join(os.environ.get("DATA", os.path.expanduser("~/data")), f"{SC}_{VID}"); os.makedirs(D, exist_ok=True)


def get(u):
    for a in range(6):
        try:
            return urllib.request.urlopen(u, timeout=120).read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
        except Exception:
            pass
    raise RuntimeError(u)


def lst(prefix):
    return get(f"{S3}/?list-type=2&delimiter=/&prefix={prefix}").decode()


# umbilicus (published, must be for this volume)
um = re.findall(r"<Key>([^<]*umbilicus[^<]*\.json)</Key>", lst(f"{SC}/representations/umbilicus/"))
um = [k for k in um if VID in k]
raw = json.loads(get(f"{S3}/{um[-1]}"))
json.dump(raw, open(f"{D}/umbilicus_published.json", "w"))
if "control_points" in raw:
    json.dump(raw, open(f"{D}/umbilicus.json", "w"))
else:
    print("umbilicus format keys:", list(raw)[:10])
# volume catalog properties
vols = re.findall(r"<Prefix>([^<]*)</Prefix>", lst(f"{SC}/volumes/"))
vol = [v for v in vols if VID in v and v.endswith(".zarr/")][0]
md = json.loads(get(f"{S3}/{vol}metadata.json"))
def find(o, k):
    if isinstance(o, dict):
        if k in o:
            return o[k]
        for v in o.values():
            r = find(v, k)
            if r is not None:
                return r
    return None
zdir, lh = find(md, "z_direction_is_top_to_bottom"), find(md, "left_handed_coordinates")
if zdir is None:                                   # bucket catalog (gzip JSON at the bucket root)
    import gzip
    cat = json.loads(gzip.decompress(get(f"{S3}/metadata.json")))
    pr = cat["samples"][SC]["volumes"][VID]["properties"]
    zdir, lh = pr.get("z_direction_is_top_to_bottom"), pr.get("left_handed_coordinates")
# lasagna
runs = sorted(re.findall(r"<Prefix>([^<]*lasagna-[^<]*/)</Prefix>", lst(f"{SC}/representations/predictions/lasagna/")))
run = [r for r in runs if VID in r][-1]
GROUP = "2"
zattrs = json.loads(get(f"{S3}/{run}{SC}_nx.ome.zarr/.zattrs"))
scale = [d for d in zattrs["multiscales"][0]["datasets"] if d["path"] == GROUP][0]["coordinateTransformations"][0]["scale"][0]
lscale = int(scale)
vox_um = float(re.findall(r"-(\d+\.\d+)um", vol)[0])
spec = dict(schema_version=1, name=SC, voxel_size_um=vox_um, normal_zarr_group=GROUP, lasagna_scale=lscale,
            paths=dict(tracks_dbm=f"tracks/{SC}_{VID}_surface_m7_L0_th0.2.dbm"))
if zdir is not None and lh is not None:
    spec.update(z_direction_is_top_to_bottom=zdir, left_handed_coordinates=lh)
json.dump(spec, open(f"{D}/spiral-scroll.json", "w"), indent=1)
jobs = []
for ch in ("nx", "ny"):
    src = f"{S3}/{run}{SC}_{ch}.ome.zarr"; dst = f"{D}/lasagna_inputs/las_008_{ch}.ome.zarr"
    os.makedirs(f"{dst}/{GROUP}", exist_ok=True)
    for m in (".zattrs", ".zgroup"):
        b = get(f"{src}/{m}")
        if b is not None:
            open(f"{dst}/{m}", "wb").write(b)
    za = json.loads(get(f"{src}/{GROUP}/.zarray")); open(f"{dst}/{GROUP}/.zarray", "w").write(json.dumps(za))
    cz, cy, cx = za["chunks"]; sz, sy, sx = za["shape"]
    zlo = max(Z0 // lscale // cz - 1, 0); zhi = min((Z1 // lscale) // cz + 2, (sz + cz - 1) // cz)
    for iz in range(zlo, zhi):
        for iy in range((sy + cy - 1) // cy):
            for ix in range((sx + cx - 1) // cx):
                jobs.append((f"{src}/{GROUP}/{iz}/{iy}/{ix}", f"{dst}/{GROUP}/{iz}/{iy}/{ix}"))


def fetch(j):
    u, f = j
    if os.path.exists(f):
        return 0
    b = get(u)
    if b is None:
        return 0
    os.makedirs(os.path.dirname(f), exist_ok=True); open(f, "wb").write(b); return 1


with ThreadPoolExecutor(48) as ex:
    n = sum(ex.map(fetch, jobs))
print(f"{SC}: umbilicus {um[-1]}; z_dir {zdir} left_handed {lh}; lasagna {run} group {GROUP} scale {lscale}; chunks fetched {n} of {len(jobs)}")
