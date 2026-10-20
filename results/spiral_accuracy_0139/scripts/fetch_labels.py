"""Fetch the public PHerc0139 inputs that are not tracks (CC BY-NC 4.0, no auth):
  - every w-labelled segment's tifxyz mesh on the fit volume 20250728140407 (9.362 um) -> data/labels/wNNN_<segid>/{x,y,z}.tif + meta.json
  - the same segments' meshes on 20260102150214 (2.399 um) for the native normal check (NATIVE list only)
  - the published umbilicus (2.399 um frame) and the catalog properties/transform of both volumes
Run from the repo root: python -I work/spiral_accuracy/scripts/fetch_labels.py"""
import gzip, json, os, re, urllib.request
from concurrent.futures import ThreadPoolExecutor
S3 = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
D = os.path.join(ROOT, "work/spiral_accuracy/data")
NATIVE = ["w025", "w031", "w037", "w043", "w049", "w055"]   # every 6th used wrap, fixed before any look


def get(u):
    for _ in range(6):
        try:
            return urllib.request.urlopen(u, timeout=300).read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
        except Exception:
            pass
    raise RuntimeError(u)


def lst(p):
    return re.findall(r"<Prefix>([^<]*)</Prefix>", get(f"{S3}/?list-type=2&delimiter=/&prefix={p}").decode())[1:]


jobs = []
for s in lst("PHerc0139/segments/"):
    sid = s.split("/")[-2]; m = re.search(r"-(w\d{3})_", sid)
    if not m:
        continue
    w = m.group(1)
    for mesh in lst(f"{s}mesh/"):
        frame = "9um" if "20250728140407" in mesh else ("2um" if "20260102150214" in mesh else None)
        if frame is None or (frame == "2um" and w not in NATIVE):
            continue
        out = os.path.join(D, f"labels_{frame}", f"{w}_{sid}")
        for f in ("meta.json", "x.tif", "y.tif", "z.tif"):
            jobs.append((f"{S3}/{mesh}{f}", os.path.join(out, f)))


def dl(j):
    u, p = j
    if os.path.exists(p):
        return
    b = get(u)
    if b is None:
        print("MISSING", u); return
    os.makedirs(os.path.dirname(p), exist_ok=True); open(p + ".part", "wb").write(b); os.replace(p + ".part", p)


with ThreadPoolExecutor(16) as ex:
    list(ex.map(dl, jobs))
um = re.findall(r"<Key>([^<]*)</Key>", get(f"{S3}/?list-type=2&prefix=PHerc0139/representations/umbilicus/").decode())
open(os.path.join(D, "umbilicus_published_2399.json"), "wb").write(get(f"{S3}/{um[-1]}"))
cat = json.loads(gzip.decompress(get(f"{S3}/metadata.json")))
props = {v: cat["samples"]["PHerc0139"]["volumes"][v]["properties"] for v in ("20250728140407", "20260102150214")}
json.dump(dict(umbilicus_key=um[-1], properties=props), open(os.path.join(D, "catalog_0139.json"), "w"), indent=1)
print("files", len(jobs), "umbilicus", um[-1])
