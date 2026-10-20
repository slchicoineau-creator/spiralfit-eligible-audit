"""Automatic umbilicus estimate from published lasagna normals (pre-registered in the work log). At each control slice z,
take the lasagna nx/ny store (group 2, 4x down; decoded (u8-128)/127), keep pixels with in-plane normal length >= 0.5
(sheets roughly parallel to z), and find the point p minimising the robust (Huber, IRLS) sum of squared distances to the
lines through each pixel along its in-plane normal (sheet normals point roughly radially, so they converge on the core).
Validation mode compares to a published umbilicus. Usage:
  python umb_estimate.py SCROLL VOLID Z0 Z1 [--validate]   -> umb_est_<SCROLL>.json (control points every 250 slices)"""
import json, os, re, sys, urllib.request, time, numpy as np
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib")); os.environ.setdefault("SB_CACHE", os.path.join(REPO, "cache", "lasagna"))
from zio import RemoteArray  # noqa
S3 = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"


def get(u):
    for a in range(8):
        try:
            return urllib.request.urlopen(u, timeout=60).read()
        except urllib.error.HTTPError:
            raise
        except Exception:
            time.sleep(2 ** a)
    raise RuntimeError(u)


def estimate(nx, ny, scale, huber=40.0):
    """nx, ny: 2D decoded in-plane normal components at 1/scale resolution. Returns (x, y) in full-res voxels."""
    L = np.hypot(nx, ny); m = L >= 0.5
    if m.sum() < 1000:
        return None
    ys, xs = np.nonzero(m); d = np.stack([nx[m], ny[m]], 1) / L[m][:, None]           # unit in-plane normal (x, y)
    P = np.stack([xs, ys], 1).astype(float) * scale
    w = np.ones(len(P))
    p = P.mean(0)
    for _ in range(30):
        # distance from p to the line through P along d = |(p-P) x d|
        r = np.abs((p[0] - P[:, 0]) * d[:, 1] - (p[1] - P[:, 1]) * d[:, 0])
        w = np.where(r <= huber, 1.0, huber / np.maximum(r, 1e-9))
        A = np.zeros((2, 2)); b = np.zeros(2)
        Pm = np.eye(2)[None] - d[:, :, None] * d[:, None, :]                            # projector orthogonal to d
        A = (w[:, None, None] * Pm).sum(0); b = (w[:, None, None] * Pm @ P[:, :, None])[..., 0].sum(0)
        pn = np.linalg.solve(A, b)
        if np.linalg.norm(pn - p) < 0.1:
            p = pn; break
        p = pn
    return float(p[0]), float(p[1])


def main():
    SC, VID, Z0, Z1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    lst = get(f"{S3}/?list-type=2&delimiter=/&prefix={SC}/representations/predictions/lasagna/").decode()
    run = [r for r in re.findall(r"<Prefix>([^<]*lasagna-[^<]*/)</Prefix>", lst) if VID in r][-1]
    za = json.loads(get(f"{S3}/{run}{SC}_nx.ome.zarr/.zattrs"))
    scale = int([d for d in za["multiscales"][0]["datasets"] if d["path"] == "2"][0]["coordinateTransformations"][0]["scale"][0])
    NX, NY = RemoteArray(f"{S3}/{run}{SC}_nx.ome.zarr", 2), RemoteArray(f"{S3}/{run}{SC}_ny.ome.zarr", 2)
    zs = list(range(Z0 - 500, Z1 + 501, 250)); pts = []
    for z in zs:
        zi = z // scale
        if not (0 <= zi < NX.shape[0]):
            continue
        nx = (NX.read(zi, zi + 1, 0, NX.shape[1], 0, NX.shape[2], workers=32)[0].astype(np.float32) - 128) / 127
        ny = (NY.read(zi, zi + 1, 0, NY.shape[1], 0, NY.shape[2], workers=32)[0].astype(np.float32) - 128) / 127
        e = estimate(nx, ny, scale)
        if e:
            pts.append(dict(x=round(e[0], 1), y=round(e[1], 1), z=z))
            print(f"{SC} z={z}: umbilicus estimate x={e[0]:.0f} y={e[1]:.0f}", flush=True)
    out = dict(control_points=pts, metadata=dict(volume=VID, method="lasagna in-plane normal convergence (IRLS Huber 40 vx), group 2"))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, f"umb_est_{SC}.json"), "w"), indent=1)
    if "--validate" in sys.argv:
        um = [k for k in re.findall(r"<Key>([^<]*)</Key>", get(f"{S3}/?list-type=2&prefix={SC}/representations/umbilicus/").decode()) if VID in k][-1]
        cp = json.loads(get(f"{S3}/{um}"))["control_points"]
        zz = np.array([c["z"] for c in cp]); err = []
        for q in pts:
            px, py = np.interp(q["z"], zz, [c["x"] for c in cp]), np.interp(q["z"], zz, [c["y"] for c in cp])
            err.append(float(np.hypot(q["x"] - px, q["y"] - py)))
        print(f"VALIDATION {SC}: errors vs published (vx) {[round(e) for e in err]} median {np.median(err):.0f}")


if __name__ == "__main__":
    main()
