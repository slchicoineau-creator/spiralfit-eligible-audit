"""Reader for villa GridStore .grid files (normal-grids xy slices): list of (N,2) int32 polylines (x, y), grid units
= full-res voxels for input-level 0 grids. Format (GridStore.cpp save/load_mmap, v>=2): 13 big-endian u32 header
(magic, version, bounds x/y/w/h, cell, num_buckets, num_paths, buckets_offset, paths_offset, json_meta_offset,
json_meta_size); the paths region [paths_offset, json_meta_offset) holds seglists back to back:
start_x, start_y (u32 BE), n (u32 BE), then n int8 = (dx, dy) pairs."""
import numpy as np


def read_grid(path):
    b = open(path, "rb").read()
    if not b:
        return []
    h = np.frombuffer(b[:52], ">u4")
    version, npaths, poff, moff = int(h[1]), int(h[8]), int(h[10]), int(h[11])
    end = moff if version >= 2 else len(b)
    out, p = [], poff
    while p < end and len(out) < npaths:
        sx, sy, n = np.frombuffer(b[p:p + 12], ">u4").astype(np.int64); p += 12
        d = np.frombuffer(b[p:p + n], np.int8).astype(np.int32).reshape(-1, 2); p += int(n)
        pts = np.vstack([[sx, sy], np.array([sx, sy]) + np.cumsum(d, 0)]) if len(d) else np.array([[sx, sy]])
        out.append(pts.astype(np.int32))
    return out


def bounds(path):
    h = np.frombuffer(open(path, "rb").read(52), ">u4")
    return dict(version=int(h[1]), x=int(h[2]), y=int(h[3]), w=int(h[4]), h=int(h[5]), cell=int(h[6]), num_paths=int(h[8]))


if __name__ == "__main__":
    import sys
    P = read_grid(sys.argv[1]); allp = np.vstack(P)
    print(bounds(sys.argv[1]), "paths", len(P), "points", len(allp), "x", allp[:, 0].min(), allp[:, 0].max(), "y", allp[:, 1].min(), allp[:, 1].max(),
          "median path len", int(np.median([len(q) for q in P])))
