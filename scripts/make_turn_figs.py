"""Figures for README sections 3-4, from results/fit_compare.json and (if present) results/tolerance_curve.json:
figs/turn_offsets_by_sector.png  modal turn offset per 45-degree sector x winding third, seed change vs umbilicus shift
figs/umbilicus_tolerance.png     turn agreement and outer sector range against the umbilicus shift, seed band shaded"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = json.load(open(os.path.join(REPO, "results", "fit_compare.json")))
panels = [("0826_8136 seed1 vs seed2", "seed 1 vs seed 2 (same inputs)"),
          ("0826_8136 umbilicus +250x", "umbilicus shifted +250 vx (2.3 mm) in x")]
fig, axes = plt.subplots(1, 2, figsize=(10, 3.2))
for ax, (key, title) in zip(axes, panels):
    sec = R[key]["sector_modal_offset"]; bands = ["inner", "middle", "outer"]; cols = [f"{45 * s:03d}" for s in range(8)]
    M = np.array([[sec[b].get(c, np.nan) for c in cols] for b in bands], float)
    im = ax.imshow(M, cmap="coolwarm", vmin=-7, vmax=7, aspect="auto")
    for i in range(3):
        for j in range(8):
            if np.isfinite(M[i, j]):
                ax.text(j, i, f"{int(M[i, j]):+d}", ha="center", va="center", fontsize=9)
    ax.set_xticks(range(8)); ax.set_xticklabels([f"{45 * s}°" for s in range(8)], fontsize=8)
    ax.set_yticks(range(3)); ax.set_yticklabels(["inner third", "middle third", "outer third"], fontsize=8)
    ax.set_title(title, fontsize=9); ax.set_xlabel("sector around the published umbilicus (0° = +x)", fontsize=8)
fig.suptitle("PHerc0826 [8136, 9136): turn offset vs the baseline fit, modal per sector (0 everywhere = same numbering)", fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(REPO, "figs", "turn_offsets_by_sector.png"), dpi=130)
print("wrote figs/turn_offsets_by_sector.png")

T = os.path.join(REPO, "results", "tolerance_curve.json")
if os.path.exists(T):
    t = json.load(open(T)); f = t["fits"]; s = t.get("summary", {})
    seeds = [f[k] for k in ("seed2", "seed3", "seed4") if k in f]
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(12, 3.4))
    for ax, get, lab in ((a1, lambda r: r["turn_share_at_mode"], "share of co-located surface at the modal turn offset"),
                         (a2, lambda r: r["outer_sector_range"], "outer-third sector range (turns; 0 = relabel)"),
                         (a3, lambda r: r["point_to_surface"]["median"], "median point-to-surface distance (vx)")):
        v = [get(r) for r in seeds]; ax.axhspan(min(v), max(v), color="#bbb", alpha=.5, label="seed 2-4 (no shift)")
        for d, m in (("x", "o-"), ("y", "s-")):
            xs = [k for k in (50, 100, 250, 500) if f"{d}{k}" in f]
            ax.plot([k * 9.362e-3 for k in xs], [get(f[f"{d}{k}"]) for k in xs], m, label=f"shift +{d}")
        if "xm250" in f:
            ax.plot([250 * 9.362e-3], [get(f["xm250"])], "k^", label="shift -x")
        ax.set_xscale("log"); ax.set_xlabel("umbilicus shift (mm)", fontsize=8); ax.set_title(lab, fontsize=8)
        ax.set_xticks([0.47, 0.94, 2.34, 4.68]); ax.set_xticklabels(["0.47", "0.94", "2.34", "4.68"], fontsize=8)
        ax.minorticks_off()
    if "share_floor" in s:
        a1.axhline(s["share_floor"], ls=":", c="k", lw=1); a2.axhline(s["range_limit"], ls=":", c="k", lw=1)
    a1.legend(fontsize=7)
    fig.suptitle(f"PHerc0826 [8136, 9136): how far can the umbilicus be off? Pre-registered tolerance: {s.get('headline', '?')} "
                 "(dotted = rule thresholds)", fontsize=9)
    plt.tight_layout(); plt.savefig(os.path.join(REPO, "figs", "umbilicus_tolerance.png"), dpi=130)
    print("wrote figs/umbilicus_tolerance.png")
