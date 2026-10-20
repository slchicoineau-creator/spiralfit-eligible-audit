"""Can blind ratings of random crops see fitter seed noise? (post-hoc re-analysis of the seed-replicated pilot)
1. Per (scroll, arm): the range of the 3 per-seed mean ratings, against a null in which the 24 crop ratings are
   shuffled across the 3 seeds (20,000 draws): is the seed-to-seed range larger than crop sampling alone gives?
2. Per-seed C - A on each scroll.
3. The same seed-1 fits rated in the pilot (12 crops) and again in the replicate (8 fresh crops).
4. Pilot C vs A, crop permutation test (12 vs 12, 50,000 draws, two-sided).
Writes results/rating_noise.json."""
import json, os, numpy as np
R = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")


def cells(key, ratings, fields):
    K = json.load(open(os.path.join(R, key)))["key"]; by = {}
    for a in json.load(open(os.path.join(R, ratings)))["answers"]:
        by.setdefault(a["id"], []).append(a["rating"])
    out = {}
    for iid, v in K.items():
        out.setdefault(tuple(v[f] for f in fields), []).append(float(np.mean(by[iid])))  # repeats averaged
    return {k: np.array(v) for k, v in out.items()}


rng = np.random.default_rng(0); res = {"seed_range_vs_crop_null": {}, "per_seed_C_minus_A": {}, "same_fit_rerated": {},
                                       "pilot_C_minus_A": {}}
S = cells("seed_rating_key.json", "trackfix_seed_ratings.json", ("scroll", "arm", "seed"))
P = cells("pilot_rating_key.json", "trackfix_pilot_ratings.json", ("scroll", "arm"))
for sc in ("PHerc0211", "PHerc0826"):
    for arm in ("A", "C"):
        L = [S[(sc, arm, s)] for s in ("s1", "s2", "s3")]; means = [x.mean() for x in L]; obs = max(means) - min(means)
        pool = np.concatenate(L); null = []
        for _ in range(20000):
            p = rng.permutation(pool); m = [p[i * 8:(i + 1) * 8].mean() for i in range(3)]; null.append(max(m) - min(m))
        null = np.array(null)
        res["seed_range_vs_crop_null"][f"{sc} {arm}"] = dict(
            seed_means=[round(m, 2) for m in means], range=round(obs, 2), crop_sd=round(float(pool.std(ddof=1)), 2),
            null_median=round(float(np.median(null)), 2), null_p95=round(float(np.quantile(null, .95)), 2),
            p=round(float((null >= obs - 1e-9).mean()), 3))
    res["per_seed_C_minus_A"][sc] = [round(S[(sc, "C", s)].mean() - S[(sc, "A", s)].mean(), 2) for s in ("s1", "s2", "s3")]
    for arm in ("A", "C"):
        res["same_fit_rerated"][f"{sc} {arm} seed1"] = dict(pilot=round(P[(sc, arm)].mean(), 2),
                                                           replicate=round(S[(sc, arm, "s1")].mean(), 2))
    a, c = P[(sc, "A")], P[(sc, "C")]; obs = c.mean() - a.mean(); pool = np.r_[a, c]; null = []
    for _ in range(50000):
        p = rng.permutation(pool); null.append(p[12:].mean() - p[:12].mean())
    res["pilot_C_minus_A"][sc] = dict(effect=round(obs, 2), p_two_sided=round(float((np.abs(null) >= abs(obs) - 1e-9).mean()), 3))
json.dump(res, open(os.path.join(R, "rating_noise.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
