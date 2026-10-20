# SPIRAL-ACCURACY-0139: result (final, human-pass version, 2026-10-14)

*Public-safe copy of `RESULT.md` (2026-10-14), for package (a)'s `results/spiral_accuracy_0139/`. It is identical to the
internal record except for this note and the data-use line under "Pre-registration and amendments" (reworded without
private names). Every number and the judge's allowed and forbidden wording are unchanged.*

This file replaces the 2026-10-13 interim `RESULT.md`. It applies conditions C1-C5 of the independent judge review
`reviews/2026-10-14_partner_judge_spiral_accuracy.md` (commit ba22d215). Every number comes from a file listed in the last
section. No interpretation is added beyond the judge's. Owner of package (a): Sebastien; he decides (prereg:6-7, 200).

## Pre-registration and amendments

| Document | Commit | Pushed | Content |
|---|---|---|---|
| `prereg/2026-10-13_partner_spiral_accuracy.md` | 25fc3143 | 10-13 10:31, before any fit | question, truth set, metric, rules |
| `..._amend1.md` | 4015ef8f | 10-13 11:29, before any 30k-step fit finished | Stage 0 record; secondary E_clean; descriptive outputs |
| `..._amend2.md` | d3533566 | 10-13 12:17, before any 30k-step fit finished | compute layout only: fits on BOTH R9700 cards; base_s1 resumed; arm B on CPU |

- **Both R9700 cards were used** (amendment 2). Card 0 ran base_s1 (resumed), x50_s1, y25_s1, base_s3, y25_s2, x25_s1. Card 1 ran
  y50_s1, base_s2, y10_s1, y50_s2, y10_s2 (`logs/driver.log`). All fits are ROCm fits and are compared only with each other.
- **base_s1 was resumed** once: paused on purpose at 12:05, resumed from iteration 8,000 (amendment 2). The judge dropped it as a
  check: E_max becomes 0.55727, and **no verdict changes at any margin**.
- No human mark, answer key, rating, private scan or team-made mesh was used. Truth = the organizers' published PHerc0139 verified wraps.

## What it measures

- Package (a) shows that the tracks-only spiral fit changes its turn numbering when the umbilicus moves 0.09-0.47 mm (fit vs fit).
- Here the same rigid shifts (10 / 25 / 50 vx = 0.094 / 0.234 / 0.468 mm, +x and +y) are scored against truth: PHerc0139, band z
  [5706, 6706), 34 labelled wraps w023-w029, w031-w044, w047-w059 (w030, w045, w046 excluded as registered).
- Primary metric **E** = share of matched labelled area (truth vertices within 6 vx of a fitted winding) that the fit gives a
  different turn number from the fit's own bulk numbering (one global offset allowed). Lower is better.
- Recipe: package recipe, 30,000 steps, N = 142 windings (`out/wind_count.json`).
- Stage 0 (a-d) and the smoke fit passed (amendment 1, section A0; `out/s0_labels.json`, `out/s0_normals.json`).

## Verdict per arm, as registered (primary E, margin 0.01)

Base fits: E 0.55834 (s1), 0.55727 (s2), 0.55716 (s3). E_max 0.55834, E_min 0.55716, spread 0.00118 (0.12 pt).
VOID: no. Seed noise too large (spread > 0.10): no. Unscorable fits (M2 < 0.2): none. Source: `out/verdict.json`.

| Arm | Shift | Fits (E) | Headline: median E(arm) - median E(base) | **m = 0.01 (registered)** | m = 0.005 | m = 0.02 |
|---|---|---|---|---|---|---|
| y10 | 0.094 mm | 0.55721, 0.55743 | +0.005 pt (`verdict.json` rounds to 0.0) | **TOLERATED** | TOLERATED | TOLERATED |
| y25 | 0.234 mm | 0.55840, 0.55912 | +0.15 pt | **TOLERATED** | TOLERATED | TOLERATED |
| y50 | 0.468 mm | 0.56633, 0.56603 | +0.89 pt | **TOLERATED** | **ACCURACY COST** | TOLERATED |
| x10 | 0.094 mm | none (0 of 2 fits) | - | **NOT RUN (0 of 2 fits; cap)** | NOT RUN | NOT RUN |
| x25 | 0.234 mm | 1 of 2: x25_s1 0.55683 | -0.04 pt (description only) | **NOT RUN (1 of 2)** | NOT RUN | NOT RUN |
| x50 | 0.468 mm | 1 of 2: x50_s1 0.55823 | +0.10 pt (description only) | **NOT RUN (1 of 2)** | NOT RUN | NOT RUN |

- **x10 has no fit.** It is missing from `out/RESULT_auto.md` and from `out/verdict.json` (those files are job-owned and were not
  edited). It is listed here as the prereg requires (prereg:162-163, 195).
- **The x single-fit E values (x25_s1, x50_s1) are description only.** They are not a verdict.
- No arm is GAIN (no arm has every fit below E_min - m).
- **Secondary E_clean** (amendment 1; 30 wraps, w034, w035, w041, w042 dropped by a label-only rule): y10 -0.03, y25 +0.04, y50
  +0.44 pt. All three are TOLERATED at 0.005, 0.01 and 0.02. x25, x50: NOT RUN (1 of 2). Primary and secondary agree at the
  registered margin.

### y50: report these two facts next to "TOLERATED" (condition C5)

- The y50 shift raised E by **+0.89 pt**. Both y50 seeds sit 0.77-0.80 pt above E_max. The full base seed range is 0.12 pt.
- At the **0.005 sensitivity margin y50 flips to ACCURACY COST** (both fits E > E_max + 0.005).
- Judge: H1's wording ("beyond seed noise", prereg:19-20) is met descriptively for y50. It is not met by the registered rule.
- Judge: about half of the y50 rise sits on the 4 label-inconsistent wraps (E_clean +0.44 pt). The main one is w034, whose modal e
  moves from 0 to -1 (`out/accuracy.json`, `modal_e_per_wrap`).

## Per-fit numbers

| Fit | Steps | E (primary) | E_clean | M2 on-sheet 4 / 6 / 8 vx | M1 right-turn share | c* | breaks |
|---|---|---|---|---|---|---|---|
| base_s1 (resumed) | 30,000 | 0.55834 | 0.52456 | 0.7133 / 0.8186 / 0.8759 | 0.3616 | 26 | 12 |
| base_s2 | 30,000 | 0.55727 | 0.52323 | 0.7105 / 0.8166 / 0.8746 | 0.3615 | 26 | 12 |
| base_s3 | 30,000 | 0.55716 | 0.52254 | 0.7136 / 0.8195 / 0.8780 | 0.3629 | 26 | 13 |
| y10_s1 | 30,000 | 0.55721 | 0.52242 | 0.7138 / 0.8199 / 0.8769 | 0.3631 | 26 | 12 |
| y10_s2 | 30,000 | 0.55743 | 0.52338 | 0.7126 / 0.8183 / 0.8759 | 0.3622 | 26 | 12 |
| y25_s1 | 30,000 | 0.55840 | 0.52331 | 0.7166 / 0.8221 / 0.8800 | 0.3630 | 26 | 13 |
| y25_s2 | 30,000 | 0.55912 | 0.52390 | 0.7197 / 0.8237 / 0.8799 | 0.3632 | 26 | 11 |
| y50_s1 | 30,000 | 0.56633 | 0.52783 | 0.7109 / 0.8151 / 0.8718 | 0.3535 | 26 | 13 |
| y50_s2 | 30,000 | 0.56603 | 0.52751 | 0.7108 / 0.8146 / 0.8714 | 0.3535 | 26 | 13 |
| x25_s1 (arm not run) | 30,000 | 0.55683 | 0.52192 | 0.7093 / 0.8166 / 0.8748 | 0.3619 | 26 | 14 |
| x50_s1 (arm not run) | 30,000 | 0.55823 | 0.52381 | 0.7115 / 0.8166 / 0.8737 | 0.3608 | 26 | 12 |
| smoke_base_s1 | **4,000 steps, sanity only, never in a verdict** | 0.75759 | 0.74567 | 0.6516 / 0.8132 / 0.8963 | 0.1971 | 25 | 11 |

Source: `out/RESULT_auto.md`, `out/accuracy.json`. The fractional residual is 0.0008 in every fit. "breaks" = number of
consecutive labelled-wrap pairs whose modal e differs.

## GPU accounting (cap, condition C2 context)

- `out/gpu_hours.json`: 13 entries, **28.96 GPU-h** (cap 30). This includes the smoke (0.32) and the sigterm segment of base_s1 (0.65).
- Rule (prereg section 8; amendment 2): no new fit starts when logged hours plus 2.4 h per running fit reach 27.
- x25_s1 started at 00:33 (23.84 logged + 2.4 = 26.24 < 27). At 01:07 (26.40 logged + 2.4 for x25_s1 = 28.80 >= 27) x10_s1, x50_s2,
  x25_s2 and x10_s2 did not start (`logs/driver.log`). The start order follows the frozen priority list.
- Arm B ran on CPU and used no GPU-h.

## E limitation (judge finding 7; condition C5)

**E is near-saturated and blind to "wrong to another wrong" renumbering.**

- The base E is 0.557 (3 seeds: 0.55716-0.55834). The unshifted fit numbers only 44 % of the matched labelled area consistently
  with the labels, while a fitted winding lies within 6 vx (56 um) of 82 % of it (M2 at 6 vx, 0.8166-0.8195).
- Its modal offset steps from -4 (w023) to +1 (w031) to -2 (w036-w044) to 0 (w047-w059): 12-13 breaks over 33 wrap steps
  (`out/accuracy.json`, `modal_e_per_wrap`, `breaks`).
- E counts only right and wrong. It does not count a change from one wrong number to another wrong number.
- Judge's per-vertex concordance on truth vertices matched by both fits. These numbers are the judge's own recomputation (script in
  the judge review appendix; the run copy was not committed). It reproduces all 11 main E values exactly.

| Pair type | e differs | right -> wrong | wrong -> right | wrong -> other wrong | wraps whose modal e changes |
|---|---|---|---|---|---|
| base vs base (3 pairs) | 5.8-6.4 % | 1.2-1.5 % | 1.2-1.5 % | 2.9-3.8 % | 0-2 |
| y10 vs bases (6) | 5.5-6.9 % | 1.4-1.6 % | 1.4-1.6 % | 2.7-3.8 % | 0-2 |
| y25 vs bases (6) | 5.9-8.4 % | 1.5-1.6 % | 1.5-1.6 % | 2.8-5.3 % | 0-4 |
| y50 vs bases (6) | 8.1-10.9 % | 2.0-2.1 % | 1.8-1.9 % | 4.3-7.0 % | 2-5 |
| x25 vs bases (3) | 6.2-6.4 % | 1.5-1.6 % | 1.5-1.6 % | 3.0-3.4 % | 0-1 |
| x50 vs bases (3) | 5.8-7.4 % | 1.5-1.6 % | 1.5-1.6 % | 2.8-4.1 % | 0-2 |

- Judge's reading: inside the labelled footprint, y50 changes about 2-5 points more of the turn numbers than seed noise does. Most
  of the extra changes go from one wrong number to another, which E does not count. Right -> wrong and wrong -> right nearly cancel.
- Against base_s1, y50's modal e changes on w025-w028 (-3 to -4; already wrong), w038 (-2 to -3; already wrong) and w034 (0 to -1;
  the only right-to-wrong flip). Checked against `out/accuracy.json`.
- "TOLERATED" here means "the share numbered consistently with the labels did not drop by more than 1 point". It does **not**
  mean "the numbering did not change".
- Label defects add a constant floor to every fit, not a bias (judge finding 6): w035 lies inside w034 (median gap ratio -0.31);
  w041 lies on w042's sheet over 64 % of the band (ratio 0.06); a sheet is probably unlabelled between w040 and w041 (ratio 2.04);
  adjacent wraps lie within 4 vx of each other on 1-24 % of their area (`out/diag_label_gaps.json`, `out/diag_label_overlap.json`).

## H2, reported twice (condition C4; descriptive only)

H2 asked: where the package's fit-vs-fit rule flags a shift, does the accuracy rule also flag it?
The package rule: a fit pair is "not tolerated" if its same-turn share at the mode is below a floor, or its outer sector range is > 1.
Base-pair shares (`out/h2_fit_compare.json`): 0.873, 0.893, 0.873.

**Disclosure.** Prereg:168-170 calls H2 "the package's rule (floor = smallest base-pair cross share at the mode ...)" and
`scripts/summarize.py:60,69` implements exactly that. The package's actual rule uses floor = lowest seed-pair share **minus 0.10**
(`followup_rescore.py:6` on PR #4, branch `judge/spiralfit-tolerance` @ 069d9547; its README:157-158). So the registered H2 is
stricter than the package rule. Also, `scripts/h2.sh:9` compares each shifted fit with **all three** bases. The prereg registered
base-s1 only. The all-bases set is a superset, and the reading is the same.

| Arm | (1) As registered: floor 0.873 (margin-free); pairs flagged vs all 3 bases | (2) Package's published rule: floor 0.773 (0.873 - 0.10); pairs flagged vs all 3 bases | (2) vs base_s1 only |
|---|---|---|---|
| y10 | 6/6 | **0/6** | 0/2 |
| y25 | 6/6 | **6/6** (5/6 at margin 0.15) | 2/2 |
| y50 | 6/6 | **6/6** | 2/2 |
| x25 | 3/3 (one fit) | **0/3** (one fit) | 0/1 |
| x50 | 3/3 (one fit) | **3/3** (one fit; range leg 2-3) | 1/1 |

- Column (1) is `out/verdict.json` `H2.not_tolerated_pairs`. Column (2) is re-scored from `out/h2_fit_compare.json` with the package's
  rule (judge finding 8); recomputed for this file and it matches the judge.
- **H2 answer (judge): no.** On 0139 the package's fit-vs-fit rule brackets +y between 0.09 and 0.23 mm, as on 0211, while the
  accuracy rule tolerates all three +y shifts. The accuracy rule does not flag where the fit-vs-fit rule does.
- Where the disagreement goes (judge): whole-fit cross share (all windings w010-w141) is y50 0.637-0.675 and x50 0.724-0.734,
  against base pairs 0.873-0.893 (`out/h2_fit_compare.json`). Inside the labelled footprint, x50 is at or near seed level (5.8-7.4 %
  vs 5.8-6.4 %) and y50 is only 2-5 points above it (table above). Most of the fit-vs-fit disagreement falls on windings the labels
  do not cover. The labels map to fit windings about w045-w086 (c* 26, e -4..+1); the package's range leg uses w098-w141.
- Both H0 clauses (prereg:21-22) hold in part (judge).
- The package-rule bracket on 0139 is **descriptive**, not a registered result, and any use must name the margin change.

## Arm B: spiral fit + the organizers' tracer as refiner (CPU, no human input)

Result: **NO CLEAR CHANGE** (`out/armb_result.json`; judge finding 9 confirms it).

- 1,015 traces seeded on base_s1 windings (8 angles, z 6206). 425 reach step 24: reach share 0.4187, above the 0.05 floor.
- Versus the fit's own vertices in the same footprint (trace minus fit footprint):

| Measure | 4 vx | 6 vx | 8 vx |
|---|---|---|---|
| On-sheet share, traces | 0.383 | 0.406 | 0.4207 |
| On-sheet share, fit footprint | 0.296 | 0.3592 | 0.4003 |
| On-sheet difference (rule: >= 0.05 at all three) | +0.087 | +0.047 | +0.020 |
| Wrong-turn share, traces / fit footprint | 0.6583 / 0.6083 | 0.6624 / 0.6266 | 0.6663 / 0.6388 |

- Wrong-turn share at 6 vx: +0.0358 (worse; the HELPS rule also needs <= +0.01). The on-sheet difference is not >= 0.05 at all three
  tolerances, so HELPS is not met. HURTS (lower by >= 0.05 at all three) is not met either.

## What may and must not be said (judge, verbatim)

Package (a) may use the **+y** result only, after Sebastien's decision, with the wording below. The text in this section is copied
unchanged from the judge review (lines 215-261).

**Allowed wording** (README, a new subsection after §4.3; keep the numbers exactly):

> **Accuracy check on a third scroll (pre-registered, SPIRAL-ACCURACY-0139).** On PHerc0139 (9.362 µm, one 1,000-slice
> band), the same tracks-only recipe was scored against 34 of the organizers' published verified wraps (w023–w059 without
> w030, w045 and w046). With the published umbilicus, the fit put a winding within 6 vx (56 µm) of 82% of the labelled area,
> but numbered only 44% of that area consistently with the labels (one global offset allowed; wrong-turn share
> 0.557–0.558 over 3 seeds; the numbering slips against the labels at 12–13 of 33 wrap steps). Rigid +y umbilicus shifts of
> 0.09 / 0.23 / 0.47 mm changed the wrong-turn share by +0.0 / +0.15 / +0.89 percentage points (2 seeds each; the 3
> unshifted seeds span 0.12 points): all inside the pre-registered 1-point margin, so the shifts are *tolerated* by that rule,
> although the 0.47 mm shift would count as an accuracy cost at a 0.5-point margin. The shifted fits still number the
> windings differently from the unshifted ones (on this scroll our fit-vs-fit rule again flags +y at 0.23 mm but not at
> 0.09 mm, descriptive; the PHerc0139 pre-registration had named a stricter floor, under which every shift is
> flagged); most of that renumbering falls outside the labelled wraps or moves a turn from one wrong number to another,
> which the wrong-turn share does not count. +x shifts were not decided (one fit per shift; compute cap). Fits ran on our
> ROCm port of the fitter and are compared only with each other; the published m7 tracks may have seen PHerc0139 labels in
> training; two label pairs (w034/w035, w041/w042) are inconsistent in this band (without them the 0.47 mm effect is +0.44 points).

Short form (SUBMIT form or Discord; one sentence):

> On a third scroll (PHerc0139, 34 published verified wraps, pre-registered), +y umbilicus shifts up to 0.47 mm raised the
> share of the labelled area given an inconsistent turn number by at most 0.9 points (inside the registered 1-point margin),
> but the unshifted fit already numbered only 44% of that area consistently, and the shifted fits still numbered the
> windings differently, mostly where this share cannot see it.

**Required edits elsewhere if the paragraph goes in:**
- README:13 scope (a): "no result here says which fit is on the right sheet; ... not an accuracy one" → add "except the
  PHerc0139 accuracy check (§4.x)".
- README:382: "These fits' placement on the sheets was not tested here" → "... not tested on 0826/0211; see §4.x for
  PHerc0139".

**Forbidden:**
- F1. Any +x accuracy claim: "both directions", "x tolerated", or any x headline. x25 and x50 have 1 of 2 fits; x10 has none.
- F2. "The umbilicus does not matter for accuracy", "0.47 mm is accurate enough", "accuracy tolerance >= 0.47 mm", or "the
  shifts did not change the numbering on the labelled wraps".
- F3. That 0139 confirms, overturns or supersedes the "< 0.47 mm" (0826) or "0.09-0.23 mm" (0211) results. Different metric,
  region, scroll and device; see above.
- F4. Pooling the 0139 ROCm fits with the 0826/0211 fits, or comparing absolute E across devices. The CUDA-vs-ROCm check is
  NOT MET under strict pooling (prereg:77-79).
- F5. Generalising to unseen scrolls or to other recipes (patches, annotations, winding model). Also: m7 leak flag, one band,
  one scroll.
- F6. Any claim about the accuracy of the outer windings (unlabelled), or significance language (evidence unit = fits,
  3 vs 2).
- F7. E_clean or the smoke fit as a headline. The H2 package-rule bracket on 0139 as a *registered* result: it is
  descriptive and must name the margin change.
- F8. Arm B as an improvement. Allowed only: "seeding the organizers' sheet tracer on the fit gave NO CLEAR CHANGE by the
  registered rule (on-sheet +8.7 / +4.7 / +2.0 points at 4 / 6 / 8 vx, below the 5-point bar at 6 and 8 vx; wrong-turn share
  +3.6 points)".

## Judge conditions applied

| Condition | Where in this file |
|---|---|
| C1 x10 listed "NOT RUN (0 of 2 fits; cap)"; x25, x50 "NOT RUN (1 of 2)" with single-fit E as description | verdict table |
| C2 both cards; all 3 amendments; base_s1 resumed, dropping it changes no verdict | "Pre-registration and amendments" |
| C3 smoke row "4,000 steps, sanity only, never in a verdict" | per-fit table |
| C4 H2 twice (registered floor 0.873; package rule floor - 0.10); misdescription and all-bases superset disclosed | "H2, reported twice" |
| C5 E limitation; y50 +0.89 pt and 0.005-margin flip next to TOLERATED | "y50" and "E limitation" |

## The two missing x fits (judge)

- Not allowed as registered: prereg section 8 is a hard cap of 30 GPU-h with no new fit after 27, and 28.96 is logged. The fits
  need amendment 3, pushed before launch.
- Low value: x25_s1 (0.55683) and x50_s1 (0.55823) already sit at or below E_max, so the completion can return only TOLERATED or
  UNDECIDED, never ACCURACY COST. E's blindness (finding 7) also applies.
- Judge's option: run them only as idle-GPU filler (about 2.6 h on both cards), labelled a post-result completion. x10 stays NOT RUN.
  Package (a) does not need them.

## Limits (judge finding 5; forbidden claim F5)

- No truth label reaches a fit: the dataset root holds only `spiral-scroll.json`, `umbilicus.json`, `lasagna_inputs/` and `tracks/`.
  The shift is applied only to the fit's copy of the umbilicus. Scoring uses the unshifted mapped umbilicus.
- Indirect channels, equal across arms: the m7 tracks may have seen 0139 labels (`04_Work_Log.md:3583`); the lasagna normals
  (resampled from the 2.399 um scan, validated in Stage 0) and the published umbilicus have unknown provenance. They do not bias a
  shift-vs-base comparison. They limit transfer to other scrolls.
- One scroll, one 1,000-slice band, ROCm only. The CUDA-vs-ROCm check is NOT MET under strict pooling (prereg:77-79).
- Evidence unit = fits (3 base vs 2 per arm). No significance is claimed.

## Files

- Prereg and amendments: `prereg/2026-10-13_partner_spiral_accuracy.md`, `_amend1.md`, `_amend2.md`.
- Judge review: `reviews/2026-10-14_partner_judge_spiral_accuracy.md` (ba22d215).
- Job-owned results (not edited): `out/RESULT_auto.md`, `out/verdict.json`. Also `out/accuracy.json`, `out/armb_result.json`,
  `out/gpu_hours.json`, `out/h2_fit_compare.json`, `out/wind_count.json`, `out/s0_*.json`, `out/diag_label_*.json`.
- Code: `scripts/` (`score_accuracy.py`, `acc_lib.py`, `summarize.py`, `h2.sh`, `run_fits_v2.sh`, `arm_b.py`, `result_md.py`).
  Logs: `logs/driver.log`, `logs/fit_*.log`. Fits and data stay local and git-ignored under `fits/`, `data/`.
