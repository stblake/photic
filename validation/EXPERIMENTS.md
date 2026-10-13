# Experiments to make the semi-analytical model work on its own

The goal is depth from imagery alone: no `DEPTHS` prior at run time. Soundings may be used to *tune* settings that apply everywhere (weights, bounds, priors on water and bottom properties), never as a per-pixel input. Every experiment is scored with the benchmark (`README.md`) against `scorecards/baseline.json`.

## Baseline (commit 3aceb7a)

- **Tier 1:** the C code implements its equations exactly.
- **Tier 2:** Lee's approximation is within 2.4% (median) of OSOAA.
- **Tier 3:** inversion fails without a prior in every scenario (+6 to +77 m), including the control scenes generated with photic's own forward model and no noise. With a prior it returns the prior.
- **Tier 4:** +19 to +44 m bias at every Exmouth site; the published product (run with an empirical prior) scored RMSE 1–4 m.

The control failure means the first problem is the inversion (error function, search, heuristics), not the physics.

## Rules

- Each experiment has a hypothesis, a prediction written before it runs, a metric, and a result recorded in the log below, including failures.
- One change per branch. Accept a change to the model core only under the README's acceptance rule.
- **Tuning uses training data only.** Tier 3: tune on clear and moderate water, test on turbid. Tier 4: tune on Marina and the gulf entrance, test on Bundegi and Tantabiddi. Report the held-out scores.
- Before trusting settings tuned at Exmouth, test them at a second region.

## Phase A: diagnose (no model changes)

Tool: `drivers/invert_lab` builds photic's `model_data` for one pixel (and its spatial window), exactly as `samodel()` does, and can evaluate `samodel_error`, run Nelder–Mead from any start with any parameters fixed, or run photic's own `samodel_optimise`. `phaseA.py` drives it on tier 3 control and OSOAA pixels.

- **A0. Reproduce.** `samodel_optimise` on single control pixels, with no LUT and no carry-over between pixels, gives the same failure as the full run. If not, the failure is in the pixel loop (LUT reuse, `start_at_previous`).
- **A1. Is the truth the minimum?**
  - Conditional profile: error against H with every other parameter at its true value.
  - Profiled: at each fixed H, re-fit the other parameters (several starts) and take the lowest error.
  - If the lowest error is not at the true H, the error function is wrong (K prior, spatial-coherence terms, weights), and no optimiser can fix it.
- **A2. Start at the truth.** Run Nelder–Mead from the true parameters with photic's step sizes. Walking away means the error function prefers something else; staying means the starting values and search are the problem.
- **A3. Give it the truth piece by piece.** Fix the water (P, G, X), or the bottom (B, q), or the depth, and invert the rest. This shows which parameters the data cannot pin down.
- **A4. Same on OSOAA scenes.** Repeat A1–A3 on tier 3's OSOAA pixels, where the forward model has its ~2–3% structural error.

## Phase B: ablate the heuristics

Run tier 3 (control, then OSOAA) with one change at a time, then the promising ones together:
- LUT shortcut off; `start_at_previous` carry-over off;
- early exit (`lowest_error < 2.5·n_scenes`) off, and the order of depth starts (currently 40 m first);
- K prior off; spatial depth and bottom coherence terms off; error weights;
- NSPATIAL 1 (9 depths and 54 bottom parameters fewer); NBOTTOMS 1;
- parameter bounds (albedo ≤ 1, P, G, X ranges) in place of `fabs()`;
- the paper's χ with the 750 nm correction; the uninitialised `h_empirical` without a prior.

## Phase C: what is achievable

- **C1. Oracle optimiser:** the same error function minimised globally in Python (`photic_ref.py`, differential evolution). If the oracle succeeds, rebuild the C search to match it; if not, the error function or the parameterisation must change.
- **C2. Information limit:** Cramér–Rao bound on depth for each water type, depth, number of scenes and noise level, and the depth beyond which the bottom is invisible. Pixels past that limit should be flagged optically deep, not given a depth.

## Phase D: tune with soundings (training sites only)

- **D1. Reflectance closure at Exmouth:** at sounding cells, fix depth to the sounding and fit the rest. A good fit means the field failure is in the inversion; a poor one points at ACOLITE, tides or the bottom spectra.
- **D2. Fit the global settings:** error weights, prior strengths, bounds and ranges of albedo, P, G and X, and a correction for Lee's ~3% high bias at κH ≈ 1–4 (tier 2). Score on held-out sites.

## Phase E (optional): a self-made prior

Invert the confident pixels, calibrate a band-ratio depth model on them, and use that as the prior for the rest. Tier 3's prior scenarios measure the potential gain.

## Log

| ID | Date | Hypothesis | Prediction | Result | Decision |
|---|---|---|---|---|---|
| A0–A4 | 2026-10-13 | The inversion fails because of the search, not the physics | The truth scores lower than photic's fit on control pixels | Partly confirmed; the error function also excludes the truth in places (A-2 to A-4) | Phase B, in the order below |
| A-extra | 2026-10-13 | The albedo starting value (4·Rrs490 ≈ 0.03–0.10) traps the search | Retrieved albedo tracks the start | Confirmed: albedo ends at its start (0.03, 0.29, 0.57) even with 50,000 evaluations | Reparameterise the bottom (B1) |
| B1 | 2026-10-13 | A bounded bottom (one albedo per region plus mixing fractions) and bounded water parameters let the search fit albedo | Albedo error shrinks and depth improves | First version: noisy tier 3 scenes improved (moderate, 2 scenes: RMSE 18.1 → 5.7 m), but with albedo bounded, depth (still `|h|`) ran away to 10^13–10^20 m and corrupted tier 4. With depth bounded to [0.05, 60] m: see "Search ceiling" | Not accepted alone |
| B1 benchmark | 2026-10-13 | (B1 with depth bounded, commit 946b44d, full benchmark) | Depth improves | Tier 3: noisy median RMSE 15.9 → 12.5 m, noise-free 30.4 → 35.2 m. Tier 4: worse at every site (bias +38 to +53 m, against +19 to +31 m), because Nelder–Mead still cannot fit and depth sits near the bound | Not accepted alone; scorecards b569a8e and 946b44d |
| B2 (lab) | 2026-10-13 | The K prior is a wall at 5 m that blocks the shallow solution | Shallow pixels recover | Clear control: photic RMSE 11.9 → 7.6 m; the truth is the global minimum at 97% of pixels (was 70%) | Keep; benchmark with B5 |
| Ceiling | 2026-10-13 | The error function is good enough if the search is good | The profiled argmin is near the truth | Confirmed for clear and moderate water (table below) | Do B5 (search) next |
| B5 | 2026-10-13 | A profile search over depth (10 common depths, then refine the 2 best) with no deep-first early exit finds the minimum | Bias near 0, RMSE near the ceiling | Bias fell from +10–20 m to −1 to +2.5 m, but RMSE was still 4–9 m. With the true depth given, Nelder–Mead still left 6–12% Rrs error (20,000 evaluations) | Branch `b5-depth-search`; the optimiser must change |
| C1 | 2026-10-13 | Bounded least squares on photic's Rrs misfit (same forward model, same window, same B5 search) recovers depth | RMSE near 1 m in clear and moderate water | Confirmed: control 0.1 m (clear) and 1.9 m (moderate); OSOAA 2.3 m and 4.4 m. Given the true depth, it fits to 0.000% in about 300 evaluations | Replace Nelder–Mead with bounded Levenberg–Marquardt in C |

### Phase A results (2026-10-13)

`phaseA.py` on 33 pixels per scenario: rows from 0.8 to 25.6 m, in sand, seagrass and coral; 2 scenes; noise-free. Each run is the full per-pixel problem photic solves: both scenes jointly, over the 3×3 window of `NSPATIAL 2` (9 depths, 27 B, 27 q, 6 P/G/X = 69 parameters). Results are in `runs/phaseA/*.json`.

| Scenario | photic (A0) bias / RMSE | from truth (A2) RMSE | truth lower than photic's fit | truth is the global minimum |
|---|---|---|---|---|
| clear, control | +11.3 / 15.0 m | 0.26 m | 24/33 | ≥ 73% |
| moderate, control | +3.3 / 14.1 m | 0.23 m | 24/33 | ≥ 61% |
| turbid, control | +18.2 / 20.3 m | 0.09 m | 33/33 | ≥ 70% |
| clear, OSOAA | +9.6 / 12.9 m | 0.22 m | 24/33 | 6% |
| moderate, OSOAA | +6.5 / 14.8 m | 0.49 m | 6/33 | 0% |
| turbid, OSOAA | +18.7 / 20.6 m | 0.24 m | 7/33 | 0% |

"Global minimum" is an upper bound: the share of pixels where no profiled fit at another depth scored lower than the truth.

- **A-1. The search fails.**
  - On control pixels, the truth scores lower than photic's fit at 24 to 33 of 33 pixels.
  - Started at the truth, Nelder–Mead stays there (RMSE 0.1–0.5 m), so the error function is fine near the truth.
  - Dropping the early exit and keeping the best of all 8 depth starts (A0x) barely helps.
  - Retrieved albedo stays at its starting value (4·Rrs490, about 0.03–0.10, against a true 0.6), even with a tight tolerance and 50,000 evaluations. The search instead pairs a dark bottom with high backscatter: at 1.5 m, X = 0.20 was fitted against a true 0.02.
  - B × q is redundant: 54 bottom parameters for 9 pixels, with only their product mattering. This leaves flat and compensating directions.
- **A-2. The K prior excludes the truth in shallow clear and moderate water.** It is active at the truth at 9/33 pixels, all shallower than 5 m. At 3.7 m in moderate water, the truth fits Rrs to 1e-9 but scores 174, against 16 for photic's fit.
- **A-3. The depth-coherence term penalises real slopes.** At the truth it is active at 10/33 pixels: the deep rows, where the threshold is 5% and the slope is about 1 m per pixel.
- **A-4. The χ → Y relation is wrong for these waters.**
  - photic estimates Y from χ = Rrs440/Rrs490 of the measured spectrum. It gives Y = 0 (clamped) in moderate and turbid water, where the truth is Y = 1; in clear water it gives Y ≈ 0.7–1.2.
  - So on the OSOAA scenes the truth has a reflectance error of 2.7% (clear), 12% (moderate) and 21% (turbid). The truth is never the minimum, and the depth that best fits with everything else true is biased by about +9 m.
  - **Benchmark gap:** tier 2 passed the true Y to photic, and the control scenes are built with photic's own Y, so neither tested this step. Tier 2 needs a variant that uses photic's Y from the simulated spectrum.
- **A-5. The spatial window helps.** With `NSPATIAL 1`, photic's bias is +25 m against +3 m on moderate control. Keep it.
- **A-6. Information limit.** In turbid water the profile is flat beyond a few metres (argmin RMSE 11 m on control), so those pixels should be flagged optically deep (C2).

### Search ceiling (2026-10-13)

Per-pixel depth RMSE (m), 33 pixels, 2 scenes, noise-free. *photic* is `samodel_optimise`. *Best start* takes the best of all 8 depth starts with no early exit. *Profiled* is a thorough search: the depth is fixed on a grid from 0.25 to 40 m and everything else is re-fitted from two starts with 20,000 evaluations; it uses the same error function as photic.

| Scenario | baseline: photic / best start / profiled | B1 (bounded): photic / profiled | B1 + B2: photic / profiled |
|---|---|---|---|
| clear, control | 15.0 / 13.4 / 0.9 | 11.9 / 1.1 | 7.6 / 0.8 |
| clear, OSOAA | 12.9 / 13.4 / 0.9 | 10.7 / 1.0 | 10.0 / 0.7 |
| moderate, control | 14.1 / 14.2 / 1.9 | 16.2 / 2.0 | 7.9 / 1.8 |
| moderate, OSOAA | 14.8 / 11.7 / 3.6 | 15.3 / 1.9 | 12.1 / 3.1 |
| turbid, control | 20.3 / 9.2 / 11.7 | 22.3 / 12.0 | 18.1 / 12.0 |
| turbid, OSOAA | 20.6 / 14.7 / 9.8 | 22.3 / 11.1 | 18.9 / 11.9 |

The B1 + B2 column was run before depth was bounded.

- **The search, not the error function, is the main failure:** about a 10 m gap in clear and moderate water. That includes the OSOAA scenes despite the χ → Y error (A-4).
- **Turbid water's ~12 m is the information limit.** The deeper half of the scene is optically deep, so it needs an optically-deep flag (C2) rather than a better search.
- **B5 moves to the front.** Use an outer search over depth (a coarse grid of common depths, re-fitting the other parameters at each), then refine all parameters from the best depth. Drop the deep-first order and the early exit.

### C1: least squares in place of Nelder–Mead (2026-10-13)

`c1_oracle.py`: photic_ref's forward model (which matches the C code exactly), photic's Rrs misfit with no K prior and no coherence terms, the B1 parameterisation, and SciPy's bounded trust-region least squares, using the B5 search. Per-pixel depth RMSE, 33 pixels, 2 scenes, noise-free:

| Scenario | B5 (Nelder–Mead) | C1 (least squares) | within tol. | Rrs misfit (median) | evaluations / pixel |
|---|---|---|---|---|---|
| clear, control | 6.1 | **0.10** | 100% | 0.000% | 2,744 |
| clear, OSOAA | 8.1 | **2.33** | 64% | 0.073% | 1,822 |
| moderate, control | 5.6 | **1.91** | 91% | 0.000% | 2,565 |
| moderate, OSOAA | 4.2 | 4.38 | 39% | 0.335% | 1,760 |
| turbid, control | 6.4 | 10.8 | 67% | 0.000% | 3,086 |
| turbid, OSOAA | 9.3 | 11.6 | 21% | 0.287% | 2,528 |

- **The optimiser was the main failure.** The physics, the misfit and the B1 parameterisation are enough to retrieve depth when they are minimised properly, and in fewer evaluations than photic uses now.
- **On OSOAA scenes, the fit is close but the depth is biased.** A misfit of 0.1–0.3% with depth errors of 2–4 m means the model reaches the data with wrong parameters, which is consistent with the χ → Y error (A-4). B3 is next after the optimiser.
- **In turbid water the deep pixels are unidentifiable** (the information limit, C2). The shallow third is retrieved.

### B6: bounded Levenberg–Marquardt (2026-10-13, branch `b6-lm`)

The LM solver (`samodel_lm`) replaces Nelder–Mead, with box-bounded parameters. The bottom is now linear: one albedo per bottom type per region, B_k in [0, 1], because softmax fractions stalled approaching a pure bottom. Bounds are handled by an active set. Phase A pixels, 2 scenes, noise-free, depth RMSE (within tolerance):

| | clear | moderate | turbid |
|---|---|---|---|
| control | 0.21 m (100%) | 1.22 m (97%) | 4.3 m (91%) |
| OSOAA | 2.8 m (48%) | 4.5 m (42%) | 11.9 m (18%) |

Tier 3, moderate water, 2 scenes, noise 0.0005, by depth: 0–2 m RMSE 0.08 m; 2–5 m 0.68 m; 5–10 m 3.4 m. At 10–30 m the depth comes out far too shallow (RMSE 12–20 m) with albedo about 0.07. Below the noise, a shallow, near-black bottom fits as well as the truth (the optically-deep ambiguity).

### B6 benchmark (commit a6d4847: B1 + B2 + B5 + LM; no optically deep test, glint offset or bottom bounds)

**Tier 4, Exmouth Gulf** (baseline → B6):

| Site | bias / RMSE | 0–2 m | 2–5 m | 10–20 m |
|---|---|---|---|---|
| Marina | +31.2 / 41.4 → **−3.9 / 4.9 m** (r 0.31 → 0.67) | 0.65 | 2.21 | 8.8 (bias −8.8) |
| Tantabiddi | +19.4 / 24.9 → **−2.0 / 6.0 m** (within tol. 5% → 45%) | 0.63 | 1.57 | 12.2 |
| Bundegi | +43.7 / 133 → **−5.4 / 8.1 m** | 0.68 | 1.80 | 15.3 |
| Gulf entrance | +29.0 / 30.2 → **−13.1 / 13.3 m** | 0.66 | 1.86 | 13.7 |

- **The first working result on real imagery.** Shallow water is good at every site: 0.6–0.7 m at 0–2 m, 1.6–2.2 m at 2–5 m.
- **Deeper water comes out too shallow**: the dark-bottom ambiguity, which C2 and B8 address (not in B6).
- Tier 4 takes 5–17 minutes per site.

**Tier 3** (median RMSE, baseline → B6):
- noise-free 30.4 → 4.1 m;
- control 37.6 → 1.8 m (within tolerance 0% → 87%);
- noisy 15.9 → 11.1 m.

**Regression: the `DEPTHS` prior scenarios, 2.9 → 9.4 m.** With LM the prior is only a starting point, and the fit leaves it where depth is unidentifiable.
- **B17:** fit the prior as a residual, (H − H_prior)/σ_prior, so it weighs in where the data say little.

### Initial depth from physics (2026-10-13, `h0_estimate.py`)

The bottom signal D_i = rrs_i − rrs_dp,i = (ρ_i/π − rrs_dp,i) exp(−M κ_i H). For the red–green pair, the bottom brightness cancels:

H0 = [ln(D_560/D_665) − ln(ρ_560/ρ_665)] / (M(κ_665 − κ_560)),

where the κ difference is dominated by known water absorption. A Stumpf-style ratio with coefficients from physics, not from soundings. Where the red signal is lost, H0 falls back to a bound from 560 nm with the brightest plausible bottom.

**Water reference.** photic's starting water fails: X = 30 a_w(640) Rrs(640) reads a bright red bottom as backscatter, which inflates rrs_dp so the bottom signal vanishes (the bottom contamination of A-4 again). With the image's own deep-water reflectance (deepest rows) and κ from the water, noise 0.0005, RMSE (within tolerance):

| | 0–2 m | 2–5 m | 5–10 m | 10–15 m |
|---|---|---|---|---|
| clear, tier 3 bottoms | 0.29 (100%) | 0.50 (100%) | 1.37 | 1.91 |
| moderate, tier 3 bottoms | 0.35 (100%) | 0.37 (86%) | 1.04 (82%) | 3.01 |
| turbid, tier 3 bottoms | 0.37 (83%) | 0.51 (100%) | 13.9 | – |
| clear, realistic dark bottoms | 3.0 | 11.4 | 13.3 | 15.5 |

- **Excellent for bright bottoms, with no optimisation.**
- **It fails for dark seagrass and coral:** their reflectance is close to or below the deep-water reflectance, so D is small or negative, and the bright-bottom fallback overestimates depth.
- **Next:**
  - make H0 class-aware (try bright and dark bottom ratios, keep the consistent one);
  - estimate κ without the bottom (red-band anchor, B7);
  - use H0 to seed the search alongside the depth grid, not to replace it;
  - then measure the speed and accuracy gain.

### Class-aware H0 and seeding the search (2026-10-13, `samodel_h0`, `SEED_MODE`)

**Estimator.** The water comes from each scene's deep-water Rrs (SCENE RINF), QAA-style:
1. u from rrs_dp;
2. at the red band a ≈ a_w, which fixes b_b;
3. Y from QAA's ratio;
4. then κ = a + b_b at every band.

For each bottom class (photic's shape, brightness within the measured class range) and each depth on a log grid, the brightness solves a linear least-squares problem: rrs_i = rrs_dp,i + (b s_i/π − rrs_dp,i) exp(−M κ_i H). The best (depth, class) wins, so dark bottoms (signal below deep water) are handled. Scenes are combined by the median after removing tides.

Prototype (`h0_estimate.py main_class`, image deep rows as RINF), class-aware H0 alone:
- realistic bottoms, 0–2 m: 0.11–0.57 m, with the class right in every case (the ratio method was off by 3 m or more);
- 5–15 m: 2–5 m RMSE, with a negative bias at depth from the QAA κ.

**Seeding.** `SEED_MODE` 1 adds H0 to the 10-depth grid; 2 uses {H0/2, H0, 2H0, 30 m}. Pixels use candidate 2's physics (Y2, bottom bounds), and RINF is the OSOAA no-bottom Rrs at 32 m. RMSE; coverage is unchanged across modes:

| | no seed | grid + H0 | H0-centred |
|---|---|---|---|
| clear, default bottoms, noise | 3.59 | 3.00 | **1.97** |
| moderate, default bottoms, noise | 4.71 | 4.67 | **2.82** |
| clear, realistic bottoms, noise | 2.08 | **0.58** | 0.59 |
| moderate, realistic bottoms, noise | 1.14 | 1.14 | **0.87** |
| moderate, default bottoms, no noise | **1.52** | 3.47 | 1.40 |
| moderate, realistic bottoms, no noise | **1.87** | 1.86 | 3.18 |
| clear, realistic bottoms, no noise | 1.68 | **1.05** | 2.03 |
| s/pixel | 0.34–0.64 | 0.39–0.61 | **0.24–0.41** |

- **With noise, H0-centred seeding is consistently better, and about 30% faster.**
- **Without noise it is mixed.** The misfit has several minima of similar size, because the forward model is slightly wrong; a different seed can reach a minimum with lower misfit but the wrong depth (grid + H0, moderate).
- **Decision:** seeding stays off unless RINF is given.
- **Next:**
  - pass RINF in tier 3 (OSOAA no-bottom Rrs);
  - estimate RINF in tier 4 from deep-water pixels (offshore of Ningaloo);
  - judge mode 2 on the full benchmark.

### Decoupling the water from depth and bottom, steps 1–2 (2026-10-13, `decouple.py`)

photic's search (profile, then refine; LM; `Y_MODEL 3`, `BOTTOM_BOUNDS 1`) on OSOAA pixels in three ways:
- **joint:** everything fitted;
- **true:** the water (P, G, X, Y per scene) held at the truth, the ceiling;
- **deep:** each scene's water fitted first to its deep-water Rrs (photic's model, no bottom), then held fixed for the shallow pixels.

The deep-water fit recovers P, G and X within about 5–40% and Y within about 30%; clear water splits P and G worst. With the bounds, the default scenes' seagrass and coral (albedo 0.6) are out of range, so only the realistic-bottom scenes are a fair test. RMSE:

| Realistic bottoms, noise 0.0005 | 2–5 m | 5–10 m | 10–15 m |
|---|---|---|---|
| clear: joint / **deep** / true | 1.09 / **0.35** / 0.28 | 2.19 / **1.22** / 0.92 | 3.08 / **1.63** / 1.77 |
| moderate: joint / **deep** / true | 0.56 / 1.00 / 0.86 | 3.32 / **1.77** / 1.65 | 7.92 / **5.23** / 5.82 |
| turbid: joint / **deep** / true | 2.05 / 5.53 / 0.48 | 6.34 / **4.03** / 3.95 | – |

| Realistic bottoms, no noise | 2–5 m | 5–10 m | 10–15 m |
|---|---|---|---|
| clear: joint / deep / true | 0.11 / 0.15 / 0.11 | 0.96 / 0.89 / **0.25** | 1.88 / 1.96 / **0.64** |
| moderate: joint / deep / true | 0.84 / **0.11** / 0.36 | 3.00 / **1.23** / 0.99 | 3.06 / 4.30 / 3.62 |

- **With noise, water from deep water gains a lot at 2–15 m**, often reaching the ceiling: a noisy pixel cannot pin down its own water, but the deep-water reference can.
- **Exceptions:** turbid 2–5 m (a few outliers) and moderate 2–5 m.
- **Without noise it is mixed.** The joint fit can bend the water to absorb the forward model's own error, while fixed water pushes that error into depth. In clear water the deep-water fit's P/G split (G 0.013 against 0.02) costs it the ceiling's gain.
- **Caveat:** tier 3 is the best case, with uniform water and a deep reference of the same water. Real lagoons can differ, hence step 3, a water shared over a tile within a scene and fitted with its pixels.

### Decoupling, step 3: water shared over a tile (2026-10-13, `tile_fit.py`, negative on tier 3)

**Method.** Per 8×8-pixel tile and per scene, one water (ln P, ln G, ln X, Y). The fit alternates two steps:
- **pixel step:** water fixed, each pixel fits depth and bottom in its 3×3 window, as photic;
- **water step:** depths and bottoms fixed, the tile's water is fitted to all 64 pixels.

Round 1 from the deep-water water reproduces decouple.py's deep result exactly. Realistic bottoms, noise 0.0005, RMSE (mean misfit):

| | 2–5 m | 5–10 m | 10–15 m | misfit |
|---|---|---|---|---|
| clear: deep-water water, fixed | 0.35 | **1.22** | **1.63** | 7.40% |
| clear: … then refitted on the tile (3 rounds) | 0.34 | 1.71 | 1.98 | 7.26% |
| clear: tile fit from photic's starting water (4 rounds) | 1.60 | 3.00 | 6.06 | 11.95% |
| moderate: deep-water water, fixed | **1.00** | **1.77** | **5.23** | 7.50% |
| moderate: … then refitted on the tile | 1.63 | 2.68 | 6.77 | 7.03% |

- **Refitting the water to the shallow pixels lowers the misfit but worsens depth.** The water bends to absorb the forward model's own 2–3% error (tier 2), and the depths pay.
- **The deep-water fit never sees the bottom term, so its water stays closer to the truth.**
- **From photic's bottom-contaminated starting water, the alternation stalls in a poor state.** Restricting the water step to pixels fitted shallower than 8 m made no difference.
- **Decision:**
  - prefer water from a deep-water reference when one exists (step 2);
  - otherwise keep the joint fit;
  - do not adopt the tile-shared water.
- **Before revisiting it:**
  - a forward model with smaller structural error (the root cause);
  - a water prior towards the deep-water estimate rather than a free refit;
  - a joint solver (Schur complement, as in bundle adjustment) instead of alternation;
  - and judge it on tier 4, where the water really varies.

### Reducing the forward model's structural error (2026-10-13)

**Where the error is** (tier 2, OSOAA, black sky). Lee's deep-water rrs has no sun-angle dependence, and OSOAA's does:
- deep-water bias +4.9% at θ_w 15°, +2.4% at 29°, −2.0% at 40°;
- shallow error follows the same pattern, more weakly (+2.0 / +0.9 / −1.2%);
- bottom albedo hardly matters (+0.6% bias at every albedo).

Scenes with different sun angles therefore get inconsistent model errors, which a multi-scene fit can only absorb through depth.

**What a perfect forward model would be worth.** The same realistic-bottom pixels, scenes from photic's own model against OSOAA (RMSE on the pixels given a depth):

| | photic's model | OSOAA |
|---|---|---|
| clear, no noise | 0.08 | 1.41 |
| moderate, no noise | 2.82 (only 12 of 48 pixels) | 1.31 |
| clear, noise | 0.54 | 0.48 |
| moderate, noise | 1.34 | 1.10 |
| turbid, noise | 0.22 | 0.54 |

- Without noise the model error costs about 1.3 m in clear water.
- **With realistic noise, the gap disappears**: noise and coverage limit depth, not the model.

**Candidate models**, against tier 2's OSOAA cases; median error, bias by sun angle (θ_w 15 / 29 / 40°):

| | clear + moderate (training) | turbid (held out) |
|---|---|---|
| Lee (photic) | 1.91%, +0.8 / +0.3 / −1.0 | **3.12%**, +4.9 / +2.3 / −1.6 |
| Lee × (1 + a(sec θ_w − 1.1)), a = 0.083 | 1.75%, +0.5 / +0.5 / +0.1 | 2.84%, +4.4 / +2.6 / −0.2 |
| Albert & Mobley (2003) as published (Eqs. 9, 13; Tables 3–4) | 3.41%, deep +8.1% | 10.6%, deep +13.2% |
| AM03, scale refitted (p1, A2) | 3.34% | 4.32%, deep +5.6% |
| AM03, 7 coefficients refitted | 1.44% | 4.74%, deep +7.6% |

- **Albert & Mobley do not beat Lee against OSOAA.**
  - Their model was fitted to Hydrolight runs of Lake Constance, with skylight, Raman scattering and fluorescence, which add light: hence the high bias.
  - Their deep-water polynomial drifts in turbid water.
  - Fully refitted, it overfits.
- **The sun-angle effect is sky-dependent.** AM03's published angular term (p5 = 0.11, Hydrolight with skylight) is much weaker than black-sky OSOAA's (the refit wants p5 ≈ 0.5). A black sky sends all light from the sun's direction and exaggerates angular effects. A correction fitted to tier 2/3 would overcorrect real imagery.
- **Decision:**
  - no change to the forward model now;
  - make the synthetic truth realistic first, with OSOAA run with an atmosphere so skylight reaches the water, before fitting any angular correction;
  - Lee's structure remains the most transferable. Refits (B10) and an angular factor are available once the truth includes skylight.
- The paper is stored in `raw/papers/`.

### Sky-lit synthetic truth (2026-10-13; correction to the section above)

**Builds.** `OSOAA_sky` is the scalar build with atmospheric scattering on (`EXPERT_MODE_FORCED_FATM_NULL .FALSE.`); `OSOAA_sky_blackocean` also has the sea source nulled, giving the surface-reflected sky and glint to subtract. The aerosol is maritime (Shettle & Fenn, relative humidity 80%), with optical thickness 0.1 at 550 nm, passed as `run_osoaa(aot=)`.

**Normalisation.** OSOAA's REFL = π I/Ed_total(level), with Ed_total direct plus diffuse at the same level (OSOAA_CALC_REFL). So REFL(0+)/π − surface is water-leaving Rrs relative to the total Ed(0+), with or without a sky. The earlier tier 2/3 numbers were correctly normalised.

**Tier 2, all 1,098 physical cases, sky against black sky:**

| | black sky | sky, AOT 0.1 |
|---|---|---|
| OSOAA rrs(0−), change from black sky | – | median −0.03% (5–95%: −0.9 to +0.8%) |
| OSOAA Rrs(0+), change from black sky | – | median −0.8% |
| Lee error, median / p95 | 2.43% / 6.68% | 2.07% / 6.55% |
| Lee bias by sun angle (θ_w 15 / 29 / 40°) | +2.2 / +1.0 / −1.2% | +2.3 / +1.0 / −1.4% |
| Lee deep-water bias by sun angle | +4.9 / +2.3 / −2.0% | +4.4 / +2.0 / −1.8% |

- **Skylight (14–33% of Ed) hardly changes the truth, and the sun-angle dependence is intact.**
- **The earlier explanation was wrong:** the black sky does not exaggerate OSOAA's angular effect.
- **The mismatch with Albert & Mobley's angular term must come from elsewhere.** The likely candidate is the particle phase function: Fournier–Forand with b̃_b = 0.018 here, against Petzold-like in Hydrolight. Their Lake Constance optics are another.
- **Race fixed:** aerosol Mie files are cached per wavelength. Parallel OSOAA runs raced to write them, so tier 2 now warms up one case per sun angle and band, and tier 3 one per band.
- **Tier 3:** `PHOTIC_TIER3_SKY=1` builds and uses sky-lit lookup tables (`lut_sky`).
- **Decision:** use the sky-lit truth from now on, since it is more realistic at no extra cost in photic.
- **Next for the forward model:** test the phase-function hypothesis (OSOAA with a Petzold-like phase function) before fitting any angular correction.

### Physics constraints (planned 2026-10-13)

Each enters the LM fit as a bound or a weighted residual (a soft prior, weight 1/σ), and is switchable. Each is tested alone: Phase A pixels with noise, then tier 3 by depth bin, then tier 4 held-out sites. Any σ is tuned on the training split only.

- **B3. Self-consistent Y.** photic derives Y from χ = Rrs440/Rrs490 of the *measured* spectrum, which in shallow water includes the bottom (A-4: Y = 0 against a true 1).
  - Compute Y instead from the model's own deep-water reflectance with QAA's relation, Y = 2.0[1 − 1.2 exp(−0.9 rrs_dp(443)/rrs_dp(560))], as a fixed point inside the forward model.
  - Variant: fit Y per scene, bounded to [0, 2.5].
  - Prediction: the OSOAA depth bias falls (moderate OSOAA RMSE 4.5 m → about 2 m); the control scenes are unchanged only if they are rebuilt with the same Y.
- **C2. Bottom detectability.**
  - From the forward model, rrs(H) − rrs_dp ≈ (ρ/π − rrs_dp) exp(−M K H). The bottom is detectable in band λ only while that contrast exceeds zσ, so H_max = max_λ ln(|ρ/π − rrs_dp|/(zσ))/(M K).
  - Report pixels whose solution exceeds H_max, or whose contrast is below the noise at every band, as optically deep ("deeper than H_max") rather than giving a depth.
  - Also report the depth uncertainty from the LM Jacobian, σ_H² = σ_r²[(JᵀJ)⁻¹]_HH, which could replace the Monte Carlo depth-σ pass.
  - Prediction: the shallow, black-bottom solutions are flagged, and RMSE on the reported pixels falls, at a stated loss of coverage.
- **B7. Red and near-infrared backscatter anchor.**
  - Where a ≈ a_w (665 nm and beyond) the bottom is invisible below a few metres, and rrs fixes b_b, hence X.
  - Add Sentinel-2 B5 (705 nm, and possibly B6 at 740 nm) as fitted bands, after renormalising the bottom spectra (which exceed 1 at 705 nm), with a larger σ to allow for atmospheric correction error.
  - Prediction: the trade-off between a dark bottom and high X disappears, and deep pixels stop collapsing shallow.
- **B8. Physical bottom library.** ρ(λ) = s Σ f_k R_k(λ), using absolute reflectances, fractions on the simplex and a brightness scale s in [0.7, 1.3]; ρ ≤ 1 at every band. A dark bottom must then have a dark spectral shape.
- **B4. Seabed continuity.** Penalise the second difference of depth across the window (a plane is free), in place of departure from the window mean (A-3).
- **B9. Water-mass priors.** Weak priors keeping P, G and X per scene within natural ranges, from the training data.

### Lee's HOPE spreadsheet (2026-10-13)

`Lee_2_QAA_HOPE_hyper.xls` (Z.-P. Lee, 2009/2014) has QAA, deep HOPE and shallow HOPE sheets, with every intermediate value. Recomputing the columns shows the same equations as photic, with these differences:

| Term | photic (Lee 1999) | Lee's spreadsheet (Lee 2002) |
|---|---|---|
| deep-water rrs | (0.084 + 0.170u)u | (0.089 + 0.125u)u |
| Rrs from rrs | 0.5 rrs/(1 − 1.5 rrs) | 0.52 rrs/(1 − 1.7 rrs) |
| D_B (bottom) | 1.04 √(1 + 5.4u) | 1.05 √(1 + 5.5u) |
| Y | χ relation on the measured Rrs440/Rrs490 | QAA's relation on the measured rrs (443/555), from the QAA sheet |
| sand | `bottom_type_sand` | `bottom_sand_Lee` (in photic's tables, unused) |

- **Identical:** a_φ, a_dg (S = 0.015), D_C, the sun path (secant 1.2) and the error (√mean d²/mean Rrs).
- **Tier 1 now checks the spreadsheet:** `photic_ref` with the spreadsheet's coefficients reproduces its modelled Rrs to 10⁻¹⁶ (deep) and 2 × 10⁻⁵ (shallow), so the equations' structure is identical. With photic's coefficients, Rrs differs by 7–9%.
- **Against exact radiative transfer (tier 2's OSOAA cases), photic's 1999 set is closer:**
  - subsurface rrs: median 1.91% vs 2.60% (clear and moderate), 3.12% vs 3.26% (turbid);
  - Rrs from rrs: 1.7% vs 6.3% (bias +1.5% vs +6.7%).
  - The 2002 set was fitted to Hydrolight with other surface and sky assumptions than OSOAA's black-sky, low-wind runs. Tier 4 decides on real data.
- **QAA's Y on the measured spectrum, Lee's own practice (`Y_MODEL 4`).** Tested against Y0 (photic now) and Y2 (self-consistent χ) on the OSOAA pixels, RMSE:

| | Y0 | Y2 | Y4 |
|---|---|---|---|
| clear, no noise, 5–10 m | 1.70 | **0.59** | 0.98 |
| clear, no noise, 10–15 m | 3.52 | **1.74** | 2.19 |
| moderate, no noise, 5–10 m | 3.12 | **2.49** | 3.05 |
| turbid, no noise, 2–5 m | 2.25 | **1.19** | 2.20 |
| moderate + noise, 5–10 m | 2.74 | **1.91** | 2.29 |
| turbid + noise, 2–5 m | 1.76 | 1.30 | **1.28** |
| clear + noise, 5–10 m | 1.27 | 1.51 | **1.17** |

- **Y2 is best overall.** Keeping the bottom out of the Y estimate halves the clear-water deep bias (10–15 m: +3.4 → +1.65 m), which explains part of the residual OSOAA bias.
- **Lee's practice (Y4) sits between Y0 and Y2.**
- **Decision:** benchmark Y2 as a second integrated candidate. Tier 3's control scenes are generated with Y0, so they will regress by construction; judge Y2 on the OSOAA scenarios and tier 4.

### Error budget (2026-10-13, `error_budget.py`)

photic's search (profile, then refine, by LM; Y fitted per scene) on tier 3's OSOAA pixels, with one component held at its true value. Depth bias (RMSE), m.

Clear, no noise:

| held true | 5–10 m | 10–15 | 15–20 | 20–30 |
|---|---|---|---|---|
| nothing | +0.59 (1.01) | +2.51 (2.53) | +3.25 (3.29) | +3.70 (3.73) |
| Y | +0.42 (0.64) | +1.82 (1.85) | +2.85 (2.87) | +3.30 (3.31) |
| water (P, G, X, Y) | +0.23 (0.25) | +0.52 (0.53) | +0.79 (0.83) | +1.31 (1.39) |
| **bottom** | **+0.12 (0.13)** | **+0.07 (0.08)** | **+0.09 (0.10)** | +0.28 (0.30) |
| Y + bottom | +0.08 (0.09) | +0.05 (0.06) | +0.08 (0.09) | +0.31 (0.32) |

Moderate, noise 0.0005:

| held true | 5–10 m | 10–15 | 15–20 |
|---|---|---|---|
| nothing | −1.43 (3.50) | −8.16 (8.40) | −12.85 (12.96) |
| water | +0.94 (1.38) | −6.02 (6.31) | −9.08 (9.41) |
| **bottom** | **−0.16 (0.52)** | **0.00 (0.90)** | +2.16 (5.10) |

- **The residual bias is the depth–albedo trade-off.** Given the bottom, depth is right to 20–30 m in clear water and to 15 m in noisy moderate water.
- **The fit absorbs the forward model's small error differently by water type:** in clear water as a brighter bottom placed deeper; in noisy moderate water as a darker bottom placed shallower.
- **Y and the water matter much less.**
- **Turbid water beyond about 5 m stays unresolvable whatever is given**, the information limit flagged by C2.
- **Anomaly:** the 0–5 m bin with the bottom held true has one outlier pixel (+2.7 m mean), to check.
- **Next:** constrain the bottom.
  - **B16, bottom coherence:** a soft prior pulling a pixel's B_k towards the window mean (the seabed locally uniform at about 60 m).
  - **B8, realistic albedo ranges per bottom type:** needs absolute library reflectances.

### B16: bottom coherence across the window (2026-10-13, `COHERENCE`, mixed)

Prior residuals (B_r,k − mean B_·,k)/σ_B pull each region's bottom albedos towards the window mean. OSOAA pixels, RMSE:

| | off | σ 0.2 | σ 0.1 | σ 0.05 |
|---|---|---|---|---|
| clear, no noise, 5–10 m | 1.70 | 0.84 | 0.63 | **0.45** |
| clear, no noise, 10–15 m | 3.52 | 2.29 | 2.14 | **1.98** |
| clear, no noise, 20–30 m | **3.77** | 5.00 | 4.86 | 4.75 |
| moderate + noise, 2–5 m | 1.03 | 0.23 | **0.13** | 0.15 |
| moderate + noise, 5–10 m | **2.74** | 2.91 | 3.55 | 3.21 |
| clear + noise, 10–15 m | **2.00** | 2.23 | 2.42 | 2.67 |

- **It helps shallow and mid depths, and hurts the deepest and noisy mid-depth bins.**
- **Why it is limited:** coherence makes a window's regions agree with each other, but the whole window can still trade brightness against depth together. The error budget showed that the *absolute* bottom brightness is what removes the bias.
- **Decision:** off by default (σ 0.2 is the mildest if it is wanted). **Next:** B8 with a bottom library built from first principles, which bounds absolute brightness per type.

### B8 step 1: bottom reflectance from first principles (2026-10-13, `bottom_rt.py`)

Physical models of submerged bottom reflectance, each a family over parameter ranges:
- **sand:** Kubelka–Munk semi-infinite layer, with k/s from mineral, iron-oxide, organic and microphytobenthos absorption;
- **seagrass:** PROSPECT-D leaves with the leaf–water relative refractive index, in a Kubelka–Munk canopy over sand;
- **coral:** zooxanthellae tissue over an aragonite skeleton.

Absorption data:
- PROSPECT-D coefficients (from the `prosail` package data);
- extracted-pigment spectra from Clementson & Wojtasiewicz (2019, CC BY 4.0), with in-vivo shifts.

Both are stored in `raw/bottom/`.

| 443 / 490 / 560 / 665 / 705 nm | model, shape (÷560) | photic's table |
|---|---|---|
| sand | 0.80 / 0.90 / 1.01 / 1.08 / 1.10 (absolute 0.47–0.64) | 0.75 / 0.88 / 1.01 / 1.09 / 1.19 |
| seagrass | 0.04 / 0.08 / 0.92 / 0.07 / 1.40 | 0.55 / 0.65 / 1.03 / 0.81 / 3.33 |
| coral | 0.09 / 0.12 / 1.11 / 0.62 / 1.41 | 0.37 / 0.52 / 1.16 / 0.80 / 2.91 |

- **Sand reproduces photic's measured shape from first principles**, with an absolute reflectance (0.4–0.7) typical of bright carbonate sand.
- **Seagrass and coral are far too dark in the blue and red**, because the assumed pigment loads and full canopy cover absorb almost everything at chlorophyll's peaks. Measured spectra include canopy gaps, the sand beneath and epiphytes.
- **The models' structure is sound, but their parameter ranges are not yet anchored.** Step 2 needs measured benthic spectra to calibrate the ranges (pigment load, cover fraction, packaging), first to see whether the families then contain photic's shapes.
- **No measured library is at hand:** sambuca_core ships none.

### B8 step 2: calibration against measured spectra (2026-10-13)

**Measured data**, both open on PANGAEA, stored in `raw/bottom/measured/`:
- **Heron Reef:** in situ underwater reflectance of 105 benthic types (Roelfsema & Phinn 2012, PANGAEA.804589).
- **Moreton Bay:** underwater seagrass and sediment (Roelfsema, Phinn & Joyce 2016, PANGAEA.864310, CC BY 3.0). Samples measured in air are excluded.
- Data quality: Heron reads 0 at 400 nm, and every spectrum is noisy beyond about 680 nm (water absorbs the red), so the fits use 410–680 nm.

**Fit:** each measured spectrum, by its family with bounded least squares. The families were widened from step 1:
- seagrass: a canopy cover fraction and epiphyte (brown) absorption;
- coral: a darker skeleton range and broad host absorption;
- sediment: mineral background absorption from bright sand to mud.

| Class (n) | Fit error, median (90th pct) | Calibrated parameters (median) | Measured R at 560 nm, 10th–90th pct |
|---|---|---|---|
| sediment (5) | **6.6%** (10.2%) | q_mineral 0.09, iron oxide 0.55, microalgae | 0.24–0.69 |
| seagrass (6) | **7.7%** (11.3%) | Cab 14 µg cm⁻², N 1.24, LAI 1.3, epiphytes (Cbrown 1.4), cover ≈ 1 | 0.05–0.10 |
| live coral (40) | 17% (29%) | chl about 23 mg m⁻² (with packaging), skeleton under tissue about 0.2, brown τ 0.34 | 0.04–0.20 |

- **The physical families reproduce real spectra.** Coral's larger error is mostly Heron's 0.01 quantisation at R ≈ 0.02–0.06.
- **Absolute brightness separates the classes:** sediment is bright, while seagrass and coral are dark. That is the constraint the error budget asks for.
- **Benchmark gap:** tier 3's truth gives every bottom albedo 0.6, so its "seagrass" has R(560) ≈ 0.6 against 0.05–0.10 measured. Realistic bounds would make tier 3's truth unreachable.
- **Step 3:**
  1. Add a tier 3 scene set with realistic albedos per class.
  2. Parameterise photic's bottom as cover fractions and a brightness per type, bounded by the measured classes, with the shapes from the calibrated families or photic's tables.
  3. Test with `error_budget.py` and the pixel harness.

### B8 step 3: bottom brightness bounded by class (2026-10-13, `BOTTOM_BOUNDS`)

**Parameterisation.** Each region has cover weights w_k in [0, 1] and a brightness b_k per type, with ρ = Σ (w_k/Σw) b_k ρ_k(λ). b_k is bounded by the in situ classes (step 2, 10th–90th percentile at 560 nm, widened): sediment 0.15–0.75, seagrass 0.03–0.15, coral 0.03–0.25.
- A pure type still lies on the bounds (the other weights 0), but "sand at albedo 0.07" is now impossible.
- The default (0) is the linear B_k in [0, 1].
- **Bug fixed:** the optically deep test's k now counts the bottom's degrees of freedom, n_regions(1 + n_bottoms), whatever the parameterisation; bounds had doubled it.

**Realistic scenes.** `tier3_scenes.ALBEDO_BY_BOTTOM` sets each bottom's albedo from the medians of the measured classes (sand 0.5, seagrass 0.07, coral 0.10). The benchmark's default scenes still use 0.6 for every bottom, which is unrealistic for seagrass and coral.

**Results** (Phase A OSOAA pixels, 2 scenes, noise 0 and 0.0005 pooled; depth and Δχ² recorded for every pixel, then thresholded). Coverage / RMSE / within tolerance:

| | Training (clear + moderate) | Turbid (held out) |
|---|---|---|
| 1×BIC, no bounds | 41% / 2.20 m / 60% | 26% / 9.9 m / 28% |
| 1×BIC, **bounds** | 41% / **1.33 m** / 73% | 26% / 9.9 m / 36% |
| 2×BIC, no bounds | 29% / 2.23 m / 62% | 8% / 1.14 m / 62% |
| 2×BIC, **bounds** | 29% / **1.21 m** / 77% | 8% / **0.55 m** / 88% |

- **The bounds cut RMSE by 40–50% at equal coverage**, and the gain holds on held-out turbid water.
- **DEEP_PENALTY stays 2.** 1×BIC adds training coverage but fails on turbid water.
- **Coverage is limited by physics:** dark seagrass (R ≈ 0.07) has little contrast against the water column, so only 14% of seagrass pixels get a depth even at 1×BIC (RMSE 0.48 m). Sand is covered to 20 m.
- **Decision:** adopt `BOTTOM_BOUNDS` in the next candidate.
- **Caveats:**
  - The class bounds come from two Australian sites; check them against Exmouth's real bottoms (tier 4).
  - The tier 3 benchmark needs a realistic-albedo scene set, otherwise its 0.6 seagrass and coral make the bounds look like a regression.

### Geometry constraints across scenes (planned 2026-10-13)

Each scene sees the same seabed through a different light path, M_s = sec θ_w,s + D sec θ_v,s, and a different water column H + τ_s. The forward model already contains both, but P, G and X are free per scene, so a change in water can absorb most of the angular signal.

- **B12. Water tied between scenes.**
  - Add prior residuals (ln P_s − ln P̄)/σ, and likewise for G and X, with σ ≈ 20–30% (tier 3's scenes differ by about ±20%).
  - The scenes' path differences then constrain κH independently of albedo.
  - Prediction: lower depth bias at 5–20 m, and more gain with 4 scenes than with 2.
- **B13. Angular-ratio depth.**
  - With the same water in two scenes s and t, the bottom contrast Δ = rrs − rrs_dp obeys ln(Δ_s/Δ_t) = −κ[(M_s − M_t)H + M_s τ_s − M_t τ_t]. The albedo cancels, so H follows in closed form.
  - First compute its Cramér–Rao precision against ΔM, κH and noise (it amplifies noise at depth). Then test it as a starting depth and as a consistency flag against the full fit.
- **B14. Direct and diffuse downwelling light.**
  - The path factor should be sec_eff = f sec θ_w + (1 − f)/μ̄_d, with f the direct fraction (from sun elevation and aerosol optical thickness) and μ̄_d ≈ 0.86 for skylight. The current model puts all the light at θ_w, which overstates the path differences between scenes at low sun.
  - The tier 3 OSOAA scenes are black-sky, so test on tier 4, with f from ACOLITE's aerosol optical thickness.
- **B15. Glint as a per-scene offset.** Restore the per-scene spectrally flat offset (`DELTA`), with its size limited by a Cox–Munk prior from each scene's sun–view geometry. Test on tier 4 (tier 3 has no glint).

### B3 result (2026-10-13, `Y_MODEL` switch in samodel.c)

Phase A pixels, 2 scenes, noise-free, depth RMSE. Control scenes are generated with Y from the measured χ (`Y_MODEL 0`), so the other variants cannot reproduce them exactly; the OSOAA scenes are the real test.

| | Y0 measured χ | Y1 QAA, self-consistent | Y2 χ, self-consistent | Y3 fitted per scene |
|---|---|---|---|---|
| clear OSOAA | 2.80 | 1.93 | 1.94 | 2.34 |
| moderate OSOAA | 4.46 | 4.36 | 4.15 | 5.03 |
| turbid OSOAA | 11.85 | 11.81 | 11.84 | 5.87 (albedo → 0: the bottom is dropped) |
| s/pixel | 0.06–0.10 | 0.3–0.7 (before caching) | 0.10–0.21 | 0.09–0.16 |

- **Prediction failed for moderate water** (4.5 → about 2 m expected; got 4.15). Self-consistent Y helps clear water (2.8 → 1.9 m) and mid depths a little.
- **QAA's relation is not a fix on its own:** it gives Y = 1.7, 0.5 and 0.27 for the clear, moderate and turbid tier 3 waters, whose true Y is 1.
- **Moderate OSOAA by depth (Y0 → Y2):**
  - 0–5 m: RMSE 0.07–0.16 m;
  - 5–10 m: 3.1 → 2.4 m (too deep, albedo 0.76–0.81);
  - 10–15 m: 1.8 m;
  - 15–30 m: 4.6–8.0 m (too shallow, albedo 0.15).
- **The 5–10 m bias is Lee's approximation.** Tier 2 found it about 3% high at κH ≈ 1–4, and the fit absorbs that by trading depth against albedo.
- **The deep error is the optically-deep ambiguity**, with the model's own 2–3% structural error acting as noise.
- **New candidate, B10:** refit Lee's shallow-water coefficients (D_uC, D_uB) to exact radiative transfer, training on tier 2's clear and moderate cases and testing on turbid and tier 3.
- **Decision:** keep the switch (default 0), and take Y2 into the final combination if it helps on the benchmark.

### C2 result (2026-10-13, `OPTICALLY_DEEP` in samodel.c)

After the search, the water is re-fitted with every region at H_MAX and no bottom (the optically deep model). The bottom counts as detected only if Δχ² = n(E_deep² − E_best²)/(100σ_rel)² > DEEP_PENALTY·k·ln n:
- k is the number of depth and bottom parameters and n the number of residuals;
- σ_rel² = (σ_noise/mean Rrs)² + ε², with ε = 0.03, the forward model's own error from tier 2.

Pixels failing the test get no depth, and coverage is reported. A plain Δχ² > 9 test flagged nothing with noise, because the shallow model's 36 extra parameters fit the noise (Δχ² ≈ 70–120 where the bottom is invisible). DEEP_PENALTY was chosen on clear and moderate pixels and checked on held-out turbid ones (Phase A pixels, 2 scenes):

| Threshold | Moderate + noise 0.0005 (training): RMSE / within tol. / coverage | Turbid + noise (held out) |
|---|---|---|
| none | 10.3–10.9 m / 29–31% / 100% | 12.4–12.6 m / 12–15% / 100% |
| 1 × k ln n (BIC) | 2.7–3.2 / 54–60% / 52–54% | 7.2–10.0 / 24–44% / 33–52% |
| **2 × k ln n** (chosen) | **1.7–1.9 / 62–71% / 44%** | **1.1–1.2 / 55–70% / 21–23%** |

Without noise, 2 × k ln n also cuts moderate water off at about 11 m (the true limit is about 15 m), and clear water at about 24 m. Tier 3 now passes the simulated noise to photic as RSIGMA. For tier 4 no noise estimate exists yet, so only ε applies; estimating noise from the image is a follow-up.

### B7 result (2026-10-13, `b7_nir.py`)

The same OSOAA scenes were run with the four usual bands and with Sentinel-2 B5 (705 nm) added; the comparison is paired (same pixels, same noise). Albedos are realistic per bottom (sand 0.6, seagrass 0.2, coral 0.25). photic now bounds each bottom type's albedo so that its reflectance is at most 1 at every fitted band (B_k ≤ 1/max ρ_k), because seagrass and coral exceed 1 at 705 nm. The scenes use 2 scenes, LM and the optically deep test (C2). RMSE (within tolerance), coverage:

| | 4 bands | 5 bands (+705 nm) |
|---|---|---|
| clear + noise, 2–5 m | 1.26 m (43%) | **0.36 m (100%)** |
| clear + noise, 10–15 m | 2.57 m (0%) | **1.21 m (100%)** |
| clear + noise, coverage | 21/48 | 24/48 |
| clear, no noise, 10–15 m | 3.52 m | 2.65 m |
| moderate + noise, 5–10 m | 2.28 m | 1.48 m |
| moderate, no noise, 5–10 m | 2.87 m (8/11 covered) | 2.06 m (10/11) |
| turbid + noise, 2–5 m | 1.50 m (33%) | **0.57 m (100%)** |

**Decision:** adopt. The 705 nm band constrains backscattering and lifts both accuracy and coverage in every water type. Tier 4 should give photic B5 as a fifth band; ACOLITE already outputs it. The tier 3 benchmark scenes stay at four bands with albedo 0.6, so a 705 nm tier 3 set would need these realistic albedos.

### B11: bottom–surface interreflection (2026-10-13, `INTERREFLECTION` switch, negative)

OSOAA's bottom response saturates, R = R0 + Aρ/(1 − cρ), with c up to 0.42 in clear water at 0.5 m, while Lee's model is linear in ρ.
- **Derivation:** c is the reflectance seen from the bottom looking up: c = r_i exp(−2 D K H) + π rrs_dp (1 − exp(−2 D K H)), the interface's internal reflectance over the round trip plus the water column's reflectance.
- **Fit to OSOAA:** clear and moderate water give r_i = 0.464 (theory about 0.48) and D = 2.07; held-out turbid water is fitted to 0.02 rms.
- **Effect on depth:** with the term in the forward model, shallow pixels get slightly worse (0–2 m bias +0.06 → −0.12 m) and the deep bias is unchanged. Lee's empirical coefficients already absorb this effect at typical albedos, so adding it double-counts.
- **Decision:** off by default.

**Unexplained so far:** the remaining OSOAA depth bias (clear 10–30 m about +3.5 m with too bright an albedo; moderate 5–10 m about +2.7 m). It survives B3, B10 and B11.

### B4, B8, B9 (2026-10-13)

- **B4, seabed continuity (`CONTINUITY` switch, negative).**
  - Second differences of depth across the 3×3 window were added as prior residuals, with σ = σ0 + σH·H̄ in three strengths (σ0 0.03–0.3 m, σH 0.015–0.15).
  - On noisy OSOAA pixels there is no meaningful effect: clear 10–15 m RMSE goes from 2.00 m to 2.29–2.53 m, and other bins change by ≤ 0.3 m.
  - The reported depth is already a weighted mean over the window, which averages the noise much as the prior would. Off by default.
- **B8, absolute bottom library: deferred.**
  - Its essential physics, ρ ≤ 1 at every band, is in through B7. The failure it targeted (shallow over a near-black bottom) is now handled by C2.
  - A full version needs absolute library reflectances, but photic has only shapes normalised near 550 nm.
- **B9, water-mass priors: deferred to Phase D.** Priors on P, G and X must come from real water, and on tier 3's synthetic water types they would only encode the synthetic values. Fit them on the Exmouth training sites.

### B10 result (2026-10-13, `LEE_*` macros in samodel.c, `b10_fit.py`)

Lee's coefficients were refitted to tier 2's OSOAA cases: 732 clear and moderate cases for training, 366 turbid for testing, physical albedo only. Subsurface rrs relative error:

| Coefficients (g0 g1 c0 c1 b0 b1) | train median / p95 | turbid (held out) median / p95 / bias |
|---|---|---|
| Lee: 0.084 0.170 1.03 2.4 1.04 5.4 | 1.91% / 5.94% | 3.12% / 7.78% / +1.84% |
| free fit: 0.083 0.149 2.95 7.2 1.01 8.3 (c0, c1 at their bounds) | 1.16% / 4.69% | 2.67% / 6.97% / −0.23% |
| rrs_dp and D_uB only: 0.083 0.147 1.03 2.4 1.024 5.25 | 1.55% / 5.33% | 2.83% / 7.48% / −0.46% |
| all within ×0.5–2: 0.083 0.147 2.06 4.8 1.015 7.30 | 1.30% / 5.02% | 2.81% / 7.27% / −0.40% |

Lee's bias by κH (+1.6 to +2.0% at κH 1–8) is removed by the refits. **The effect on depth is small.** On Phase A OSOAA pixels, noise-free:
- moderate 5–10 m: bias +2.76 → +2.58 m (constrained fit) or +2.32 m (×0.5–2 fit);
- clear 10–30 m: bias +3.4 to +3.8 m, almost unchanged.

The free fit pushes the water-column elongation (c0) to about 3, which is not physical. **Decision:** not adopted. Lee's coefficients stay the default; the refits stay available through the macros. The remaining depth bias is not mainly from Lee's coefficients.

### Phase B order (from Phase A)

1. **B1. Bottom parameterisation and starts.** One albedo per region, bounded to [0, 1], plus n_b − 1 mixing fractions, replacing B × q. Start albedo mid-range, with several starts. Bounds on P, G and X replace `fabs()`.
2. **B2. K prior:** remove it, or replace it with weak bounded priors on P, G and X.
3. **B3. Y:** fit Y per scene, or use QAA's ratio estimate. Add the tier 2 variant first.
4. **B4. Depth coherence:** penalise departure from a local plane rather than from the mean, or lower its weight.
5. **B5. Optimiser:** after B1 to B4, compare Nelder–Mead with a bounded or global method on the reduced problem; drop the LUT shortcut if it costs accuracy.
