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

### Phase B order (from Phase A)

1. **B1. Bottom parameterisation and starts.** One albedo per region, bounded to [0, 1], plus n_b − 1 mixing fractions, replacing B × q. Start albedo mid-range, with several starts. Bounds on P, G and X replace `fabs()`.
2. **B2. K prior:** remove it, or replace it with weak bounded priors on P, G and X.
3. **B3. Y:** fit Y per scene, or use QAA's ratio estimate. Add the tier 2 variant first.
4. **B4. Depth coherence:** penalise departure from a local plane rather than from the mean, or lower its weight.
5. **B5. Optimiser:** after B1 to B4, compare Nelder–Mead with a bounded or global method on the reduced problem; drop the LUT shortcut if it costs accuracy.
