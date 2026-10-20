# Future experiments: a relocatable, learned, physics-constrained bathymetry model

A plan, not yet started (written 2026-10-20).

**The goal** is a model trained once, on ICESat-2 depths and soundings from many sites, that can then be moved to a new site and run on imagery alone, with no local soundings, ICESat-2 or calibration of any kind. At every pixel it should give:
- depth, with σ_H;
- the water's attenuation κ_s(λ);
- the bottom reflectance ρ(λ).

It must improve on Lyzenga's log-linear model, which has to be refitted at every site.

**How this fits the rules.** `EXPERIMENTS.md` allows soundings to tune settings that apply everywhere, never as a per-pixel input at run time. A network trained across sites is such a setting, so this plan is a candidate route to the sounding-free goal itself. It is scored against DROP, the current sounding-free model. Every experiment follows the same rules: a hypothesis, a prediction written before it runs, a metric, and a recorded result, including failures.

## Why Lyzenga cannot be the answer, and what to keep from it

At the Arabian Gulf palm box, joint Lyzenga with `POWER` (six scenes × four bands), fitted to ICESat-2 there, scores **90.0% hits** and RMSE 1.21 m on the test half. DROP scores 19.4% hits, with no correct depth beyond 5 m. That gap is the size of the prize. But Lyzenga's strength is local:
- **Its coefficients are the site's water and bottom.** The slope is 1/(κM), and the band weights cancel the site's bottom-brightness changes. Moved to another water, another atmosphere or another sensor, they are wrong. H0(c) measures by how much.
- **It saturates near 16 m and reads 0–2 m too deep.** Regression dilution compresses the depths, and `POWER` is a one-parameter patch for it.
- **It needs every scene at a cell, is ill-conditioned** (25 coefficients, |R| ratio 91), **and gives depth only.**

**What to keep:**
- Lyzenga's log features;
- its joint use of all scenes, which is worth 10 points over a composite;
- the evidence that labelled depths fix what imagery alone could not: the attenuation's absolute scale.

## What relocation requires

1. **Inputs that mean the same at every site.**
   - No coordinates and no site identity.
   - Reflectance is expressed against each scene's own deep water (rrs∞_s, by DROP's deep-pixel rule) and noise (σ_n, from the deep pixels).
   - The scene's tide comes from a global tide model (GOT4.10, or FES), which is not a sounding.
   - Sun and view geometry and the sensor's band set are inputs.
2. **Physics in the loop.** A learned mapping from reflectance to depth will meet waters and atmospheres outside its training set. A forward model inside the network lets it check itself and adapt at a new site:
   - fit the scene's water and atmospheric offset to that scene's own reflectance at test time (H5);
   - flag, as out of distribution, cells it cannot reconstruct.
3. **Diversity in training:**
   - many sites, stratified by water clarity, bottom type, latitude and sun, aerosol and haze, tide range and sensor;
   - synthetic scenes from photic's tier 3 machinery for pretraining, where κ and ρ are known.
4. **Abstain where the bottom is not visible,** as DROP does, instead of returning the training set's mean depth.

## The physics Lyzenga linearises

With DROP's forward model (Lee's, EXPERIMENTS.md Phase F), for scene s, band λ and water column Z_s = H + h_s (tide):

  rrs_sλ = rrs∞_sλ [1 − exp(−κ_sλ M^C_s Z_s)] + exp(−κ_sλ M^B_s Z_s) ρ(λ)/π

With M^C ≈ M^B = M, this gives

  ln(rrs_sλ − rrs∞_sλ) = ln(ρ(λ)/π − rrs∞_sλ) − κ_sλ M_s Z_s

**The crux for relocation is κ's absolute scale.**
- Within one site's imagery, κ and H trade off against each other. EXPERIMENTS.md G1: neither the per-pixel misfit, the pooled misfit nor the tides fix the scale on these scenes.
- QAA estimates it from rrs∞, but ACOLITE's red/NIR residual inflates it, which is DROP's depth ceiling.
- A model trained on labelled depths at many sites can learn how κ's scale relates to what the image shows: rrs∞'s spectrum, the shape of the shallow-water spectra, the NIR/SWIR residual, the geometry. That is a learned replacement for QAA, including ACOLITE's typical residual.
- **Whether that relation transfers to sites it has not seen is the central hypothesis,** tested by H4 and H7. If it fails, no amount of network capacity rescues relocation, and the model falls back to needing local labels.

## Training and test data: many sites

**A site pipeline.** Generalise the gulf tooling (`tier4_fetch_atl24.py`, `tier4_gulf.py screen/select/acolite/reference`, `gulf_lyzenga.py prepare`) into a site-agnostic `sites.py add <name> <box>`. For each site it:
- fetches ATL24;
- screens and selects scenes by haze, then sun, then tide;
- runs ACOLITE on the regular-grid output used for the whole-scene gulf;
- builds the reference with the cross-pass consistency filter;
- prepares the bands;
- deletes the L1 data afterwards, leaving about 1–3 GB per site.

**Sites.** Choose them before any model is trained, stratified across:
- **water clarity:** a Kd(490) climatology, from oceanic to turbid coastal;
- **bottom:** sand, seagrass, coral and mixed;
- **latitude and sun elevation;**
- **aerosol:** desert dust such as the Gulf, maritime, continental;
- **tide range;**
- **sensors:** Landsat 8/9 first, Sentinel-2 later.

Record each site's strata in a table in this file. The first aim is 15–30 training sites; H7 decides how many are enough.

**Truth beyond ICESat-2.** Sites with surveyed soundings or bathymetric lidar give truth at depths ICESat-2 under-samples, and full maps rather than tracks:
- Exmouth (WA DoT);
- candidates to verify for availability and licence: national hydrographic and lidar archives such as AusSeabed and NOAA's coastal lidar.

Reserve the best-surveyed sites as test regions.

**Held-out regions, fixed now and never used for any choice:**
- **Exmouth:** soundings. It has been used to tune DROP, but never by this model, so it stays unseen for it.
- **The Arabian Gulf whole scene:** ICESat-2, including the palm box. It is excluded from training entirely.
- **At least two more regions** from different strata, chosen with the training sites.

Score each held-out region once per frozen model, with the prediction written here first.

**Synthetic pretraining.**
- Scale up tier 3: random IOPs (P, G, X, Y), bottoms drawn from the measured library, depths of 0–30 m, geometries, scene counts and tides, red/NIR offsets δ of 0–5 × 10⁻³ (the G0 plan), and noise.
- Use OSOAA's lookup tables, so the training set carries the real structural error of Lee's model (~2.4%) rather than matching the decoder exactly.
- κ, ρ and H are known at every pixel.

## The model ladder

Each rung must beat the one below it in **leave-region-out cross-validation** over the training sites before the next is built. The ladder separates what learning across sites buys from what the physics buys.

**H0. Baselines, on the same cells and splits.**
- (a) **Local Lyzenga:** joint `POWER`, fitted to each site's own labels. This is the upper reference, "locally calibrated". In the new harness it must reproduce `cv.json`'s 83.9% at the palm box to within 0.5 points.
- (b) **DROP default:** the relocatable reference to beat.
- (c) **Relocated Lyzenga:** coefficients pooled over the other sites and applied unchanged. This measures how badly a local empirical model relocates.

**H1. Relocatable regression on normalised features, depth only.**
- **Features:**
  - ln(rrs − rrs∞_s) and ln rrs, each normalised by σ_n;
  - log band ratios;
  - rrs∞_s's own spectrum and the NIR/SWIR residual;
  - tide, sun and view zenith.
- **Models:** gradient-boosted trees (LightGBM) and a small MLP.
- **Loss:** Huber on depth, weighted to balance the depth bins and the sites.
- **Hypothesis:** the site-relative features carry enough of the water's state for depth to transfer.
- **Prediction:** beats H0(c) and DROP in leave-region-out hits, but stays well below H0(a).
- **Watch:** regression dilution (report slope and r) and the trees' inability to extrapolate past the deepest labels.

**H2. A set encoder over scenes.**
- **Architecture:** a shared per-scene network, taking each scene's features with its tide, geometry and band set, then masked pooling over the scenes present (Deep Sets).
- **Why:** it takes any number of scenes, any mix of sensors and missing scenes, as a relocatable model must.
- **Prediction:** coverage near 100%, with hits falling smoothly as scenes are dropped.

**H3. Spatial context.**
- **Architecture:** a small CNN over patches of about 7×7 to 15×15 cells, generalising `RADIUS`.
- **Hypothesis:** neighbours constrain the water and the bottom, and reduce noise.
- **Splits:** with leave-region-out there is no label leakage between train and test. Within a site, masking labels near the scored tracks is needed only for H0(a).

**H4. Physics decoder: depth, attenuation and bottom together (the core experiment).**
- **Encoder:** H3's encoder outputs a latent state per cell.
  - **Depth:** H below MSL, with Z_s = H + h_s.
  - **Bottom:** mixture weights over the library's classes, plus a brightness bounded by class (B8). ρ(λ) is a bounded combination of measured shapes, not a free spectrum.
  - **Water:** for each scene, κ_s(λ) from IOPs in Lee's parameterisation as photic uses it: a = a_w + a_φ(P) + a_g(G), b_b = b_bw + X(440/λ)^Y. This fixes κ's spectral shape physically and keeps κ ≥ a_w + b_bw. P, G and X are a smooth field (coarse grid, upsampled), so the water cannot bend to the bottom (the failures of decoupling step 3 and of B12).
  - **Deep water:** rrs∞ comes from the same IOPs (Lee: (g₀ + g₁u)u), so the deep water and the attenuation agree.
  - **Offset:** a per-scene offset δ_s(λ), unattenuated and smooth in space, with its spectral shape penalised (PROFILE_QAA mode 3, Phase G).
- **Decoder:** the forward model, differentiable (PyTorch or JAX). Before training, it must match `photic_ref.py` and photic's DROP equation to 1e-6 in rrs (the tier 1 discipline).
- **Loss:**
  - **Reconstruction:** χ² of every scene's rrs at every water cell, with σ² = σ_n² + (ε rrs)². This is self-supervised, and it is the only loss available at a new site.
  - **Depth:** Huber on labelled cells.
  - **Priors:** bounds, smoothness of the water and of δ, and one bottom shared across scenes.
  - **Pretraining:** on the synthetic scenes first, with κ and ρ supervised too; then fine-tune on the real sites.
- **Hypotheses:**
  - (a) training through the physics makes depth transfer better than H3;
  - (b) κ and ρ are recovered;
  - (c) the encoder predicts κ's absolute scale at unseen sites.
- **Predictions:**
  - (a) leave-region-out hits beat H3 by ≥2 points, mostly beyond 10 m;
  - (b) on held-out synthetic scenes, κ within 10% and ρ brightness within 15% where the bottom is visible;
  - (c) at held-out sites, the encoder's κ(490) is within 20% of the κ that the site's own labels imply (H0(a)'s slope).
- **Ablations, one at a time:**
  - no reconstruction loss (H4 → H3);
  - no synthetic pretraining;
  - free per-pixel κ (expect degeneracy);
  - no δ;
  - one bottom class;
  - M^C = M^B.

**H5. Test-time adaptation at a new site, with no labels.**
- At the new site, keep the encoder's depth and bottom as a prior.
- Refit each scene's water field and δ_s to that site's reflectance by the reconstruction χ² alone, then re-infer H and ρ.
- Optionally fine-tune the encoder's last layers on the same loss.
- **Hypothesis:** adaptation corrects the water and atmosphere where they differ from training, without giving back κ's scale. The prior on the scale is what the labels taught; the χ² alone cannot move it (G1).
- **Prediction:** gains at the held-out sites least like the training set.
- **Check:** the scale must not drift; compare κ before and after.

**H6. Uncertainty and abstention.**
- **σ_H:** from a deep ensemble (5 seeds) combined with the decoder's sensitivity ∂rrs/∂H.
- **Optically deep:** where the sensitivity is below the noise, the cell is optically deep (DROP's T_deep), and the model abstains.
- **Out of distribution:** a reduced reconstruction χ² well above 1, or ensemble disagreement, flags the cell and lowers its confidence.
- **Target** (EXPERIMENTS.md Phase G): 60–75% of covered cells within ±σ_H, at held-out sites.

**H7. Learning curve: how many sites are enough.**
- Train on 1, 2, 4, 8, 16, … sites, as random subsets stratified by water type with several draws, and score leave-region-out each time.
- If held-out hits are still rising at the largest set, add sites. If they are flat, the model is data-sufficient, and the remaining gap is the model or the inputs.
- Report the curve by stratum: a stratum with few sites (e.g. turbid) may need its own.

## Validating κ and ρ without field truth

1. **Synthetic scenes** (held-out IOPs, bottoms and geometries): score depth, κ and ρ against the truth at every cell, including OSOAA's structural error and the red offsets. **Gate:** κ and ρ are reported for real scenes only if they pass here.
2. **Leave one scene out.** Fit with five scenes, then predict the sixth's reflectance from the latent H and ρ, refitting only that scene's water and δ from its own image. Score by the reduced χ², against DROP. This tests the whole latent state on unseen data.
3. **Tides.** At fixed H and ρ, the reflectance change between scenes of different tide is set by κM. Compare the predicted and observed slopes, binned by depth.
4. **κ against the labels.** At each held-out site, the κ the site's own depths imply (H0(a)) is a measured check that the model never saw (H4c).
5. **Independent water estimates:** over deep water, compare with ocean-colour Kd(490) on the same day (Sentinel-3 OLCI or MODIS; the G1a plan).
6. **Bottom plausibility:** ρ inside the measured library's envelope (`bottom_spectra_insitu.nc`), and class maps compared qualitatively with a public habitat map.

## Splits and held-out discipline

- **The primary split is by region:** leave-region-out cross-validation over the training sites.
  - Splits within a site, by pass as in `gulf_lyzenga.py`, serve only the local baseline H0(a).
  - Sites whose boxes or scenes overlap go in the same fold.
- **Held-out regions** (above) are used for nothing until a frozen model is scored on them, once. A second look needs a new region.
- **The reference's own errors:** apply the cross-pass consistency filter at every site, and use a robust depth loss. ICESat-2 thins out where attenuation limits it, so the surveyed test sites carry the deep evaluation.
- **Bootstrap over regions** (and over passes within them), not over cells, for every interval.

## Metrics

- **Depth at held-out regions:**
  - hits (within 0.5 m + 10%) over all reference cells, the headline;
  - RMSE, bias, r and slope;
  - all by depth bin (0–2, 2–5, 5–10, 10–20, 20–30 m), with coverage.
- **Two gaps per region:**
  - over DROP: the relocatable improvement;
  - to local Lyzenga, H0(a): what calibrating locally would still add.
- **σ_H:** the share of covered cells within ±σ_H, by depth bin.
- **κ and ρ:**
  - synthetic: κ relative error by band; ρ brightness error and spectral angle;
  - real sites: H4(c)'s κ ratio, and the leave-one-scene-out χ².

## Gates

- **GR1 (H1–H3):** in leave-region-out cross-validation, the best model beats both DROP and relocated Lyzenga (H0(c)), with the bootstrap 95% interval for hits above 0.
- **GR2 (H4–H5):**
  - synthetic κ and ρ meet H4(b);
  - leave-region-out depth at least matches H3;
  - H4(c) holds at a majority of the cross-validation regions;
  - if (c) fails, relocation of κ's scale is not supported. Report that, and narrow the goal (e.g. one model per water type, or a few labels per site).
- **GR3 (held out, once per region):** the frozen model, with no local labels:
  - beats DROP's default at every held-out region with no material RMSE regression;
  - closes at least half the gap to local Lyzenga.
  - **The target:** match local Lyzenga's hits at 0–10 m, the Phase G target restated for a relocated model.
- **GR4 (H7):** the learning curve is flat at the training-set size used. If it is not, add sites before claiming a final result.

## Implementation

- **Prototype in Python** under `validation/`:
  - `sites.py`, the site pipeline;
  - `ml_depth.py`, the models and the harness, reading the same prepared bands, references and splits as `gulf_lyzenga.py`, so H0 is a like-for-like check.
- **Dependencies:** add PyTorch (or JAX) and LightGBM to the venv on the new machine, and record them in the README.
- **Scale and storage:**
  - about 30 M water cells per Landsat scene area, × 6 scenes × 4 bands;
  - train on random patches (10⁵–10⁶ per epoch) and infer tile by tile;
  - storage is about 1–3 GB per site after deleting L1, so 30 sites take tens of GB;
  - a GPU is needed from H3 on.
- **Determinism:** fixed seeds and deterministic kernels; record seeds in the log.
- **Into photic, only after GR3:**
  - a `MODEL` that reads the network's weights and the bottom library from files (data in files, not C arrays);
  - the encoder's inference, plus the H5 adaptation, which is a small least-squares fit of each scene's water and δ through photic's existing forward model;
  - a synthetic regression test in `make test`;
  - CNN inference in C is the hardest part. Keep a variant with a per-cell MLP and `RADIUS` features as the first deliverable.

## Risks

- **κ's scale may not be predictable from imagery.** If it depends on things the image does not show (e.g. vertical structure, or a sensor- or region-specific atmospheric residual), relocation fails at the scale, as for DROP. H4(c) and H7 test this early, before the C port.
- **Atmospheric correction differs by region.**
  - ACOLITE's residual over desert dust is not that over maritime aerosol. A model can learn one region's residual as "water".
  - The explicit δ_s(λ), test-time adaptation (H5) and aerosol strata in the training set are the defences.
  - A G6-style joint atmosphere fit would remove the dependence on ACOLITE altogether.
- **Training-set bias.** ICESat-2 sees clear water best, so turbid sites are under-represented. Stratify, and report per stratum.
- **Few deep labels:** tracks give thousands of cells per site and few beyond 15 m. The reconstruction loss, synthetic pretraining and the surveyed sites carry the deep range.
- **Degeneracy in H4:** κ, H and ρ brightness can trade off. Coarse water fields, bounded classes and the synthetic gate are the checks. If needed, fix P, G and X from the deep pixels (DROP's approach) and learn only a correction.
- **Leakage between sites:** overlapping boxes, shared scenes or the same reef system in two folds. Group them.

## Relation to the rest of the project

- **This is now the project's main line** (decision of 2026-10-20). DROP is frozen at its default settings and no longer developed. It remains in photic for two reasons:
  - it is the relocatable baseline (H0(b)) that every gate must beat;
  - its forward equation is H4's decoder, and photic's DROP code is the reference the decoder is checked against.
- **Phase G's findings are design inputs here,** not open work:
  - the atmospheric offset δ_s(λ);
  - κ's absolute scale as the cause of the depth ceiling;
  - aerosol varying across the image.
- **It replaces the open decision G1b(iii)** (ICESat-2 as DROP's κ-scale anchor at each site). Labels are used once, globally, in training, never at the site modelled, which is within `EXPERIMENTS.md`'s rules.
- **Fallback if GR2 finds that κ's scale relocates but the network does not:** feed the learned κ-scale predictor into DROP's exact depth profile. That would be the only reason to change DROP again.
- **The whole-scene gulf run** (EXPERIMENTS.md, "In progress") is the first site through the new pipeline, and is reserved as a held-out region.

## Order and rough cost

1. **Site pipeline and site selection, written into this file before any training.** The whole-scene gulf is the first test of the pipeline. About a week, then download and ACOLITE time per site.
2. **H0 and H1,** leave-region-out on the first 5–10 sites. A few days, on CPU.
3. **Synthetic pretraining set, and the differentiable forward model checked against `photic_ref.py`.** About a week.
4. **H2 and H3.** About a week.
5. **H4,** synthetic then real, **and H5.** Two to three weeks. GR2 decides whether relocation of κ's scale is supported.
6. **H6, then H7** (more sites as needed).
7. **GR3 on the held-out regions, once.** Then the C port and a paper.
