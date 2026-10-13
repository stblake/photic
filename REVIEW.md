# photic code review (October 2026)

**Scope:** every C source file in `model/`, about 23k lines of project code. The vendored `linenoise.c` and `asa047.c` were reviewed only where photic calls them. The review is read-only: no source files were changed.

**Method:**
- Line-by-line reading.
- Compilation with Apple clang 17 and GCC 16 at `-Wall -Wextra`.
- Small standalone test drivers for the I/O and numeric helpers, some run under AddressSanitizer.
- Physics checked against `paper/photic.tex`.

**Status update:** the build items (priority 1–2 and §1) the physics items (priority 3–5) and the data-corruption items (priority 6–8, 11–12) the memory-safety items (priority 9–10, plus most of §3 critical/high) and the ALIGN items (priority 13; §5 #5–7, #9a) were fixed after this review; items marked **[fixed]** are resolved, see the git history. The code now builds with Apple clang 17 and GCC 16.

**Labels:** severity is critical / high / medium / low. **[verified]** means the finding was reproduced by a test program or by compiling. Items marked *uncertain* depend on intent that the code doesn't record.

---

## 0. Priority list

Fix these first; each is either a build blocker or silently corrupts scientific output.

| # | Area | Issue |
|---|---|---|
| 1 | build | `OK`/`COMMENT`/`EXIT`/`STOP` are undefined, so `bam.c` does not compile at HEAD. They were deleted in `c0038fe`; the original values were 1000/1001/1002/1005. [verified] |
| 2 | build | Globals are defined in `common.h`, which needs `-fcommon`; GCC ≥ 10 gives 122 link errors. `typedef int bool` fails under C23 (GCC ≥ 15). The makefile has MacPorts paths hard-coded. [verified] |
| 3 | physics | **[fixed]** `HTIDE` has no effect: `samodel_Rrs` uses `H`, not `H + tide` (samodel.c:2873, 2931, 2934). Multi-temporal fits ignore tidal differences. [verified] |
| 4 | physics | **[fixed]** Jerlov interpolation weights are swapped (jerlov.c:385-386), and the water-type fraction is reversed for decreasing ratio tables (jerlov.c:333). [verified] |
| 5 | physics | **[fixed]** The QAA coefficient `h2 = 0.469` should be `-0.469` (secchi.c:165). [verified] |
| 6 | I/O | **[fixed]** `compress_2d` starts from `grmax = FLT_MIN` and uses `else if`, so the true maximum can be missed and wraps in the 16-bit pack. A constant grid gives scale 0, then NaN. This affects every NetCDF write (nc.c:286-316). [verified: 10 → 1.99979] |
| 7 | I/O | **[fixed]** Bilinear interpolation is wrong on the last row and column; bicubic never reads the last row or column (interp.c:73, 124). This affects UNITE, INTERPOLATE, PANSHARPEN, OVERLAY and SRTM. [verified] |
| 8 | memory | **[fixed]** `smooth()` reads `[j-1]` at `j=0` (smoothutils.c:96-121). [verified: ASan heap-buffer-overflow] |
| 9 | memory | **[fixed]** `tokenize()` has no bounds on `MAX_ARGS` or `MAX_ARG_SIZE` (bam.c:6597), so it overflows the global `parsed`. |
| 10 | memory | **[fixed]** Empty names match unused slots in every grid and scene lookup. Missing arguments therefore cause NULL dereferences, use-after-free and double free (e.g. `DELETE` with no argument). |
| 11 | process | **[fixed]** `PROCESS ADD`: the inverted `!= 0` tests mean the default output is all zeros and `COEFFICIENTS` is ignored (process.c:488-502). [verified] |
| 12 | process | **[fixed]** `PROCESS UNITE`: `\|\| method_median_no_spval` tests the constant 6 (always true), so nodata is mixed into MEAN/MIN/EXP and the MEDIAN path writes `scene_list[-1]` (process.c:2359, 2379). [verified] |
| 13 | align | **[fixed]** `ALIGN` resamples soundings with replacement by default, which drops ~37% of them; `MAXSOUNDINGS > n` overflows (align.c:252-266). |
| 14 | OpenMP | `samodel()` is one big `omp parallel` block, so setup, the σ pass, `depth *= -1` and the file writes all run on every thread. It is latent only because OpenMP is off. |

---

## 1. Build and infrastructure

- **critical** `bam.c:287,397,399,6327-6332,6588-6594`: the macros are undefined (see the priority list). Restore them in `bam.h`.
- **high** `common.h:316-460`: globals are defined in the header. Move them to `extern` declarations, with a single definition in `common.c`.
- **high** `common.h:58`: `typedef int bool` is invalid in C23. Use `<stdbool.h>`, or add `-std=gnu17` to the makefile.
- **medium** `makefile`: paths are hard-coded to `/opt/local`, so Homebrew NetCDF is not found. There are no `-Wall` flags, no rules for the standalone utilities, `clean` uses `rm` without `-f`, and the targets have no dependencies, so they always rebuild.
- **medium** `filllakes.c:10`: includes a non-existent `filllakes.h`.
- **low** None of `io.h`, `nc.h`, `interp.h`, `smoothutils.h`, `filllakesutils.h`, `srtm.h` or `process.h` has an include guard. `MEMCHECK` is defined twice (`common.h`, `io.h`). `smooth.h` declares an `inline verify` that has no definition.
- **low** `usgs2grid.sh`: EPSG:32639 is hard-coded although the comment says "pseudo mercator"; `$xyz2nc` is never set; `cut -c1-43` assumes a fixed filename length; `mv *.nc` moves every NetCDF file in the directory.
- **low** About 140 unused or set-but-unused variables.

Standalone utilities (they build with `common.c nc.c io.c` plus a helper, and `-fcommon`):

| Tool | Status |
|---|---|
| xyz2nc | builds (make rule) |
| asc2nc, grid2grid, smooth, smoothspval, emptygrid, srtm | build, but have no make rule |
| filllakes | needs the include fixed |

---

## 2. Core semi-analytical model (`samodel.c`)

### Critical
1. **OpenMP region design** (L305-1758). Every thread runs all of the following:
   - setup and sorting;
   - allocation of the 16000-entry LUT, about 400 MB per thread;
   - `srand`, plus `rand()` in the depth-σ pass (not thread-safe);
   - writes to the shared `depth_sigma`;
   - `depth *= -1` (L1488), so with an even number of threads the sign comes out wrong;
   - all the NetCDF writes.

   In addition, `schedule(dynamic)` defeats the sorted-order hot start. Restructure the code so that only the pixel loop is parallel before ever enabling `-fopenmp`.
2. **[fixed]** **The tide is ignored** in the forward model (L2867-2873, L2931, L2934). `md->depth = H + tide` is dead code: it is overwritten in `samodel_optimise` (L2000-2017).

### High
3. `md->h_empirical` is uninitialised when no DEPTHS grid is given, and stale when the grid is nodata at a pixel (L1826, L1138, L968-971). It controls the "sand only if > 8 m" decision.
4. `best` is used uninitialised if every `nelmin` result is NaN or ≥ 1e4 (L2385, L2413-2416). NaN can be reached when Rrs440 = Rrs490 = 0, through `chi = 0/0` (L2898).
5. `md->K` beyond `n_bands[k_scene]` is uninitialised `malloc` memory, but it is scanned by `array_min_double2(md->K, n_scenes, max_n_bands, 0)` (L1994, L2705) and copied into the LUT.
6. The output `model_error` grid stores `md->Rrs_error` rather than `md->model_error`. For LUT hits, the value comes from a different pixel (L1121 vs L1072).

### Medium
7. `bottom_albedo` and K are taken from the *last* region evaluated, not from `origin` (L2877-2881, L2915).
8. The high-K penalty `K_min > 0.7` is evaluated outside the scene loop, so only the last scene counts (L2728).
9. In the depth-σ pass, `trial_depths[k_sample]` is left unset when `n_regions == 0` (L1430). The pass also hot-starts across unrelated random pixels (L1454).
10. Precedence bugs: `(int) max_depth_reached/depth_interval` casts before dividing, so pixels deeper than that floor (or than 30 m) get σ = 0 (L1385, L1387).
11. `NBOTTOMS` is not validated:
    - above 8 it overflows `bottom_types` (8 rows);
    - above 12 it overflows `B[12]`, `q[12]` and `bottom_type_indexes[12]`;
    - 0 gives `q_sum = 0`, so the results are NaN.

    Separately, `NSMOOTH = 0` silently skips every pixel.
12. The band order coastal/blue/green/red is hard-coded as indices 0-3. With fewer than 4 bands the code reads out of bounds (L360-363, L718-739, L1153, L2803).
13. Unmodelled pixels output depth 0 instead of nodata. The `bottom_*` defaults are −9999, but the files are written with the input band's nodata value. `extract_Rrs_data` does not reject negative reflectance in neighbouring pixels, even though the origin pixel check does (L939 vs L3003).
14. There are no parameter bounds; the paper claims user-defined ranges. B is unbounded, so `1 - 1.5*rrs` can reach ≤ 0. κ is clamped to [0, 2.5] (L2913), which is not in the paper.
15. LUT entries leak their inner arrays (L1718). `free_model_data` (L3143) is unused, and would be wrong for array elements.
16. Discrepancies between the code and the paper. Decide which is authoritative:

    | Item | Code | Paper |
    |---|---|---|
    | χ | `Rrs440/Rrs490` (L2898) | (R440−R750)/(R490−R750), which is commented out at L2897 |
    | P start | 490/550 | 440/550 |
    | B start | 4.0·mean Rrs490 | 0.4·Rrs(490) |
    | Δ | disabled (`DELTA 0`) | included |
    | Error weights | 80/15/10/15, with extra B and K terms | 0.85/0.15 |
    | H start depths | 8 | 17 |
    | LUT | 16000 entries, matched by mean band difference | 256 entries, matched by SAM |
    | Region size | (2n−1)² (L2971, L2981), but allocation uses (2n+1)² | (2r+1)² |
    | a₀/a₁ | Lee's coefficients | CSIRO coefficients from Appendix A |

    The alternative `a0_Blake`/`a1_Blake` coefficient arrays have 36 entries instead of 41, so enabling them would read out of bounds.

### Low
17. `PLOTTING=1` does not compile: there is a stray `}` at L1364. [verified] `LUT=0` reads uninitialised variables (L1008, L1104), and `ALL_BOTTOMS=0` leaves `best_*` unset.
18. Possible divisions by zero or NaN:
    - in `model_spectral_diff_mean` (L2833);
    - in the IOD calculation (L2058);
    - `acos > 1` in the sort index (L762), which feeds NaN into `qsort`;
    - `q_sum = 0` (L2035, L2495);
    - a zero simplex step when a start value is 0.
19. The K-penalty thresholds for H < 2 and H < 4 are identical (L2713 vs L2717). *uncertain*
20. Dead code:
    - `samodel_display` (wrong parameter layout);
    - `polcyn_depth` (declared, never defined);
    - `spectral_weights`, `bottom_coral_dekker`;
    - `if (true || …)` (L2547);
    - an unreachable `k_shallow > n_shallow_pre` check (L808).
21. Type mixing: `approx_equal` (float) is used on doubles; `vec_mean2_double` returns float; progress timing uses `clock()`.
22. `samodelgraphics.c` sizes its plots from `n_bands[0]` for every scene, and the x-axis is fixed at 400-700 nm.

## 2b. Findings from the synthetic regression test (`make test`)

These were found after the original review. The test generates noise-free scenes from known depths with photic's own forward model and tables, so the true parameters reproduce the data exactly. Instrumenting `samodel_error` confirmed that the error at the true parameters is 0.

1. **critical — the model does not recover depth from its own noise-free forward model.**
   - The test uses two scenes, H 1.5–8.25 m, a sand bottom with B 0.6, and turbid IOPs (K_min 0.30/0.35 m⁻¹).
   - The optimiser converges to local minima with about 15% Rrs misfit, low albedo (B ≈ 0.06) and depths far too deep (bias +50 m, RMSE 52 m).
   - Simplifying to a single-pixel region (`NSPATIAL 1`) and a single bottom (`NBOTTOMS 1`) does not help.
   - Likely contributors: Nelder–Mead starting points (8 H starts from 40 m down to 1 m; B started at 4·mean Rrs490); `reqmin`/`konvge`; and no parameter bounds (positivity only via `fabs`).
   - *To investigate before trusting results in new areas.*
2. **high — a hard-coded turbid-water prior in `samodel_error` (samodel.c ~L2693-2718).**
   - For H < 5 m, a K_min below about 0.19–0.28 m⁻¹ adds a penalty of 100·(1/(0.01+K) − 1/(0.01+K₀))².
   - At K = 0.05 m⁻¹ (clear water) that is about 1.4×10⁴, which dwarfs the Rrs term (in percent). So clear, shallow water is pushed to deep solutions.
   - This prior is not in the paper. It suits the turbid Pilbara/Gulf sites, but biases clear-water sites.
3. **high — the solution is very sensitive to the error-term weights.**
   - Changing `weight_depth` from 15 to 15.5 (+3%) moves the modelled depths by up to 12.6 m and changes P/G/X in almost every pixel.

## 2c. Findings from the physics validation benchmark (`validation/`, October 2026)

See `validation/README.md` and `validation/scorecards/baseline.json`.

**Tier 1: the equations are implemented exactly.**
- photic's `samodel_Rrs` matches an independent Python implementation of the paper's equations bit for bit over 51,840 cases (relative difference 0).
- It matches CSIRO's SAMBUCA to 2e-15 given the same water backscatter. With SAMBUCA's own b_bw formula the difference is ≤ 0.7%.

Deliberate departures from the paper, as the relative change in Rrs if switched to the paper's form:

| Departure | Median | 95th percentile | Max |
|---|---|---|---|
| **a₀/a₁:** the code uses Lee's tables, but the paper says the CSIRO-derived Appendix A tables are used | 0.6% | 13% | **72%** |
| **χ:** the code omits the paper's 750 nm correction | 0.04% | 1.8% | 2.3% |

The Y and κ clamps are never active on the grid. The tables are within 1–3% of their sources (a_w vs Pope & Fry, limited by photic's 10 nm table; b_bw vs Lee's formula to 0.1%).

**Tier 2: Lee's semi-analytical model is a good approximation to exact radiative transfer.** This is against OSOAA (scalar, no atmosphere, sun glint removed with a black-ocean run), over 1,098 cases in the physical range: 10 depths, 3 water types, 3 bottoms, 3 sun angles and 5 Sentinel-2 bands.
- Subsurface rrs: median 2.4%, 95th percentile 6.7%, max 11%.
- Above-surface Rrs: median 2.4%, 95th percentile 7.0%, max 11.6%.
- photic is systematically about 3% high where κH ≈ 1–4, and the bias depends on sun angle (+2% at 20°, −1% at 60°).

New finding:

1. **high: bottom albedo can exceed 1.** photic's seagrass and coral spectra are normalised to 1 at 550 nm but reach about 3.3 and 2.9 at 705 nm, so B·ρ_table exceeds 1 (it is about 2.0 for B = 0.6). That is physically impossible, and there Lee's Rrs conversion overshoots by up to 42%. B is unbounded (§2 #14), so the inversion can reach such states. These 162 cases are excluded from the tier 2 statistics and reported separately.

**Tier 3: the inversion does not retrieve depth on its own.** Synthetic 32×32 scenes are built from OSOAA (lookup-table interpolation error ≤ 0.15%).
- **Without a depth prior, photic fails in every scenario:** clear, moderate and turbid water; 1, 2 or 4 scenes; with or without noise. Bias is +6 to +77 m and at most 15% of pixels are within ±(0.5 m + 10%). This includes the control generated with photic's own forward model, so it is optimiser and heuristic failure, not model-structure error.
- **With a ±20% depth prior (`DEPTHS`), photic returns roughly the prior.** Bias is about 0 and RMSE 2–3 m, but its errors correlate 0.8–0.99 with the prior's own errors, so it adds little or no skill (`rmse_skill_vs_prior` −0.68 to +0.21; negative is worse than the prior). The scorecard now reports the prior's own score and `rmse_skill_vs_prior`.
- **Retrieved albedo is about 0.05 instead of 0.6 in every scenario.**
- **The run-time LUT shortcut copies results between pixels.** photic copies a previous LUT entry whenever the spectra match within an adaptive threshold and the depth and K are within 10% of the previous pixel's; that rule is in §2 of samodel.c. In the quick scenario about 100% of pixels were copied rather than optimised.

**Tier 4: Exmouth Gulf against WA DoT singlebeam surveys.** Six Sentinel-2 scenes (2019–2024) were processed with ACOLITE, with tides from the Exmouth 02 gauge, at four sites. Scenes covering less than half as much of a site as the best scene are dropped. Run without a depth prior, the current code fails at every site (baseline scorecard, `validation/scorecards/baseline.json`):

| Site | photic bias | RMSE | within tol. | r | published bias | RMSE | within tol. | r |
|---|---|---|---|---|---|---|---|---|
| Marina | +31.2 m | 41.4 m | 5% | 0.31 | −0.11 m | 0.97 m | 68% | 0.96 |
| Bundegi | +43.7 m | 133 m | 2% | 0.23 | +0.64 m | 3.7 m | 42% | 0.78 |
| Tantabiddi | +19.4 m | 24.9 m | 5% | −0.17 | not covered | | | |
| Gulf entrance | +29.0 m | 30.2 m | 0% | 0.05 | +1.94 m | 2.2 m | 50% | 0.96 |

- The build-fix-only commit `a0eb70f` gives the same result (bias +29.3 m vs +29.0 m), so this is not a regression from the fixes in this review.
- The **published photic Exmouth Gulf product** (README) scores well on the same cells (right-hand columns); it does not cover Tantabiddi.
- So either the published runs used an empirical depth prior and settings this benchmark does not reproduce, or the code that produced them differs from the code in this repository. Establishing which is the first task for improving the core.

---

## 3. Interpreter (`bam.c`)

### Critical
1. **[fixed]** `tokenize` (L6597-6613): no bounds checks. A run of about 130 spaces, or a token of 128 or more characters, overflows `parsed`.
2. **[fixed]** `DELETE` (L5556-5605): the pointer is not set to NULL after `free`, and an empty name matches a free slot. The result is a double free, and `n_grids_in_use` can go negative.
3. **[fixed]** Empty-name lookups (`parse_old_grid`, `parse_old_scene`, and inline lookups at L835, 1962, 2788, 2856, 4030, 5556, 5656, 6119, …) resolve to unallocated slots. All lookups should also require `allocated_grids[n]`.

### High
4. **[fixed]** `COMPUTE K … COORDINATES` frees the `bpl_*` buffers even though they are only allocated in REGION mode (L2587-2591).
5. `COMPUTE LSM DEEP` uses the uninitialised `lon`/`lat` instead of `londeep`/`latdeep`, so DEEP never works (L2714-2727).
6. **[fixed]** `READ COMMANDS`: the `while(!feof)` loop does not check `fgets`, so the last line runs twice (L6268-6278). The same bug is in `srtm.c:99`.
7. **[fixed]** `READ GRID`: the name is registered before the file-existence check, which leaves a poisoned slot (L6474-6489).
8. **[fixed]** `READ SOUNDINGS`: the old soundings are freed before the file is checked, leaving a dangling pointer (L6365-6375).
9. **[fixed]** `free(mparams)` does not set the pointer to NULL, and the polcyn model never allocates it. Running lyzenga, then polcyn, then any model gives a double free or use-after-free (L3977-3980, L4652, L4947, L5160-5214).
10. **[fixed]** `MODEL stumpf_log_ratio` is accepted, but `zdepth` stays NULL, so it crashes (L4444-4467).
11. **[fixed]** Orphaned "allocated" slots with NULL arrays are left after early returns (L1700/1715/1754, 3043-3086, 4071-4423, 4147/4185). Exit cleanup then crashes.
12. **[fixed]** Loops bounded by `MAX_GRIDS` (265) index `parsed[k]` (max 128) (L1345, L4024).
13. Unchecked user-coordinate indexing in `COMPUTE K` line mode (L2156-2187) and `COMPUTE LSM REGION` (L2660-2686).
14. `COMPUTE LSB` continues with uninitialised `land_index`/`shallow_index` when they are missing (L2759-2801).
15. `COMBINE`: `R_blue` is never set (L786 is commented out). `best_scene = -1` then indexes `scene_data[-1]` (L790-799).

### Medium
16. **[fixed]** A NULL return from `linenoise` on EOF (Ctrl-D) leads to `strlen(NULL)`. The `> MAX_LINE` check is off by one, and `raw_line` leaks (L344-354).
17. `COMBINE` resets `Rinf_mean`/`nobs` inside the scene loop, and averages `ps_in[0]` instead of `ps_in[m]` (L580-598).
18. `model_error()` divides the global `weight_local_constraints` by `n_model_bands` on every call, so the weight decays towards 0 (L5355).
19. `COMPUTE PERCENTILE 100` reads `vec[k]`, one past the end (L1922).
20. `COMPUTE K REGION` leaves the tail of `Ri`/`Rj` uninitialised (L2203-2295).
21. The empirical models use stale global `land_index`/`shallow_index` when LAND/SHALLOW are omitted (L4079-4423).
22. `BOUNDARY` sampling loops forever (`while(true)`) if no valid pixels exist (L4335-4375).
23. The "grid table full" guards leave `ps` uninitialised. L2992 checks room for +8 grids but allocates 10.
24. `COMPARE`: the nodata check happens before `fabs` in one place and after it in another. `NEGATESOUNDINGS` is a no-op (L957-980).
25. `DELETE` does not update the scene `band_indexes` or the global slot indices, so they become dangling.
26. `SCENE` errors leave a half-built scene, and there is no scene DELETE.
27. `strcat`/`sprintf` into the fixed label buffers can overflow (L2480-2490, L2569).

### Low
- `fopen` and `cpgopen` failures are unchecked, and several buffers leak.
- Several NaN results come from divisions by zero (L1094, 1116, 4696, 5304).
- `READ COMMANDS` can recurse without limit.
- `EXIT` and `STOP` are prefix matches.
- `SET MODEL OPTIMISATION` succeeds but also prints "unknown command".
- `PRINT LINEWIDTH` prints the wrong value.
- `HELP` is an empty stub.
- `random_search()` is empty, and `model_optimisation` is never read.
- Several functions are declared but never defined.
- The argument loop in `main` does nothing.
- Usage text disagrees with the parser in several places.

---

## 4. Grid processing (`process.c`)

### Critical
1. **[fixed]** `ADD`: the inverted tests (L488-502) are described in the priority list. [verified]
2. **[fixed]** `UNITE`: the `|| method_median_no_spval` bug (L2359, L2379) is described in the priority list. [verified]

### High
3. ~~`INTERPOLATE REGION <grid>`: `ps_region` is declared `bool`, so the extent always comes from slot 0 or 1 (L3144, L3226).~~ **Not a bug** while `bool` is `typedef int`, but it *would* become one under `_Bool`/C23. That is one reason `bool` must stay an `int`.
4. `POLYGONS`: a missing polygon file causes an infinite loop, because `continue` runs without `k++` (L998-1001).
5. **[fixed]** `UNITE`: the default extent starts at `elon = 0`, so grids west of Greenwich produce an enormous output (L2112, L2280).
6. `lyzenga_extinction_cutoff`: the float loop runs 101 times, so `partial_indx > k` and reads out of bounds; it also divides by zero when k = 0 (L3738-3742). This is reached by `SHALLOW THRESHOLD AUTOMATIC`. [verified]
7. `OVERLAY`: bicubic interpolation extrapolates outside the layer, so a small layer overwrites the whole base grid (L1304-1318).
8. `SHALLOW GRIDS`: the bands are not checked for matching dimensions or for nodata (L3655-3683).
9. `MNDWI` is actually NDWI: it uses NIR, where MNDWI uses SWIR. The automatic threshold is a radiance value compared against an index in [−1, 1], the NIR band's nodata is not checked, and it can divide by zero (L1667-1720).
10. **[fixed by the per-command sweep]** Grids are registered before validation in ADD, RELATIVEERROR, ERROR, OVERLAY, MASK and SHALLOW. CLOUD, NOISE, DEGLINT and PANSHARPEN leave stale names.

### Medium
11. `DEGLINT`: divides by zero when nobs = 0 or var = 0, and NIR nodata is not checked. It subtracts relative to mean(NIR), where Hedley et al. (2005) use min(NIR). *uncertain*
12. `SHALLOW`: pixels above max (land, cloud) are counted as deep water, it divides by zero when nobs = 0, and GRIDS is only recognised at `parsed[9]`.
13. `SMOOTH WEIGHT` is read with `atoi`, so fractional weights are truncated (L3413).
14. `RELATIVEERROR` multiplies by \|d\| when \|d\| < 1, and returns 0 when \|d\| < 0.1 (L695-700). *uncertain*
15. `MEDIAN` with an all-nodata window reads `vec[0]` uninitialised and divides by zero (L221-244).
16. `POLYGONS`: a mismatch between the FILL and polygon counts uses uninitialised fill values. EXTERIOR mode ignores RANGE outside the bounding box.
17. `LAND`/`MNDWI`: the `coords_present` test is inverted (L1440, L1654).
18. `extent_eps = 0.005°`, about 550 m, is also used for the cellsize tolerance, so grids offset by many pixels pass the commensurability check.

### Low
- `FILL` has no guard for a full grid table.
- Empty arguments match empty slots.
- `UNITE` with no inputs, or RESOLUTION ≤ 0, is not rejected.
- `MEDIAN` skips one extra edge row and column.
- `MASK`/`CLOUD` RADIUS semantics: r = 1 is a single pixel, and CLOUD's RADIUS is never used.
- `PANSHARPEN` can divide by zero, and `fabs` discards negative values.
- `lyzenga_extinction_cutoff` has dead code and prints 100 diagnostic lines per call.
- Pixel mapping uses `ncols` where it should use `ncols-1`.
- There is heavy duplication: about 17 copies of the lookup boilerplate and 7 of the commensurability check. `process.h` is included twice, and the header comments are stale.

---

## 5. Optics, alignment, refine, plotting

### High
1. **[fixed]** `jerlov.c:385-386`: the interpolation weights are swapped. [verified]
2. **[fixed]** `jerlov.c:333`: the water-type fraction is reversed for decreasing ratio tables, which is the normal case for Landsat-8 band pairs. Example: the code gives 0.098 where about 0.90 is expected.
3. `jerlov.c:28-49`: the table mixes Austin & Petzold (1984) values for types I–C1 with Jerlov (1976) values for C3–C9. It is not monotonic at C1→C3, so a single ratio can match two intervals.
4. **[fixed]** `secchi.c:165`: `h2` should be `-0.469`. [verified]
5. **[fixed]** `align.c:252-266`: the default resampling with replacement; overflow when MAXSOUNDINGS > n; division by zero when n = 0.
6. **[fixed]** `align.c:396-398`: nodata `z_model` entries are transformed, then missed by the exact-equality nodata check at :627.
7. **[fixed]** `align.c:316-318`: `NEGATESOUNDINGS` flips the sign once per grid, so it is a no-op for an even number of grids.
8. `refine.c:242-244`: specifying only LAND or only SHALLOW makes the whole output nodata. The mask dimensions are not checked.
9a. **[fixed]** `align.c:177-210` (found while fixing #5): the soundings-file branch reused the argument index `k` as a loop variable, so every option after `SOUNDINGS <file>` (CHARTDATUM, SPVAL, MAXSOUNDINGS, NEGATEGRID, NEGATESOUNDINGS) was silently ignored. These options are now read. `ALIGN` also now rejects an empty input-grid list.
9. `plot.c:191-258`: `ii` is reused as both coordinate index and window index, so multi-line LOGSCATTER plots read the wrong coordinates.

### Medium
10. In plot and write, row and column mapping uses `nrows` where it should use `nrows-1`, with no bounds check (plot.c:229-241, 529-538; write.c:199-213).
11. A vertical line gives `grad = inf`, so NaN (plot.c:214, 519; write.c:189).
12. `REFINE SCALE` with `SHAPE` ≠ 1 does not map to [dmin, dmax] (refine.c:264-268).
13. `REFINE` and `ALIGN` register their outputs before validation.
14. `ALIGN` weights:
    - the most common depth gets weight 0;
    - a single sounding gives NaN;
    - a sounding of 0 divides by zero;
    - the kernel width is hard-coded at 1 m.
15. `ALIGN`: `weighted_params` can be left uninitialised, and `ifault` is not checked.
16. `ALIGN`: `CHARTDATUM` and `SPVAL` are parsed but ignored. The code minimises mean absolute relative error ×100, whereas the paper describes a weighted RMS. The paper's `exp(+…)` for W_i looks like a typo.
17. `secchi.c`: uses max(Rrs) instead of Rrs at the band of minimum Kd. The Kd(530) weights sum to 0.95. Wavelengths are 481/554/656 where Landsat-8 uses 482/561/655. *uncertain*
18. `secchi.c`: no guards against zero, negative or NaN Rrs, so `sqrt`/`log10` can be given negative arguments.
19. Histograms with an explicit region write out of bounds for points outside that region (graphics.c:86-89, 247-251).
20. The GMT palette defines only colour indices 17-36, but up to index 80 is used (graphics.c:107-141).
21. The `L > Lsm + 1.0` threshold is a magic number that depends on input units (jerlov.c:107). *uncertain*

### Low
- Palette index off-by-ones: the log histogram reaches 81, and RYGB reaches 83.
- A constant grid divides by zero in `plot_array`.
- The REGION clip branches are identical, so the clip does nothing.
- Image cells are drawn half a cell off.
- The NOAA land colour is set but never used.
- The title `_` replacement mutates `parsed`.
- `compute_k` has no bounds check.
- `compute_k_from_ratio` is unused.
- The header and definition array sizes disagree (`-Warray-parameter`).
- `jerlov_water_type_str` ignores `str_len`.
- Several leaks on `cpgopen` failure.
- The axis-label `strcat` can overflow.
- `log(0)` occurs at Lsm in `plot.c`.

---

## 6. Utilities and I/O (`common.c`, `io.c`, `nc.c`, `interp.c`, `write.c`, `smoothutils.c`, `filllakesutils.c`, standalone tools)

### High
1. **[fixed]** `nc.c:286-296`: the `compress_2d` maximum bug. [verified]
2. **[fixed]** `nc.c:305,316`: a constant or all-nodata grid gives scale 0, so 0/0 → NaN, which is then converted to short (undefined behaviour).
3. **[fixed]** `smoothutils.c:96-121`: an out-of-bounds read at the edges. [verified, ASan]
4. `io.c:66,147`: the header token buffer `str[32]` is read with an unbounded `fscanf("%s")`.
5. `io.c:229-237`: `read_asc` does not check counts. Extra values write to a negative row, too few leave the grid uninitialised, and a non-numeric token loops forever.
6. **[fixed]** `interp.c:73-74`: bilinear interpolation on the last row or column returns column n−2. There is no bounds check, and `nx = 1` reads index −1. [verified]
7. **[fixed]** `interp.c:124-127`: bicubic `msize = ny-1`, so the last row and column are never used. [verified: 100 → 0]
8. `write.c:199-213`: `WRITE LINE` reads out of bounds.
9. `write.c:285-296`: `WRITE GRID … REGION` writes nothing.

### Medium
10. `write.c:270,306`: `.asc` output is silently written as NetCDF.
11. `nc.c:219-229`: `read_nc` uses untyped `nc_get_var`, so a file with double coordinates or float data overflows the buffer. The variable type, dimension order and latitude direction are never checked. It only reads files written by photic.
12. `nc.c:206,209`: `MEMCHECK(lons)` checks the pointer-to-pointer, so the check never fires.
13. `io.c:262`: `write_asc` uses `%.2f`, which destroys reflectance and K values (0.0123 → 0.01).
14. Georeferencing: `xllcorner` is used as the coordinate of the first cell centre, so NetCDF coordinates are half a cell south-west of the true centres. `xyz2nc` uses centres, so the tools are inconsistent with each other.
15. Comparing a float cell with a double nodata value fails for values that are not exactly representable, e.g. −9999.9 (`smoothspval.c:57,86`, `srtm.c:127,131`). [verified]
16. `filllakesutils.c:111,113`: the east and north distance weights ignore the box offset.
17. `xyz2nc.c:85-97`: the spacing-detection loop has no bound. Cell placement uses `floor()`, collisions overwrite earlier points, and dx ≠ dy is not handled.
18. `io.c:285,299,352,440`: `char c = fgetc()` is an EOF bug where `char` is unsigned (ARM Linux). `read_csv` adds a phantom row when the file ends with a newline.
19. `fopen` is never checked in `io.c` or `srtm.c`. Every tool copies `argv` values into 512-byte buffers with unbounded `strcpy`.
20. **[fixed]** `common.c:1449-1458`: when the table is full, `parse_new_grid` and `parse_new_scene` leave `*ps_out` unset and then `strcpy` through it.
21. `common.c:307-330`: `interp_1d` with descending X uses uninitialised values when x is out of range.

### Low
- `linear_fit` r is NaN for constant input.
- Standard deviations use the population form (divide by n).
- `random_in_range(min == max)` divides by zero.
- `random_norm` can take `log(0)`.
- `compare_ints` can overflow.
- Functions return inconsistent values on empty input (`BIG`, 0, or nodata).
- `vec_mean2_double` returns float.
- `malloc` is unchecked in several places.
- `smoothed_index_2d` uses a (2r−1) window.
- `srtm` processes the last file twice (`feof` loop).
- `write_xyz` always appends.

`point_in_polygon` is correct (a standard crossing test). [verified]

---

## 7. Cross-cutting recommendations

1. **Make it build anywhere.** Restore the macros. Make the globals `extern` with definitions in `common.c`. Use `<stdbool.h>`. Use `nc-config` and pkg-config in the makefile, add `-Wall -std=gnu17`, and add rules for the standalone tools.
2. **Make the registry safe.** Write a single `find_grid(name)` that rejects empty names and checks `allocated_grids`, and a single `register_grid(name, …)` that is called only once the array exists. Set pointers to NULL on free. Use these instead of the ~30 inline copies.
3. **Write a nodata policy.** Write one `is_nodata(v, spval)` helper that handles NaN, and use it everywhere.
4. **Fix the physics bugs** in the priority list (tide, Jerlov, QAA), then reconcile each `samodel.c` vs paper difference (§2.16): decide whether the code or the paper is authoritative and document it.
5. **Write a regression harness.** A small synthetic scene (forward-model Rrs from known H, B and IOPs → `MODEL Lee_Semi_Analytical` → compare) would catch most physics regressions. Add unit drivers for `interp`, `nc` round-trip, `smooth` and `read_asc` under ASan.
6. **OpenMP:** move setup, the σ pass and the writes out of the parallel region. Use per-thread RNG state (`rand_r` or similar) and a static schedule over pixel chunks, which keeps the hot-start ordering.
