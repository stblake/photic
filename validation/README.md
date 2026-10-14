# photic physics validation benchmark

A repeatable, scored benchmark for photic's semi-analytical (SA) model. Every change to the model core should be judged against it. It has four tiers, and each isolates one source of error:

| Tier | Question | How | Runtime |
|---|---|---|---|
| 1 | Does the C code implement the equations? | `drivers/forward_dump` (photic's `samodel_Rrs`) against `photic_ref.py` (written from the paper), plus SAMBUCA (CSIRO) and the published IOP tables | seconds |
| 2 | How good is Lee's semi-analytical approximation? | photic's forward model against exact radiative transfer (OSOAA: scalar, no atmosphere, sun glint removed with a black-ocean run) for 1,260 cases: depth, water type, bottom, sun angle and band | ~30 min the first time, then cached |
| 3 | Can the inversion retrieve the truth? | Synthetic 32×32 multi-scene images built from OSOAA lookup tables (with noise, tides, an optional depth prior), plus a control generated with photic's own forward model | OSOAA tables ~5 min per scene (cached), then photic only |
| 4 | Does it work on real imagery? | Exmouth Gulf: six Sentinel-2 L1C scenes → ACOLITE → photic, against WA DoT singlebeam surveys (and the published photic product) at four sites | ACOLITE ~5 min per scene (cached), then photic only |

## Running

```sh
cd validation
PY=~/photic_dev/validation_data/.venv/bin/python
$PY bench.py run --tier 1 --tier 2 --tier 3 --tier 4   # add --quick for the day-to-day subset of tier 3
$PY bench.py score                                      # writes scorecards/<commit>.json and .md
$PY bench.py compare baseline <commit>                  # what changed
```

The benchmark names photic's spectral data in every command file it writes (`spectra_inputs.py`), so the baselines do not move when photic's defaults do:
- `PHOTIC_BENCH_BOTTOM`: the bottom library (default `model/data/bottom_hydrolight.csv`).
- `PHOTIC_BENCH_WATER`: `legacy` (default; the 10 nm Pope & Fry, Morel and Lee tables), `default` (photic's own), or a list of water CSVs.
- `PHOTIC_BENCH_SRF=1`: tier 4 scenes get their sensor's spectral response.

Scorecards record the setting.

`bench.py run --model Lee_Profile --cid <name>` runs tiers 3 and 4 with the profile model (`MODEL Lee_Profile`, EXPERIMENTS.md Phase F) instead of `Lee_Semi_Analytical`.

`make` builds `drivers/forward_dump` against `../model`. Results for each commit go to `$PHOTIC_VALIDATION_DATA/runs/<commit>/` (default `~/photic_dev/validation_data`); only scorecards are kept in git.

## Quick A/B tests (`ab.py`)

To iterate on the model, compare two variants on a few hundred out-of-sample points in about a minute, before running the full benchmark:

```sh
$PY ab.py HEAD .                  # the working tree against the last commit
$PY ab.py . .:-DLUT=0             # a compile-time switch, without editing the code
$PY ab.py . .:-DY_MODEL=2 --split train --seed 2
$PY ab.py HEAD .~Lee_Profile      # the profile model (MODEL Lee_Profile) against photic's SA model
```

- **Variants:** a git ref or `.` (the working tree), optionally followed by `:` and extra compiler flags. Each build is cached under `ab/build/`, and its photic must support `MODEL ... POINTS`.
- **Points:** a fixed, depth-stratified sample (`--n 300`, `--seed`). `--split test` (the default) takes them from tier 3's turbid water and tier 4's Bundegi and Tantabiddi; `--split train` from the training data.
- **How it runs:** both variants run on the full scenes with `POINTS`, so only the sampled cells are modelled.
- **Not tested:** behaviour that depends on modelling neighbouring cells in sequence, such as the lookup table and carry-over. Those need the full benchmark.
- **Report:** paired per tier:
  - coverage;
  - RMSE and the RMSE difference with a 95% bootstrap interval;
  - share within tolerance and wins;
  - RMSE by depth;
  - σ calibration and run time.
- **Verdict:** a change counts as better or worse only if the interval excludes zero and the RMSE moves by at least `--min-effect` (2%).

This is a screen, not the acceptance test below. Use it to choose what to benchmark. Reusing one seed for many decisions overfits that sample, so confirm a promising change with another seed or split.

## Diagnosing the inversion

`EXPERIMENTS.md` is the experiment plan and log for making the model work without a depth prior. `profile_model.py` is the Python prototype of the profile model (Phase F), run on `ab.py`'s points against a photic variant. `paper_profile.py` runs the profile-model paper's experiments (ablation, scene subsets) and writes its figures (`paper/msab_figures/`). `phaseA.py` runs photic's own error function and optimiser on single tier 3 pixels through `drivers/invert_lab`, against the known truth.

## Data and tools (all under `$PHOTIC_VALIDATION_DATA`, not in git)

The steps are idempotent: anything already present is reused.

- `.venv`: Python 3.14 with numpy, scipy, xarray, netCDF4, rasterio, pyproj, GDAL 3.13 (JPEG-2000), sambuca_core, pyTMD and pypdf.
- `tools/OSOAA` (github.com/CNES/RadiativeTransferCode-OSOAA, v2.0), built with gfortran, in two patched copies:
  - `tools/OSOAA_scalar`: `CTE_POLAR_SWITCHED_OFF 1`, `EXPERT_MODE_FORCED_FATM_NULL .TRUE.`;
  - `tools/OSOAA_blackocean`: also `EXPERT_MODE_FORCED_FSEA_NULL .TRUE.` (surface reflection only, used to remove sun glint).
- `tools/acolite` (github.com/acolite/acolite).
- `raw/tables/pope97.txt`: Pope & Fry (1997), from omlc.org.
- `raw/lee/Lee_2_QAA_HOPE_hyper.xls`: Lee's own QAA and HOPE spreadsheet (supplied by the user; not redistributed). Tier 1 checks that `photic_ref` with the spreadsheet's coefficients reproduces its modelled Rrs, and reports the effect of photic's coefficients.
- `raw/s2/*.SAFE`: Sentinel-2 L1C tile 49KHR, downloaded from the public `gcp-public-data-sentinel-2` bucket (`tier4_fetch_s2.py`).
- `raw/l8/`: Landsat 8 OLI L1TP (Collection 1, path 115 row 75), six scenes from the public `gcp-public-data-landsat` bucket (`tier4_fetch_l8.py`). `tier4_prepare.py acolite_l8` then `sites_l8` builds `derived/tier4/sites_l8` (~31 m cells); run tier 4 on it with `PHOTIC_TIER4_SITES=sites_l8`.
- `raw/bottom/measured`: in situ underwater reflectance libraries from PANGAEA (CC BY 3.0): Heron Reef 804589, Moreton Bay 864310, Fiji 872506, Cook Islands 872505, Glovers Reef 824861 (`bottom_library.py fetch`). `bottom_library.py pack` puts every spectrum, on 1 nm bins with its class and half, into `model/data/bottom_spectra_insitu.nc` (in git, 0.6 MB). `bottom_library.py build` then writes `model/data/bottom_measured.csv` from the library half of each class; tier 3's `_insitu` scenarios draw their bottoms from the held-out half.
- `raw/bathymetry/exmouth`: WA DoT Exmouth point surveys (`EX.zip`; singlebeam; x/y/z in GDA94 MGA zones 49/50; metres relative to AHD; CC BY 4.0) and the 5 m marina composite.
- `raw/tides`: WA DoT Exmouth 02 tide gauge, 1997–2025, 5-minute data, cm above LAT, AWST (CC BY 4.0).
- `raw/published`: the photic Exmouth Gulf 30 m product linked from the README.
- `raw/tides/got`: the GOT4.10c global ocean tide model (NASA GSFC, open), read with pyTMD, for sites without a gauge. At the Exmouth gauge (2019, hourly) it has r = 0.96 but underestimates the range (sd 0.34 against 0.53 m; residual 0.2 m RMS).

## Acceptance rule for changes to the model core

1. Work on a branch.
2. Run all four tiers, `bench.py score`, then `bench.py compare baseline <commit>`.
3. Tier 1's equations must still pass. If the change alters the forward model on purpose, update `photic_ref.py` to match the new intended equations.
4. Accept the change only if the depth metrics in tiers 3 and 4 (RMSE, bias and `frac_within_tol` overall and by depth bin) improve with no regression elsewhere.
5. Then update `model/tests/baseline` (`make test-baseline`) and `scorecards/baseline.json` in the same commit, explaining why the numbers moved.
