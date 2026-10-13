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

`make` builds `drivers/forward_dump` against `../model`. Results for each commit go to `$PHOTIC_VALIDATION_DATA/runs/<commit>/` (default `~/photic_dev/validation_data`); only scorecards are kept in git.

## Diagnosing the inversion

`EXPERIMENTS.md` is the experiment plan and log for making the model work without a depth prior. `phaseA.py` runs photic's own error function and optimiser on single tier 3 pixels through `drivers/invert_lab`, against the known truth.

## Data and tools (all under `$PHOTIC_VALIDATION_DATA`, not in git)

The steps are idempotent: anything already present is reused.

- `.venv`: Python 3.14 with numpy, scipy, xarray, netCDF4, rasterio, pyproj, GDAL 3.13 (JPEG-2000), sambuca_core, pyTMD and pypdf.
- `tools/OSOAA` (github.com/CNES/RadiativeTransferCode-OSOAA, v2.0), built with gfortran, in two patched copies:
  - `tools/OSOAA_scalar`: `CTE_POLAR_SWITCHED_OFF 1`, `EXPERT_MODE_FORCED_FATM_NULL .TRUE.`;
  - `tools/OSOAA_blackocean`: also `EXPERT_MODE_FORCED_FSEA_NULL .TRUE.` (surface reflection only, used to remove sun glint).
- `tools/acolite` (github.com/acolite/acolite).
- `raw/tables/pope97.txt`: Pope & Fry (1997), from omlc.org.
- `raw/s2/*.SAFE`: Sentinel-2 L1C tile 49KHR, downloaded from the public `gcp-public-data-sentinel-2` bucket (`tier4_fetch_s2.py`).
- `raw/bathymetry/exmouth`: WA DoT Exmouth point surveys (`EX.zip`; singlebeam; x/y/z in GDA94 MGA zones 49/50; metres relative to AHD; CC BY 4.0) and the 5 m marina composite.
- `raw/tides`: WA DoT Exmouth 02 tide gauge, 1997–2025, 5-minute data, cm above LAT, AWST (CC BY 4.0).
- `raw/published`: the photic Exmouth Gulf 30 m product linked from the README.

## Acceptance rule for changes to the model core

1. Work on a branch.
2. Run all four tiers, `bench.py score`, then `bench.py compare baseline <commit>`.
3. Tier 1's equations must still pass. If the change alters the forward model on purpose, update `photic_ref.py` to match the new intended equations.
4. Accept the change only if the depth metrics in tiers 3 and 4 (RMSE, bias and `frac_within_tol` overall and by depth bin) improve with no regression elsewhere.
5. Then update `model/tests/baseline` (`make test-baseline`) and `scorecards/baseline.json` in the same commit, explaining why the numbers moved.
