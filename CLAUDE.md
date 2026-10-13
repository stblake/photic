# CLAUDE.md

Guidance for working on **photic**, a physics-based, semi-analytical, shallow-water satellite-derived bathymetry (SDB) model written in C. It adapts Lee et al.'s HOPE model (1998/1999) from hyperspectral to multispectral imagery (Landsat 8 / Sentinel-2, atmospherically corrected with ACOLITE). It does this by assuming that depth and bottom composition are constant across multiple scenes (up to tide) and that water optical properties are constant across a small spatial region around each pixel. The method is described in `paper/photic.tex` (arXiv:2002.02298), which is the reference for the intended physics.

## Repository layout

- `model/`: all source code (~26k lines of C, flat directory, single `makefile`).
- `paper/`: LaTeX source and figures for the paper (`paper/makefile` runs pdflatex).
- `*.jpg` at the root: preview images of published SDB datasets, referenced from `README.md`.
- There are no tests, no CI, and no example command files or sample data in the repo.
- Commit history is entirely "Add files via upload" (GitHub web uploads), so git history explains nothing; diff between commits to see what changed.

## Building

```sh
cd model
make                  # photic and xyz2nc
make utils            # asc2nc grid2grid smooth smoothspval emptygrid srtm filllakes
make WARN=-Wall       # with warnings
make clean
```

Build notes:

- **Compilers:** builds warning-free with Apple clang 17. With Homebrew GCC 16 (`make CC=gcc-16`), the only warning is a `-Wstringop-overflow` in `align.c:437`. It is real: `ALIGN` does not reject zero input grids.
- **Dependencies:**
  - NetCDF flags come from `nc-config` (Homebrew here).
  - PGPLOT is expected under MacPorts `/opt/local`. Override with `make PGPLOT_DIR=...`, e.g. for giza. It links `-lcpgplot -lpgplot -lX11`.
- **Language standard:** must be compiled as **C17 or earlier**. The makefile passes `-std=gnu17`. `common.h` uses `typedef int bool` and `#error`s under C23. Do not switch to `<stdbool.h>`/`_Bool`: some code stores non-0/1 values (e.g. grid slot indices) in `bool` variables.
- **Globals:** `common.h` *declares* the globals (`extern`) and `common.c` *defines* them, with their initial values. When adding a global, add both.
- **Header dependencies:** the makefile builds object files with `-MMD -MP`, so header changes trigger the right rebuilds. Build artefacts are ignored via `model/.gitignore`.
- **OpenMP is not enabled:** there is no `-fopenmp`, and `#include <omp.h>` is commented out. Do not enable it without restructuring `samodel()` (see Gotchas).
- **Parser return codes:** `OK`/`COMMENT`/`EXIT`/`STOP`, defined near the top of `bam.c`, are what `parse_input()` returns.

## Running

`photic` is an interactive interpreter (REPL via linenoise) over a registry of named in-memory grids.

- `./photic` starts the REPL; `./photic cmds.bam` runs a command file, then drops into the REPL unless the file ends with `EXIT`.
- Command syntax:
  - Tokens are UPPERCASE, case-sensitive and positional, split on **single spaces** (double spaces create empty tokens; tabs are not separators).
  - Optional keywords must appear in the order each handler expects.
  - `//` starts a comment line. `STOP` ends the current command file. `EXIT` quits immediately (no cleanup). `%` and `%1`..`%9` recall previous lines.
- Files are written to the current directory as side effects: `history.bam`, `k.bam`, `model.bam`, `output.csv`, and the SA-model outputs `modelled_*.nc` and `<scene>_K_*.nc`.

Main commands (dispatch in `bam.c:run_command`):

| Command | Where | Purpose |
|---|---|---|
| `READ GRID <name> <file>` / `READ SOUNDINGS` / `READ BOTTOM` / `READ COMMANDS` | bam.c | Load rasters (.nc/.asc), xyz soundings, bottom spectra CSV, command files |
| `WRITE GRID\|LINE\|POINTS` | write.c | Output grids (always NetCDF), transects, point samples |
| `PROCESS <op> IN ... OUT ...` | process.c | Grid ops: INTERPOLATE, SMOOTH, SHALLOW, PANSHARPEN, NOISE, DEGLINT, UNITE (mosaic), CLOUD, MASK, LAND, OVERLAY, POLYGONS, ERROR, RELATIVEERROR, MNDWI, FILL, ADD, MEDIAN |
| `SCENE NAME s BANDS ... WAVELENGTHS ... THETAVIEW v THETASUN w [HTIDE h] ...` | bam.c | Define a multi-band scene for the SA model |
| `MODEL Lee_Semi_Analytical IN <scenes> OUT <10 grids> [DEPTHS g]` | bam.c → samodel.c | **The core physics-based model** |
| `MODEL lyzenga_log_linear\|polcyn_log_ratio ...` | bam.c | Older empirical models fitted to soundings (simulated annealing) |
| `MODEL Lee_Kd_LS8\|Lee_Secchi_LS8\|Kd_490_Mueller\|Chl_Morel ...` | bam.c, secchi.c | Band-ratio water-quality products |
| `COMPUTE AREA\|LSM\|LSB\|BOTTOM\|K\|PERCENTILE` | bam.c, jerlov.c | Deep-water signal, attenuation (Jerlov), areas, etc. |
| `COMPARE GRID g SOUNDINGS f` | bam.c | Validation statistics vs soundings |
| `ALIGN`, `REFINE` | align.c, refine.c | Blend/calibrate depth grids against soundings; clip, scale and mask |
| `PLOT GRID\|LINE\|LOGSCATTER` | plot.c, graphics.c | PGPLOT graphics |
| `SET <key> <val>`, `PRINT <key>`, `DELETE <grid>`, `COPY`, `COMBINE` | bam.c | Globals, inspection, registry management |

## Architecture

- **`common.h`/`common.c`**: constants (`MAX_GRIDS` 265, `MAX_SCENES` 32, `MAX_ARGS`/`MAX_ARG_SIZE` 128, `BAM_NODATA` −99999), the structs `geogrid`, `scene` and `model_data`, all global state, and helpers (statistics, allocation, `approx_equal`, grid and scene lookup through `parse_old_grid`/`parse_new_grid`/`parse_*_scene`, `commensurate_grids`).
- **Grid registry**: parallel global arrays `gridded_data[]`, `grid_names[]` and `allocated_grids[]`, plus `n_grids_in_use`, with per-slot parameter arrays (`Lsm`, `Lsb`, `Alpha`, `thetas`, …). Handlers read arguments from the global `parsed[][]` and look grids up by linear `strcmp` over names.
- **`bam.c`**: `main`, REPL, command-file runner, tokenizer, dispatcher, and most `run_*` handlers (SCENE, MODEL, COMPUTE, COMPARE, SET/PRINT, empirical models). `#define BAM_VERSION "1.95"`.
- **`samodel.c`**: the semi-analytical model.
  - `samodel()` sets up IOP spectra, sorts pixels into a processing order, loops over pixels with a 16000-entry LUT, then estimates depth σ by re-modelling with noise and writes NetCDF outputs.
  - `samodel_optimise()` and `samodel_optimise_one_bottom_combination()` call Nelder–Mead (`nelmin` in `asa047.c`) on `samodel_error()`.
  - `samodel_Rrs()` is the Lee forward model.
  - Compile-time switches sit at the top of the file: `LUT=1 ALL_BOTTOMS=1 DELTA=0 SAM=0 RESTART=0 LOCALLY_AVERAGE=1`. Only these defaults build cleanly; `PLOTTING=1` does not compile.
  - The parameter vector is `[H (N_r) | B (N_r·N_b) | q (N_r·N_b) | P,G,X per scene]`, with B, P, G and X scaled ×100. Positivity is enforced only through `fabs()`, and there are no bounds.
  - Arrays are indexed `[region][scene][band]`, padded to `max_n_bands`, and the padding is uninitialised.
- **`jerlov.c`**: Jerlov water-type attenuation table and interpolation. **`secchi.c`**: Lee 2016 Landsat-8 Kd/Secchi (QAA-style).
- **`io.c`** reads and writes ESRI `.asc`, `.xyz` and `.csv` files, choosing the format by substring `.nc`/`.asc` in the filename. **`nc.c`** handles photic's own NetCDF layout: dimensions `lons`/`lats`, variable `bathy(lats,lons)` stored as **16-bit packed shorts** (lossy). It does not read arbitrary NetCDF files.
- **`interp.c`**: bilinear and bicubic interpolation, plus `grid2gridinterp`. **`smoothutils.c`**, **`filllakesutils.c`**: filters.
- Vendored third-party code; don't modify casually: `linenoise.c/.h` (antirez) and `asa047.c/.h` (Burkardt's Nelder–Mead).

## Conventions

- **Grid orientation**: `array[row][col]`; row 0 is **south**, column 0 is **west** (`.asc` is flipped on read and write). Grids are node-registered: `lat = slat + i*cellsize`, `lon = wlon + j*cellsize`, and `elon = wlon + (ncols-1)*cellsize`. Coordinates are geographic degrees (WGS84).
- **Argument order** differs between helpers: `allocate/copy_float_array_2d(…, nrows, ncols)`, but `smooth`, `fill_lakes` and the interpolators take `(ncols, nrows)`.
- **Depth sign**: bathymetry outputs are **negative down**. Empirical-depth inputs to the SA model are also expected negative.
- **Band order**: the SA model and Secchi code assume scene bands 0–3 are coastal, blue, green, red.
- **"Radius" r** in smoothing and spatial-region code means a (2r−1)-wide window, so r=1 is a single pixel.
- **Nodata comparisons** are inconsistent: some code uses exact `==`, some `approx_equal` (a *relative* tolerance, so nodata 0 needs exact equality), and some treats values ≤ 0 as invalid. NaN is never handled. When editing a function, match the convention it already uses.
- **Adding a command handler**: follow the existing `run_*` pattern of reading `parsed[k]` and looking grids up. When creating an output grid, set the name, `allocated_grids[ps]=true` and `n_grids_in_use++` only **after** validation succeeds and the array is allocated. Many existing handlers get this wrong.
- Style: 2-space indentation in newer files, tabs in older ones (`samodel.c`, `common.c`). Comments are `//` section headers in sentence case ending with a full stop. Australian/British spelling (`normalise`, `colour`).

## Gotchas / known hazards

- **Empty names match unused slots.** Unused grid and scene slots have name `""`, so a missing argument resolves to an unallocated slot (NULL array, use-after-free after `DELETE`). Lookups should also check `allocated_grids[n]`.
- **`tokenize()` has no bounds checks** on `MAX_ARGS`/`MAX_ARG_SIZE`.
- **The OpenMP design in `samodel()` is unsafe.** The whole function is one `omp parallel` region: setup, `rand()`, the depth-σ pass, `depth *= -1` and file writes all run on every thread. It is latent only because OpenMP is off.
- **Tidal height (`HTIDE`) is parsed but has no effect** on the SA forward model (`samodel_Rrs` uses `H`, not `H + tide`).
- **`write_nc` packs to 16-bit**, and its max-finding is buggy (`compress_2d`: `grmax = FLT_MIN`, `else if`), so the first cell's value can wrap.
- **`WRITE GRID … REGION` writes nothing**, and `.asc` output is silently NetCDF.
- `REVIEW.md` (October 2026) lists every known defect with file:line references and severity. Consult it before changing a module, and update it when fixing items.

## Verification

There is no test suite. To check a change:
1. Build with both compilers and no new warnings: `make clean && make all utils`, then `make clean && make CC=gcc-16 all utils`. For a stricter check, use `make WARN="-Wall -Wextra -Wno-unused-parameter"`; the baseline has ~140 unused-variable warnings and some `-Wsometimes-uninitialized`, so don't add to them.
2. Smoke test: a command file that reads a small `.asc`, runs `PROCESS SMOOTH`, writes and re-reads `.nc`, `PRINT GRIDSTATS`, then `EXIT`. Run it from a scratch directory, because photic writes `history.bam` to the current directory.
3. For numerical helpers, write a small standalone driver in a scratch directory that links the relevant `.c` files, and run it under `-fsanitize=address,undefined`.
4. End-to-end runs need user-supplied ACOLITE-corrected imagery and command files, which are not in the repo; ask the user for them.
