# Spectral data for the semi-analytical model

The semi-analytical model (`MODEL Lee_Semi_Analytical`, `samodel.c`) reads its spectral data from the CSV files here. They are read by `spectra.c` and are not compiled into photic.

| Data | Default | Replace with |
|---|---|---|
| Water tables a_w, b_bw, a_0, a_1 | `aw_pope_fry_1997.csv`, `bbw_morel_1974.csv`, `aphi_lee_1994.csv` | `READ WATER <file>`; each file replaces only the columns it holds |
| Bottom reflectance library | a built-in bright sand (HydroLight ooid sand, brightness 0.15–0.75) | `READ BOTTOM <file.csv>`, or `READ BOTTOM <spectra.nc> [CLASSES c1 c2 ...] [ALL]` |
| Sensor spectral response | none: tables are taken at the band wavelengths | `SCENE ... SRF <file>`: tables are averaged over each band's response |

The default water files are looked up in `$PHOTIC_DATA` if it is set, otherwise in this directory (its path is compiled into `spectra.o`). `PRINT SPECTRA` lists the tables in use and where they came from.

## File format

- Lines starting with `#` are comments; put the data's source there.
- The first other line is the header: `wavelength`, then one name per column.
- Each data row starts with a wavelength in nm, increasing down the file, then one value per column.
- Values between rows are interpolated linearly. Outside the file's range they are extrapolated linearly, so cover the bands you use.

Water files have columns named `aw`, `bbw`, `a0` and/or `a1` (units 1/m; `a0` and `a1` are dimensionless, with a_φ = (a0 + a1 ln P) P and P = a_φ(440)). Each table may have its own wavelength grid.

Bottom library files have one column per bottom type (up to 12), in the order photic numbers them (`modelled_bottom_type` is 1 for the first column). The first three are reported as the `sand`, `seagrass` and `coral` outputs. Reflectances may be absolute or relative; photic normalises each type to 1 at 550 nm. Two optional labelled rows give each type's brightness bounds, as reflectance at 550 nm:

```
brightness_min,0.156,0.042,...
brightness_max,0.796,0.159,...
```

Without them a type's bounds are 0–1. `NBOTTOMS` (default 3) sets how many types, from the first column, the model uses; it is capped at the library's size.

A NetCDF file of individual spectra, like `bottom_spectra_insitu.nc`, can also be read: `READ BOTTOM bottom_spectra_insitu.nc [CLASSES sand seagrass ...] [ALL]`. It gives one type per class:
- the shape is the median of the class's spectra, each first smoothed with an 11 nm boxcar onto 410–700 nm every 5 nm and normalised at 550 nm;
- the bounds are the 5th–95th percentiles of their reflectance at 550 nm over the whole class.

By default the types come from each class's library half (`half` = 0), in the order of the file's `classes` attribute. That matches `bottom_measured.csv` to 4×10⁻⁴. `CLASSES` picks and orders the types; `ALL` uses both halves, but that defeats tier 3's held-out test. The file needs the variables `wavelength` (nm), `reflectance(spectrum, wavelength)`, `class(spectrum)` (strings) and `half(spectrum)`.

Spectral response files (`srf/`) have one column per band, holding the relative response. A scene band uses the column whose response-weighted centroid is nearest its `WAVELENGTHS` value (within 20 nm). Its wavelength then becomes that centroid, which the CDOM and particle backscattering terms use.

## Files

| File | Contents | Source |
|---|---|---|
| `aw_pope_fry_1997.csv` | a_w, 400–800 nm every 10 nm (default) | Pope & Fry (1997), sampled every 10 nm; 710–800 nm unsourced |
| `aw_pope_fry_1997_full.csv` | a_w, 380–727.5 nm every 2.5 nm, then 730–800 nm | Pope & Fry (1997), full published table; 730–800 nm as above |
| `bbw_morel_1974.csv` | b_bw of seawater (default) | Morel (1974), via Lee's HOPE spreadsheet |
| `bbw_zhang_2009.csv` | b_bw of seawater at 20 °C, salinity 35, 400–800 nm every 5 nm | Zhang, Hu & He (2009); 12–15% below Morel |
| `aphi_lee_1994.csv` | a0, a1 (default) | Lee (1994), via Lee's HOPE spreadsheet |
| `aphi_bunbury_csiro.csv` | a0, a1, 400–750 nm every 5 nm | CSIRO measurements near Bunbury, WA (photic paper, Appendix A) |
| `bottom_hydrolight.csv` | 8 types normalised at 550 nm: HydroLight averages | HydroLight bottom files (Mobley); bounds for sand, seagrass, coral from in situ spectra |
| `bottom_measured.csv` | 7 types from 350 in situ spectra: sand, seagrass, coral, dark sediment, brown algae, green algae, turf/rubble | Five PANGAEA reef libraries (CC BY 3.0); built by `validation/bottom_library.py` |
| `bottom_spectra_insitu.nc` | All 520 spectra of the five libraries on 1 nm bins (400–800 nm), with target names, sources, photic's class and its library/held-out half. This is the data behind `bottom_measured.csv` and tier 3's `_insitu` truth; `READ BOTTOM` reads it | Five PANGAEA reef libraries (CC BY 3.0; DOIs in the file); `validation/bottom_library.py pack` |
| `bottom_other.csv` | 6 spectra formerly in `samodel.c` | Lee et al. (2001), Sagawa et al. (2012), Pinnel (Bolivar) |
| `srf/S2A_MSI.csv`, `S2B_MSI.csv`, `S2C_MSI.csv` | Sentinel-2 MSI B1–B5 | ESA S2 SRF v4.0 (2024), via ACOLITE; `validation/srf_files.py` |
| `srf/L8_OLI.csv`, `L9_OLI.csv` | Landsat 8 OLI and Landsat 9 OLI-2 B1–B4 | NASA/USGS RSR, via ACOLITE; `validation/srf_files.py` |

Each file's header comment gives the full reference.
