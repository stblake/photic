# _photic_ - A Physics-Based, Satellite-Derived Bathymetry Model


## About

_photic_ is a C implementation of the physics-based, satellite-derived bathymetry model as described in (Blake, 2020). 

## Model

_photic_ is a semi-analytical model for bathymetry, water turbidity and bottom composition; which is primarily based on the physics-based model, HOPE, of (Lee, 1998)(Lee, 1999). Unlike the model of Lee, which was originally designed to use hyperspectral imagery, our model is specifically designed to use multispectral satellite imagery. In particular, we adapt to the greatly decreased spectral resolution by introducing temporal and spatial assumptions on the depth and water turbidity.

## Features 


## Open Data from _photic_

DISCLAIMER! This should be obvious, but these data from _photic_ should not be used for navigation or for any other purpose relating to safety at sea or where commercial losses could occur. 

Each of these studies used multiple LANDSAT 8 images.  

### Western Australia - Pilbara Coast from Exmouth to Whim Creek 

The modelled satellite-derived bathymetry comprised $17060$ $\text{km}^2$ at a horizontal resolution of 30m. The modelled gridded dataset has dimensions 12662 by 16077 and is projected in the WGS84 (EPSG:4326) coordinate system. The dataset is a 60 MB compressed geotiff, which can be downloaded [here](https://drive.google.com/file/d/1bmgpVht8NhCzgVcrFxv0umlSzyGk7NlD/view?usp=sharing). 

<p align="center">
  <img src="pilbara_coast_viridis.jpg" width="900" title="Pilbara Coast - Exmouth to Whim Creek Satellite-Derived Bathymetry">
</p>

### Western Australia - Exmouth Gulf

The modelled satellite-derived bathymetry comprised $4175$ $\text{km}^2$ at a horizontal resolution of 30m. The modelled gridded dataset has dimensions 3930 by 2858 and is projected in the WGS84 (EPSG:4326) coordinate system. The dataset is a 14 MB compressed geotiff, which can be downloaded [here](https://drive.google.com/file/d/10xRTwMh0dvOyOHNgeuCu6dTsOyj10Qw5/view?usp=sharing). 

<p align="center">
  <img src="exmouth_gulf_30m_viridis.jpg" width="550" title="Exmouth Gulf Satellite-Derived Bathymetry">
</p>

### Western Australia - Murion Island

The modelled satellite-derived bathymetry comprised $257$ $\text{km}^2$ at a horizontal resolution of 30m. The modelled gridded dataset has dimensions 1040 by 305 and is projected in the WGS84 (EPSG:4326) coordinate system. The dataset is a 675 KB compressed geotiff, which can be downloaded [here](https://drive.google.com/file/d/1PAQU4ZpIxWniSXUFjDP8znQfUGpBv114/view?usp=sharing).  (The image below has been rotated 90 degrees.)

<p align="center">
  <img src="murion_island_viridis.jpg" width="900" title="Murion Island Satellite-Derived Bathymetry">
</p>

### UAE - Abu Dhabi to Dubai

The modelled satellite-derived bathymetry comprised $3138$ $\text{km}^2$ at a horizontal resolution of 30m. The modelled gridded dataset has dimensions 3370 by 4700 and is projected in the WGS84 (EPSG:4326) coordinate system. The dataset is a 11.5 MB compressed geotiff.

<p align="center">
  <img src="dubai_to_abu_dhabi_viridis.jpg" width="900" title="UAE - Abu Dhabi to Dubai Satellite-Derived Bathymetry">
</p>

### UAE - Mubarraz Island to Abu Dhabi 

The modelled satellite-derived bathymetry comprised $11408$ $\text{km}^2$ at a horizontal resolution of 30m. The modelled gridded dataset has dimensions 4535 by 7476 and is projected in the WGS84 (EPSG:4326) coordinate system. The dataset is a 38 MB compressed geotiff.

<p align="center">
  <img src="abu_dhabi_to_mubarraz_island_viridis.jpg" width="900" title="UAE - Abu Dhabi to Mubarraz Island Satellite-Derived Bathymetry">
</p>

### Qatar - Eastern Qatar to Sir Baniyas Island

The modelled satellite-derived bathymetry comprised $14716$ $\text{km}^2$ at a horizontal resolution of 30m. The modelled gridded dataset has dimensions 5774 by 6905 and is projected in the WGS84 (EPSG:4326) coordinate system. The dataset is a 47 MB compressed geotiff.

<p align="center">
  <img src="H_qatar_east_viridis.jpg" width="900" title="Eastern Qatar to Sir Baniyas Island Satellite-Derived Bathymetry">
</p>

### Qatar, Bahrain and Saudi Arabia

The modelled satellite-derived bathymetry comprised $12232$ $\text{km}^2$ at a horizontal resolution of 30m. The modelled gridded dataset has dimensions 7361 by 7344 and is projected in the WGS84 (EPSG:4326) coordinate system. The dataset is a 43 MB compressed geotiff.

<p align="center">
  <img src="qatar_west_30m_viridis.jpg" width="800" title="Qatar, Bahrain, Saudi Arabia Satellite-Derived Bathymetry">
</p>

## Installation & Dependencies 

### NetCDF
NetCDF developer libraries can be installed via brew with `brew install netcdf`. 

### PGPLOT
PGPLOT can be installed (on a mac) via macports with `port install pgplot`, or compiled from source, which can be downloaded from

https://sites.astro.caltech.edu/~tjp/pgplot/

### Giza
Alternatively, Giza can be used as a drop-in replacement for the ageing PGPLOT library:

https://github.com/danieljprice/giza

### Linenoise
_photic_ uses the _linenoise_ library for the repl:

https://github.com/antirez/linenoise

### XQuartz (OSX)
For interactive graphics (on a mac) we need to install XQuartz:

https://www.xquartz.org/

### GDAL
Converting from LANDSAT and Sentinel-2 imagery into a _photic_-readable gridded file requires gdal (https://gdal.org/), which can be installed via brew with `brew install gdal`. 

### ACOLITE
All experiments with _photic_ have used the excellent atmospheric correction software, ACOLITE (https://odnature.naturalsciences.be/remsem/software-and-data/acolite): 

https://github.com/acolite/acolite


## References

(Blake, 2020) https://arxiv.org/abs/2002.02298

(Lee, 1998) Z. Lee et al, ``Hyperspectral remote sensing for shallow waters. I. A semianalytical model'', \textit{Applied Optics}, vol. 37, no. 27, pp. 6329-6338, 1998

(Lee, 1999) Z. Lee et al, ``Hyperspectral remote sensing for shallow waters: 2. Deriving bottom depths and water properties by optimization'', \textit{Applied Optics}, vol. 38, no. 18, pp. 3831-3843, 1999

