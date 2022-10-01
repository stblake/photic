# _photic_ - A Physics-Based, Satellite-Derived Bathymetry Model


## About

_photic_ is a C implementation of the physics-based, satellite-derived bathymetry model as described in (Blake, 2020). 

## Model


## Features 


## Open Data from _photic_

DISCLAIMER! This should be obvious, but these data from _photic_ should not be used for navigation or for any other purpose relating to safety at sea or where commercial losses could occur. 

Each of these studies used multiple LANDSAT 8 images.  

### Western Australia - Exmouth to Whim Creek 

The modelled satellite-derived bathymetry comprised $17060$ $\text{km}^2$ at a horizontal resolution of 30m. 

<p align="center">
  <img src="pilbara_coast_viridis.jpg" width="1000" title="Pilbara Coast - Exmouth to Whim Creek Satellite-Derived Bathymetry">
</p>

### Western Australia - Exmouth Gulf

The modelled satellite-derived bathymetry comprised $4175$ $\text{km}^2$ at a horizontal resolution of 30m. 

<p align="center">
  <img src="exmouth_gulf_30m_viridis.jpg" width="550" title="Exmouth Gulf Satellite-Derived Bathymetry">
</p>

### Western Australia - Murion Island

The modelled satellite-derived bathymetry comprised $257$ $\text{km}^2$ at a horizontal resolution of 30m. (The image below has been rotated 90 degrees.)

<p align="center">
  <img src="murion_island_viridis.jpg" width="1000" title="Murion Island Satellite-Derived Bathymetry">
</p>

### UAE - Abu Dhabi to Dubai

The modelled satellite-derived bathymetry comprised $3138$ $\text{km}^2$ at a horizontal resolution of 30m. 

<p align="center">
  <img src="dubai_to_abu_dhabi_viridis.jpg" width="1000" title="UAE - Abu Dhabi to Dubai Satellite-Derived Bathymetry">
</p>

### UAE - Mubarraz Island to Abu Dhabi 

The modelled satellite-derived bathymetry comprised $11408$ $\text{km}^2$ at a horizontal resolution of 30m. 

<p align="center">
  <img src="abu_dhabi_to_mubarraz_island_viridis.jpg" width="1000" title="UAE - Abu Dhabi to Mubarraz Island Satellite-Derived Bathymetry">
</p>

### Qatar - Eastern Qatar to Sir Baniyas Island

The modelled satellite-derived bathymetry comprised $14716$ $\text{km}^2$ at a horizontal resolution of 30m. 

<p align="center">
  <img src="H_qatar_east_viridis.jpg" width="1000" title="Eastern Qatar to Sir Baniyas Island Satellite-Derived Bathymetry">
</p>

### Qatar, Bahrain and Saudi Arabia

The modelled satellite-derived bathymetry comprised $12232$ $\text{km}^2$ at a horizontal resolution of 30m. 

<p align="center">
  <img src="qatar_west_30m_viridis.jpg" width="1000" title="Qatar, Bahrain, Saudi Arabia Satellite-Derived Bathymetry">
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
