[![DOI](https://zenodo.org/badge/144193557.svg)](https://zenodo.org/badge/latestdoi/144193557)
[![Build Status](https://travis-ci.com/scivision/maptran.svg?branch=master)](https://travis-ci.com/scivision/maptran)
[![Build status](https://ci.appveyor.com/api/projects/status/rtmoplrumvonsscs?svg=true)](https://ci.appveyor.com/project/scivision/maptran)

# Maptran
Modern Fortran 3D coordinate conversions for geospace ecef enu eci.
Similar to Python 
[PyMap3D](https://github.com/scivision/pymap3d).

## Install

Requires Fortran 2008 compiler, such as `gfortran`, `ifort`, PGI, `nagfor`, `flang`, Cray, IBM XL, etc.

```sh
cd bin

cmake ..

cmake --build .
```

Optionally, verify Fortran functionality:
```sh
ctest -V
```

This creates `libmaptran.so` or similar, a shared library with compile-time polymorphism enabled by configuring Fortran preprocessor with one of:

* `-Drealbits=32`
* `-Drealbits=64`
* `-Drealbits=128`

## Usage

The modern Fortran API is simple like PyMap3D and Matlab Mapping Toolbox.
`elemental` procedures are used throughout to enable seamless support of scalar or array coordinate inputs. 

```fortran
use maptran

call geodetic2ecef(lat,lon,alt, x,y,z)
call geodetic2aer(lat,lon,alt, observer_lat, observer_lon, observer_alt)
```

### Functions

Popular mapping toolbox functions ported to Fortran include the
following, where the source coordinate system (before the "2") is
converted to the desired coordinate system:
```
aer2ecef  aer2enu  aer2geodetic  aer2ned
ecef2aer  ecef2enu  ecef2enuv  ecef2geodetic  ecef2ned  ecef2nedv
enu2aer  enu2ecef   enu2geodetic
geodetic2aer  geodetic2ecef  geodetic2enu  geodetic2ned
ned2aer  ned2ecef   ned2geodetic
azel2radec radec2azel
lookAtSpheroid
```

Abbreviations:

-   [AER: Azimuth, Elevation, Range](https://en.wikipedia.org/wiki/Spherical_coordinate_system)
-   [ECEF: Earth-centered, Earth-fixed](https://en.wikipedia.org/wiki/ECEF)
-   [ENU: East North Up](https://en.wikipedia.org/wiki/Axes_conventions#Ground_reference_frames:_ENU_and_NED)
-   [NED: North East Down](https://en.wikipedia.org/wiki/North_east_down)

### Caveats

* Atmospheric effects neglected in all functions.
* Planetary perturbations and nutation etc. not fully considered.

