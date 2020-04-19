# Maptran 3D

[![DOI](https://zenodo.org/badge/144193557.svg)](https://zenodo.org/badge/latestdoi/144193557)
![Actions Status](https://github.com/geospace-code/maptran3d/workflows/ci/badge.svg)

Modern Fortran 3D coordinate conversions for geospace ecef enu eci.
Similar to Python
[PyMap3D](https://github.com/scivision/pymap3d)
and Matlab / GNU Octave
[Matmap3d](https://github.com/geospace-code/matmap3d).

## Install

Requires Fortran 2008 compiler, such as `gfortran`, `ifort`, PGI, `nagfor`, `flang`, Cray, IBM XL, etc.
Use CMake or Meson to build the suite, which creates `libmaptran.a` or similar.
Compile-time polymorphism enabled by configuring with one of:

* `-Drealbits=32`
* `-Drealbits=64`

The large real values typical of map coordinates can lead to large error with 32-bit reals.
64-bit real is the default.

### CMake

```sh
ctest -S setup.cmake -VV
```

### Meson

```sh
meson build

meson test -C build
```

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

* [AER: Azimuth, Elevation, Range](https://en.wikipedia.org/wiki/Spherical_coordinate_system)
* [ECEF: Earth-centered, Earth-fixed](https://en.wikipedia.org/wiki/ECEF)
* [ENU: East North Up](https://en.wikipedia.org/wiki/Axes_conventions#Ground_reference_frames:_ENU_and_NED)
* [NED: North East Down](https://en.wikipedia.org/wiki/North_east_down)

### Caveats

* Atmospheric effects neglected in all functions.
* Planetary perturbations and nutation etc. not fully considered.
