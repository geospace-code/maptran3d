[![Build Status](https://travis-ci.com/scivision/maptran.svg?branch=master)](https://travis-ci.com/scivision/maptran)
[![Build status](https://ci.appveyor.com/api/projects/status/rtmoplrumvonsscs?svg=true)](https://ci.appveyor.com/project/scivision/maptran)

# Maptran
Modern Fortran 3D coordinate conversions for geospace ecef enu eci.
Similar to Python [PyMap3D](https://github.com/scivision/pymap3d).

## Install

Requires Fortran 2003+ compiler, such as `gfortran`, `ifort`, PGI, `nagfor`, `flang`, etc.

```sh
cd bin

cmake ..

cmake --build .
```

Optionally, verify Fortran functionality:
```sh
ctest -V
```

## Usage

The modern Fortran API is simple like PyMap3D and Matlab Mapping Toolbox.
`elemental` procedures are used throughout to enable seamless support of scalar or array coordinate inputs. 
Default precision is `real64`, set at the top of `maptran.f90`.

```fortran
use maptran

call geodetic2ecef(lat,lon,alt, x,y,z)
call geodetic2aer(lat,lon,alt, observer_lat, observer_lon, observer_alt)
```
