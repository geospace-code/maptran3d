submodule (maptran) enu
implicit none

contains

module procedure geodetic2enu
!! geodetic2enu    convert from geodetic to ENU coordinates
!!
!! ## Inputs
!!
!! *  lat,lon, alt:  ellipsoid geodetic coordinates of point under test (degrees, degrees, meters)
!! *  lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
!! *  spheroid: Ellipsoid parameter struct
!! *  angleUnit: string for angular units. Default 'd': degrees
!!
!! ## outputs
!!
!! *  e,n,u:  East, North, Up coordinates of test points (meters)

real(wp) x1,y1,z1,x2,y2,z2, dx,dy,dz


call geodetic2ecef(lat,lon,alt,x1,y1,z1,spheroid,deg)
call geodetic2ecef(lat0,lon0,alt0,x2,y2,z2,spheroid,deg)

dx = x1-x2
dy = y1-y2
dz = z1-z2

call ecef2enuv(dx, dy, dz, lat0, lon0, east, north, up, deg)

end procedure geodetic2enu


module procedure enu2geodetic
! enu2geodetic   convert from ENU to geodetic coordinates
!
! Inputs
! ------
!  East, North, Up: coordinates of point(s) (meters)
! lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
! spheroid: Ellipsoid parameter struct
! deg: .true. degrees
!
! outputs
! -------
! lat,lon,alt: geodetic coordinates of test points (degrees,degrees,meters)

real(wp) :: x,y,z

call enu2ecef(east, north, up, lat0, lon0, alt0, x, y, z, spheroid, deg)
call ecef2geodetic(x, y, z, lat, lon, alt, spheroid,deg)

end procedure enu2geodetic


module procedure aer2enu
! aer2enu  convert azimuth, elevation, range to ENU coordinates
!
! Inputs
! ------
! az, el, slantrange: look angles and distance to point under test (degrees, degrees, meters)
! az: azimuth clockwise from local north
! el: elevation angle above local horizon
! angleUnit: string for angular units. Default 'd': degrees
!
! Outputs
! -------
! e,n,u:  East, North, Up coordinates of test points (meters)

real(wp) :: r, a, e
logical :: d

d=.true.
if (present(deg)) d = deg

a = az
e = el
if (d) then
  a = radians(a)
  e = radians(e)
endif

!> Calculation of AER2ENU
up = slantRange * sin(e)
r = slantRange * cos(e)
east = r * sin(a)
north = r * cos(a)

end procedure aer2enu


module procedure enu2aer
!! enu2aer   convert ENU to azimuth, elevation, slant range
!!
!! ## Inputs
!!
!! *  e,n,u:  East, North, Up coordinates of test points (meters)
!! *  deg: .true. degrees
!!
!! ## outputs
!!
!! *  az, el, slantrange: look angles and distance to point under test (degrees, degrees, meters)
!! *  az: azimuth clockwise from local north
!! *  el: elevation angle above local horizon

real(wp) :: r, e, n, u
logical :: d

real(wp), parameter :: tolerance = 1e-3_wp !< 1mm precision

!> singularity fixes
e = east
if (abs(e) < tolerance) e = 0

n = north
if (abs(n) < tolerance) n = 0

u = up
if (abs(u) < tolerance) u = 0

r = hypot(e, n)
slantRange = hypot(r, u)

!> radians
elev = atan2(u, r)
az = modulo(atan2(e, n), 2 * pi)

d=.true.
if (present(deg)) d = deg

if (d) then
  elev = degrees(elev)
  az = degrees(az)
endif

end procedure enu2aer

end submodule enu
