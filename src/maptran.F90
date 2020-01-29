module maptran
use, intrinsic:: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
use, intrinsic:: iso_fortran_env, only: real32, real64
#if REALBITS==128
use, intrinsic :: iso_fortran_env, only : real128
#endif

implicit none
private

#if REALBITS==32
integer,parameter :: wp=real32
#elif REALBITS==128
integer,parameter :: wp=real128
#else
integer, parameter :: wp=real64
#endif

type,public :: Ellipsoid
  real(wp) :: SemimajorAxis, Flattening, SemiminorAxis
end type

real(wp), parameter :: pi = 4._wp * atan(1.0_wp)

type(Ellipsoid), parameter, public :: wgs84Ellipsoid = &
      Ellipsoid(SemimajorAxis=6378137._wp, &
                SemiminorAxis=6378137._wp * (1._wp - 1._wp / 298.2572235630_wp), &
                Flattening = 1. / 298.2572235630_wp)

public :: wp,ecef2geodetic, geodetic2ecef, aer2enu, enu2aer, aer2ecef, ecef2aer, &
          enu2ecef, ecef2enu, aer2geodetic, geodetic2enu, enu2uvw,&
          geodetic2aer,enu2geodetic,degrees,radians, anglesep, &
          lookAtSpheroid

interface ! aer
module elemental subroutine ecef2aer(x, y, z, lat0, lon0, alt0, az, el, slantRange, spheroid, deg)
real(wp), intent(in) :: x,y,z, lat0, lon0, alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: az,el, slantRange
end subroutine ecef2aer

module elemental subroutine aer2ecef(az, el, slantRange, lat0, lon0, alt0, x,y,z, spheroid, deg)
real(wp), intent(in) :: az,el, slantRange, lat0, lon0, alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: x,y,z
end subroutine aer2ecef

module elemental subroutine geodetic2aer(lat, lon, alt, lat0, lon0, alt0, az, el, slantRange, spheroid, deg)
real(wp), intent(in) :: lat,lon,alt, lat0, lon0, alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: az, el, slantRange
end subroutine geodetic2aer

module elemental subroutine aer2geodetic(az, el, slantRange, lat0, lon0, alt0, lat1, lon1, alt1, spheroid, deg)
real(wp), intent(in) :: az, el, slantRange, lat0, lon0, alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: lat1,lon1,alt1
end subroutine aer2geodetic

end interface

interface ! ecef

module elemental subroutine ecef2geodetic(x, y, z, lat, lon, alt, spheroid, deg)
real(wp), intent(in) :: x,y,z
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: lat, lon
real(wp), intent(out), optional :: alt
end subroutine ecef2geodetic

module elemental subroutine geodetic2ecef(lat,lon,alt, x,y,z, spheroid, deg)
real(wp), intent(in) :: lat,lon !< not value due to ifort segfault bug
real(wp), intent(in) :: alt
real(wp), intent(out) :: x,y,z
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
end subroutine geodetic2ecef

module elemental subroutine enu2ecef(e, n, u, lat0, lon0, alt0, x, y, z, spheroid, deg)
real(wp), intent(in) :: e,n,u,lat0,lon0,alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: x,y,z
end subroutine enu2ecef

module elemental subroutine ecef2enu(x, y, z, lat0, lon0, alt0, east, north, up, spheroid, deg)
real(wp), intent(in) :: x,y,z,lat0,lon0,alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: east,north,up
end subroutine ecef2enu

module elemental subroutine ecef2enuv(u, v, w, lat0, lon0, east, north, up, deg)
real(wp), intent(in) :: u,v,w
real(wp), value :: lat0,lon0
logical, intent(in), optional :: deg
real(wp), intent(out) :: east, north, up
end subroutine ecef2enuv

end interface

interface ! enu

module elemental subroutine geodetic2enu(lat, lon, alt, lat0, lon0, alt0, east, north, up, spheroid, deg)
real(wp), intent(in) :: lat, lon, alt, lat0, lon0, alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: east, north, up
end subroutine geodetic2enu

module elemental subroutine enu2geodetic(east, north, up, lat0, lon0, alt0, lat, lon, alt, spheroid, deg)
real(wp), intent(in) :: east, north, up, lat0, lon0, alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: lat, lon, alt
end subroutine enu2geodetic

module elemental subroutine aer2enu(az, el, slantRange, east, north, up, deg)
real(wp), intent(in) :: az,el  !< not value due to ifort segfault bug
real(wp), intent(in) :: slantRange
logical, intent(in), optional :: deg
real(wp),intent(out) :: east, north,up
end subroutine aer2enu

module elemental subroutine enu2aer(east, north, up, az, elev, slantRange, deg)
real(wp),intent(in) :: east, north, up
logical, intent(in), optional :: deg
real(wp), intent(out) :: az, elev, slantRange
end subroutine enu2aer

end interface

contains


elemental subroutine lookAtSpheroid(lat0, lon0, h0, az, tilt, lat, lon, srange, spheroid, deg)
!! lookAtSpheroid: Calculates line-of-sight intersection with Earth (or other ellipsoid) surface from above surface / orbit
!!
!! ## Inputs
!! * lat0, lon0: latitude and longitude of starting point
!! * h0: altitude of starting point in meters
!! * az: azimuth angle of line-of-sight, clockwise from North
!! * tilt: tilt angle of line-of-sight with respect to local vertical (nadir = 0)
!!
!! ## Outputs
!! * lat, lon: latitude and longitude where the line-of-sight intersects with the Earth ellipsoid
!! * d: slant range in meters from the starting point to the intersect point
!!
!!  Values will be NaN if the line of sight does not intersect.
!!  Algorithm based on:
!!  https://medium.com/@stephenhartzell/satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6 Stephen Hartzell

real(wp), intent(in) :: lat0, lon0, h0, az, tilt
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: lat, lon, srange

type(Ellipsoid) :: ell
real(wp) :: a,b,c,east,north,up,u,v,w,x,y,z, el, val, radical, magnitude
logical :: d
real(wp):: nan

nan = ieee_value(0._wp, ieee_quiet_nan)

ell = wgs84Ellipsoid
if (present(spheroid)) ell = spheroid

d=.true.
if (present(deg)) d = deg

 a = ell%SemimajorAxis
 b = a
 c = ell%SemiminorAxis

if (d) then
  el = tilt - 90._wp
else
  el = tilt - pi / 2
endif

call aer2enu(az, el, 1._wp, east, north, up, deg=d) ! fixed 1 km slant range
call enu2uvw(east,north,up, lat0,lon0, u,v,w, deg=d)
call geodetic2ecef(lat0,lon0,h0, x,y,z, spheroid=ell, deg=d)

val = -a**2 * b**2 * w * z - a**2 * c**2 * v * y - b**2 * c**2 * u * x

radical = a**2 * b**2 * w**2 + a**2 * c**2 * v**2 - a**2 * v**2 * z**2 + 2 * a**2 * v * w * y * z - &
          a**2 * w**2 * y**2 + b**2 * c**2 * u**2 - b**2 * u**2 * z**2 + 2 * b**2 * u * w * x * z - &
          b**2 * w**2 * x**2 - c**2 * u**2 * y**2 + 2 * c**2 * u * v * x * y - c**2 * v**2 * x**2

magnitude = a**2 * b**2 * w**2 + a**2 * c**2 * v**2 + b**2 * c**2 * u**2

!> Return nan if radical < 0 or d < 0 because LOS vector does not point towards Earth
if (radical > 0) then
  srange = (val - a * b * c * sqrt(radical)) / magnitude
else
  srange = nan
endif

if (srange < 0) srange = nan

call ecef2geodetic(x+srange*u, y+srange*v, z+srange*w, lat, lon, spheroid=ell, deg=d)

end subroutine lookAtSpheroid


elemental subroutine enu2uvw(east,north,up, lat0,lon0, u,v,w, deg)
!! # enu2uvw   convert from ENU to UVW coordinates
!!
!! ## Inputs
!!
!! * e,n,up:  East, North, Up coordinates of point(s) (meters)
!! * lat0,lon0: geodetic coordinates of observer/reference point (degrees)
!! * deg: ,true. degrcees
!!
!! ## outputs
!!
!! * u,v,w:   coordinates of test point(s) (meters)

real(wp), intent(in) :: east,north,up
real(wp), value :: lat0,lon0
real(wp), intent(out) :: u,v,w
logical, intent(in), optional :: deg

real(wp) :: t
logical :: d

d=.true.
if (present(deg)) d = deg

if (d) then
  lat0 = radians(lat0)
  lon0 = radians(lon0)
endif


t = cos(lat0) * up - sin(lat0) * north
w = sin(lat0) * up + cos(lat0) * north

u = cos(lon0) * t - sin(lon0) * east
v = sin(lon0) * t + cos(lon0) * east

end subroutine enu2uvw


elemental real(wp) function anglesep(lon0,lat0,lon1,lat1)
! angular separation between two points on sphere
! all input/output in DEGREES

real(wp), intent(in) :: lat0,lon0,lat1,lon1
real(wp) :: la0,lo0,la1,lo1

la0 = radians(lat0)
lo0 = radians(lon0)
la1 = radians(lat1)
lo1 = radians(lon1)

anglesep = degrees(2 * asin(sqrt(haversine(la1 - la0) + &
                   cos(la0) * cos(la1) * haversine(lo1 - lo0))))

end function anglesep


elemental real(wp) function haversine(theta)
!! theta: angle in RADIANS
real(wp), intent(in) :: theta

haversine =  (1 - cos(theta)) / 2._wp

end function haversine


elemental real(wp) function degrees(rad)
real(wp), intent(in) :: rad

degrees = 180._wp / pi * rad
end function degrees


elemental real(wp) function radians(deg)
real(wp), intent(in) :: deg

radians = pi / 180._wp * deg
end function radians

end module
