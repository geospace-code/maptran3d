submodule (maptran) sphere
implicit none

contains

module procedure lookAtSpheroid
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

end procedure lookAtSpheroid



end submodule sphere