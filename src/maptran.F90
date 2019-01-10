module maptran
use, intrinsic:: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
use, intrinsic:: iso_fortran_env, only: real32, real64, real128
implicit none
private

#if REALBITS==32
integer,parameter :: wp=real32
#elif REALBITS==64
integer,parameter :: wp=real64
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
          enu2ecef, ecef2enu, aer2geodetic, geodetic2enu,&
          geodetic2aer,enu2geodetic,degrees,radians, anglesep, &
          lookAtSpheroid

contains


pure elemental subroutine lookAtSpheroid(lat0, lon0, h0, az, tilt, lat, lon, srange, spheroid, deg)
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


elemental subroutine ecef2geodetic(x, y, z, lat, lon, alt, spheroid, deg)
!! convert ECEF (meters) to geodetic coordintes
!!
!! based on:
!! You, Rey-Jer. (2000). Transformation of Cartesian to Geodetic Coordinates without Iterations.
!! Journal of Surveying Engineering. doi: 10.1061/(ASCE)0733-9453

real(wp), intent(in) :: x,y,z
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: lat, lon
real(wp), intent(out), optional :: alt

real(wp) :: ea, eb, r, E, u, Q, huE, Beta, eps, sinBeta, cosBeta
type(Ellipsoid) :: ell
logical :: d, inside

ell = wgs84Ellipsoid
if (present(spheroid)) ell = spheroid

ea = ell%SemimajorAxis
eb = ell%SemiminorAxis

r = sqrt(x**2 + y**2 + z**2)

E = sqrt(ea**2 - eb**2)

! eqn. 4a
u = sqrt(0.5 * (r**2 - E**2) + 0.5 * sqrt((r**2 - E**2)**2 + 4 * E**2 * z**2))

Q = hypot(x, y)

huE = hypot(u, E)

!> eqn. 4b
Beta = atan2(huE / u * z, hypot(x, y))

!> final output
if (abs(beta-pi/2) <= epsilon(beta)) then !< singularity
  lat = pi/2
  cosBeta = 0._wp
  sinBeta = 1._wp
elseif (abs(beta+pi/2) <= epsilon(beta)) then !< singularity
  lat = -pi/2
  cosBeta = 0._wp
  sinBeta = -1._wp
else
  !> eqn. 13
  eps = ((eb * u - ea * huE + E**2) * sin(Beta)) / (ea * huE * 1 / cos(Beta) - E**2 * cos(Beta))
  Beta = Beta + eps

  lat = atan(ea / eb * tan(Beta))
  cosBeta = cos(Beta)
  sinBeta = sin(Beta)
endif

lon = atan2(y, x)

! eqn. 7
if (present(alt)) then
  alt = hypot(z - eb * sinBeta, Q - ea * cosBeta)
  
  !> inside ellipsoid?
  inside = x**2 / ea**2 + y**2 / ea**2 + z**2 / eb**2 < 1._wp
  if (inside) alt = -alt
endif


d=.true.
if (present(deg)) d = deg

if (d) then
  lat = degrees(lat)
  lon = degrees(lon)
endif

end subroutine ecef2geodetic


elemental subroutine geodetic2ecef(lat,lon,alt, x,y,z, spheroid, deg)
!! # geodetic2ecef   
!! convert from geodetic to ECEF coordiantes
!!
!! ## Inputs
!!
!! * lat,lon, alt:  ellipsoid geodetic coordinates of point(s) (degrees, degrees, meters)
!! * spheroid: Ellipsoid parameter struct
!! * deg: .true. degrees
!!
!! ## outputs
!!
!! * x,y,z:  ECEF coordinates of test point(s) (meters)

real(wp), intent(in) :: lat,lon !< not value due to ifort segfault bug
real(wp), intent(in) :: alt
real(wp), intent(out) :: x,y,z
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg

real(wp) :: N, sinLat, cosLat, cosLon, sinLon, lt, ln
type(Ellipsoid) :: ell
logical :: d

d = .true.
if (present(deg)) d = deg

ell = wgs84Ellipsoid
if (present(spheroid)) ell = spheroid

lt = lat
ln = lon
if (d) then
  lt = radians(lat)
  ln = radians(lon)
endif

!> Radius of curvature of the prime vertical section
N = radius_normal(lt, ell)

!! Compute cartesian (geocentric) coordinates given  (curvilinear) geodetic coordinates.

!> singularities.  Benchmark shows nearly zero runtime impact of these if statements for any real precision
if (abs(lt) <= epsilon(lt)) then
  cosLat = 1._wp
  sinLat = 0._wp
elseif (abs(lt-pi/2) <= epsilon(lt)) then
  cosLat = 0._wp
  sinLat = 1._wp
elseif (abs(lt+pi/2) <= epsilon(lt)) then
  cosLat = 0._wp
  sinLat = -1._wp
else
  cosLat = cos(lt)
  sinLat = sin(lt)
endif

if (abs(ln) <= epsilon(ln)) then
  cosLon = 1._wp
  sinLon = 0._wp
elseif (abs(ln-pi/2) <= epsilon(ln)) then
  cosLon = 0._wp
  sinLon = 1._wp
elseif (abs(ln+pi/2) <= epsilon(ln)) then
  cosLon = 0._wp
  sinLon = -1._wp
elseif (abs(ln+pi) <= epsilon(ln) .or. abs(ln-pi) <= epsilon(ln)) then
  cosLon = -1._wp
  sinLon = 0._wp
else
  cosLon = cos(ln)
  sinLon = sin(ln)
endif


x = (N + alt) * cosLat * cosLon
y = (N + alt) * cosLat * sinLon
z = (N * (ell%SemiminorAxis / ell%SemimajorAxis)**2 + alt) * sinLat

end subroutine geodetic2ecef


elemental subroutine aer2geodetic(az, el, slantRange, lat0, lon0, alt0, lat1, lon1, alt1, spheroid, deg)
!! aer2geodetic  convert azimuth, elevation, range of target from observer to geodetic coordiantes
!!
!! ## Inputs
!!
!! *  az, el, slantrange: look angles and distance to point under test (degrees, degrees, meters)
!! *  az: azimuth clockwise from local north
!! *  el: elevation angle above local horizon
!! *  lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
!! *  spheroid: Ellipsoid parameter struct
!! *  deg: .true.: degrees
!!
!! ## Outputs
!!
!! *  lat1,lon1,alt1: geodetic coordinates of test points (degrees,degrees,meters)
  
real(wp), intent(in) :: az, el, slantRange, lat0, lon0, alt0 
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: lat1,lon1,alt1

real(wp) :: x,y,z

call aer2ecef(az, el, slantRange, lat0, lon0, alt0, x, y, z, spheroid, deg)

call ecef2geodetic(x, y, z, lat1, lon1, alt1, spheroid, deg)
end subroutine aer2geodetic


elemental subroutine geodetic2aer(lat, lon, alt, lat0, lon0, alt0, az, el, slantRange, spheroid, deg)
! geodetic2aer   from an observer's perspective, convert target coordinates to azimuth, elevation, slant range.
!
! Inputs
! ------
! lat,lon, alt:  ellipsoid geodetic coordinates of point under test (degrees, degrees, meters)
! lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
! spheroid: Ellipsoid parameter struct
! angleUnit: string for angular units. Default 'd': degrees, otherwise Radians
!
! Outputs
! -------
! az, el, slantrange: look angles and distance to point under test (degrees, degrees, meters)
! az: azimuth clockwise from local north
! el: elevation angle above local horizon

real(wp), intent(in) :: lat,lon,alt, lat0, lon0, alt0 
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: az, el, slantRange

real(wp) :: east,north,up


call geodetic2enu(lat, lon, alt, lat0, lon0, alt0, east,north,up, spheroid, deg)
call enu2aer(east, north, up, az, el, slantRange, deg)
  
end subroutine geodetic2aer



elemental subroutine geodetic2enu(lat, lon, alt, lat0, lon0, alt0, east, north, up, spheroid, deg)
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

real(wp), intent(in) :: lat, lon, alt, lat0, lon0, alt0 
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: east, north, up

real(wp) x1,y1,z1,x2,y2,z2, dx,dy,dz


call geodetic2ecef(lat,lon,alt,x1,y1,z1,spheroid,deg)
call geodetic2ecef(lat0,lon0,alt0,x2,y2,z2,spheroid,deg)

dx = x1-x2;
dy = y1-y2;
dz = z1-z2;

call ecef2enuv(dx, dy, dz, lat0, lon0, east, north, up, deg)

end subroutine geodetic2enu


elemental subroutine enu2geodetic(east, north, up, lat0, lon0, alt0, lat, lon, alt, spheroid, deg)
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

real(wp), intent(in) :: east, north, up, lat0, lon0, alt0 
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: lat, lon, alt

real(wp) :: x,y,z

call enu2ecef(east, north, up, lat0, lon0, alt0, x, y, z, spheroid, deg)
call ecef2geodetic(x, y, z, lat, lon, alt, spheroid,deg)

end subroutine enu2geodetic


elemental subroutine aer2ecef(az, el, slantRange, lat0, lon0, alt0, x,y,z, spheroid, deg)
! aer2ecef  convert azimuth, elevation, range to target from observer to ECEF coordinates
!
! Inputs
! ------
! az, el, slantrange: look angles and distance to point under test (degrees, degrees, meters)
! az: azimuth clockwise from local north
! el: elevation angle above local horizon
! lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
! spheroid: Ellipsoid parameter struct
! deg: .true. degrees
!
! outputs
! -------
! x,y,z: Earth Centered Earth Fixed (ECEF) coordinates of test point (meters)

real(wp), intent(in) :: az,el, slantRange, lat0, lon0, alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: x,y,z

real(wp) :: x0,y0,z0, e,n,u,dx,dy,dz


!> Origin of the local system in geocentric coordinates.
call geodetic2ecef(lat0, lon0, alt0,x0, y0, z0, spheroid,deg)
!> Convert Local Spherical AER to ENU
call aer2enu(az, el, slantRange, e, n, u,deg)
!> Rotating ENU to ECEF
call enu2uvw(e, n, u, lat0, lon0, dx, dy, dz,deg)
!> Origin + offset from origin equals position in ECEF
x = x0 + dx
y = y0 + dy
z = z0 + dz

end subroutine aer2ecef


elemental subroutine ecef2aer(x, y, z, lat0, lon0, alt0, az, el, slantRange, spheroid, deg)
! ecef2aer  convert ECEF of target to azimuth, elevation, slant range from observer
!
! Inputs
! ------
! x,y,z: Earth Centered Earth Fixed (ECEF) coordinates of test point (meters)
! lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
! spheroid: Ellipsoid parameter struct
! deg: .true.: degrees
!
! Outputs
! -------
! az, el, slantrange: look angles and distance to point under test (degrees, degrees, meters)
! az: azimuth clockwise from local north
! el: elevation angle above local horizon

real(wp), intent(in) :: x,y,z, lat0, lon0, alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: az,el, slantRange

real(wp) :: east, north, up

call ecef2enu(x, y, z, lat0, lon0, alt0, east, north, up, spheroid, deg)
call enu2aer(east, north, up, az, el, slantRange, deg)

end subroutine ecef2aer


elemental subroutine aer2enu(az, el, slantRange, east, north, up, deg)
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

real(wp), intent(in) :: az,el  !< not value due to ifort segfault bug
real(wp), intent(in) :: slantRange
logical, intent(in), optional :: deg
real(wp),intent(out) :: east, north,up

real(wp) :: r, a, e
logical :: d

d=.true.
if (present(deg)) d = deg

a = az
e = el
if (d) then
  a = radians(az)
  e = radians(el)
endif

!> Calculation of AER2ENU
up = slantRange * sin(e)
r = slantRange * cos(e)
east = r * sin(a)
north = r * cos(a)

end subroutine aer2enu


elemental subroutine enu2aer(east, north, up, az, elev, slantRange, deg)
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

real(wp),intent(in) :: east, north, up
logical, intent(in), optional :: deg
real(wp), intent(out) :: az, elev, slantRange

real(wp) :: r, e, n, u
logical :: d

real(wp), parameter :: tolerance = 1e-3_wp !< 1mm precision

!> singularity fixes
e = east
if (abs(e) < tolerance) e = 0._wp

n = north
if (abs(n) < tolerance) n = 0._wp

u = up
if (abs(u) < tolerance) u = 0._wp


r = hypot(e, n)
slantRange = hypot(r, u)

!> radians
elev = atan2(u, r)
az = modulo(atan2(e, n), 2._wp * pi)

d=.true.
if (present(deg)) d = deg

if (d) then
  elev = degrees(elev)
  az = degrees(az)
endif
  
end subroutine enu2aer


elemental subroutine enu2ecef(e, n, u, lat0, lon0, alt0, x, y, z, spheroid, deg)
! enu2ecef  convert from ENU to ECEF coordiantes
!
! Inputs
! ------
! e,n,u:  East, North, Up coordinates of test points (meters)
! lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
! spheroid: Ellipsoid parameter struct
! angleUnit: string for angular units. Default 'd': degrees
!
! outputs
! -------
! x,y,z: Earth Centered Earth Fixed (ECEF) coordinates of test point (meters)
real(wp), intent(in) :: e,n,u,lat0,lon0,alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: x,y,z   

real(wp) :: x0,y0,z0,dx,dy,dz           
          

call geodetic2ecef(lat0, lon0, alt0, x0, y0, z0, spheroid, deg)
call enu2uvw(e, n, u, lat0, lon0, dx, dy, dz, deg)

 x = x0 + dx
 y = y0 + dy
 z = z0 + dz
end subroutine enu2ecef



elemental subroutine ecef2enu(x, y, z, lat0, lon0, alt0, east, north, up, spheroid, deg)
!! ecef2enu  convert ECEF to ENU
!!
!! ## Inputs
!!
!! * x,y,z: Earth Centered Earth Fixed (ECEF) coordinates of test point (meters)
!! * lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
!! * spheroid: Ellipsoid parameter struct
!! * angleUnit: string for angular units. Default 'd': degrees
!!
!! ## outputs
!!
!! * e,n,u:  East, North, Up coordinates of test points (meters)
  
real(wp), intent(in) :: x,y,z,lat0,lon0,alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real(wp), intent(out) :: east,north,up

real(wp) :: x0,y0,z0

call geodetic2ecef(lat0, lon0, alt0, x0,y0,z0, spheroid,deg)
call ecef2enuv(x - x0, y - y0, z - z0, lat0, lon0, east, north, up, deg)
end subroutine ecef2enu


elemental subroutine ecef2enuv(u, v, w, lat0, lon0, east, north, up, deg)
!! ecef2enuv convert *vector projection* UVW to ENU
!!
!! ## Inputs
!!
!! * u,v,w: meters
!! * lat0,lon0: geodetic latitude and longitude (degrees)
!! * deg: .true. degrees
!!
!! ## Outputs
!!
!! * east,north,Up:  East, North, Up vector

real(wp), intent(in) :: u,v,w
real(wp), value :: lat0,lon0
logical, intent(in), optional :: deg
real(wp), intent(out) :: east, north, up

real(wp) :: t
logical :: d

d=.true.
if (present(deg)) d = deg

if (d) then
  lat0 = radians(lat0)
  lon0 = radians(lon0)
endif

t     =  cos(lon0) * u + sin(lon0) * v
east  = -sin(lon0) * u + cos(lon0) * v
up    =  cos(lat0) * t + sin(lat0) * w
north = -sin(lat0) * t + cos(lat0) * w
end subroutine ecef2enuv


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


elemental real(wp) function radius_normal(lat,E)

real(wp), intent(in) :: lat
type(Ellipsoid), intent(in) :: E

!> singularity  pi/2 issue is inherent to real32
if (abs(lat) <= epsilon(lat)) then
  radius_normal = E%SemimajorAxis
else
  radius_normal = E%SemimajorAxis**2 / sqrt( E%SemimajorAxis**2 * cos(lat)**2 + E%SemiminorAxis**2 * sin(lat)**2 )
endif

end function radius_normal


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
