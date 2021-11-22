submodule (maptran) ecef

implicit none (type, external)

contains

module procedure ecef2geodetic
!! convert ECEF (meters) to geodetic coordintes
!!
!! based on:
!! You, Rey-Jer. (2000). Transformation of Cartesian to Geodetic Coordinates without Iterations.
!! Journal of Surveying Engineering. doi: 10.1061/(ASCE)0733-9453

real :: ea, eb, r, E, u, Q, huE, Beta, eps, sinBeta, cosBeta
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
  cosBeta = 0
  sinBeta = 1
elseif (abs(beta+pi/2) <= epsilon(beta)) then !< singularity
  lat = -pi/2
  cosBeta = 0
  sinBeta = -1
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
  inside = x**2 / ea**2 + y**2 / ea**2 + z**2 / eb**2 < 1
  if (inside) alt = -alt
endif


d=.true.
if (present(deg)) d = deg

if (d) then
  lat = degrees(lat)
  lon = degrees(lon)
endif

end procedure ecef2geodetic


module procedure geodetic2ecef
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

real :: N, sinLat, cosLat, cosLon, sinLon, lat, lon
type(Ellipsoid) :: ell
logical :: d

d = .true.
if (present(deg)) d = deg

ell = wgs84Ellipsoid
if (present(spheroid)) ell = spheroid

lat = llat
lon = llon
if (d) then
  lat = radians(lat)
  lon = radians(lon)
endif

!> Radius of curvature of the prime vertical section
N = radius_normal(lat, ell)

!! Compute cartesian (geocentric) coordinates given  (curvilinear) geodetic coordinates.

!> singularities.  Benchmark shows nearly zero runtime impact of these if statements for any real precision
if (abs(lat) <= epsilon(lat)) then
  cosLat = 1
  sinLat = 0
elseif (abs(lat-pi/2) <= epsilon(lat)) then
  cosLat = 0
  sinLat = 1
elseif (abs(lat+pi/2) <= epsilon(lat)) then
  cosLat = 0
  sinLat = -1
else
  cosLat = cos(lat)
  sinLat = sin(lat)
endif

if (abs(lon) <= epsilon(lon)) then
  cosLon = 1
  sinLon = 0
elseif (abs(lon-pi/2) <= epsilon(lon)) then
  cosLon = 0
  sinLon = 1
elseif (abs(lon+pi/2) <= epsilon(lon)) then
  cosLon = 0
  sinLon = -1
elseif (abs(lon+pi) <= epsilon(lon) .or. abs(lon-pi) <= epsilon(lon)) then
  cosLon = -1
  sinLon = 0
else
  cosLon = cos(lon)
  sinLon = sin(lon)
endif

x = (N + alt) * cosLat * cosLon
y = (N + alt) * cosLat * sinLon
z = (N * (ell%SemiminorAxis / ell%SemimajorAxis)**2 + alt) * sinLat

end procedure geodetic2ecef


module procedure enu2ecef
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

real :: x0,y0,z0,dx,dy,dz

call geodetic2ecef(lat0, lon0, alt0, x0, y0, z0, spheroid, deg)
call enu2uvw(e, n, u, lat0, lon0, dx, dy, dz, deg)

 x = x0 + dx
 y = y0 + dy
 z = z0 + dz
end procedure enu2ecef


module procedure ecef2enu
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

real :: x0,y0,z0

call geodetic2ecef(lat0, lon0, alt0, x0,y0,z0, spheroid,deg)
call ecef2enuv(x - x0, y - y0, z - z0, lat0, lon0, east, north, up, deg)
end procedure ecef2enu


module procedure ecef2enuv
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

real :: t, lat0, lon0
logical :: d

d=.true.
if (present(deg)) d = deg

lat0 = llat0
lon0 = llon0
if (d) then
  lat0 = radians(lat0)
  lon0 = radians(lon0)
endif

t     =  cos(lon0) * u + sin(lon0) * v
east  = -sin(lon0) * u + cos(lon0) * v
up    =  cos(lat0) * t + sin(lat0) * w
north = -sin(lat0) * t + cos(lat0) * w
end procedure ecef2enuv


elemental real function radius_normal(lat,E)

real, intent(in) :: lat
type(Ellipsoid), intent(in) :: E

!> singularity  pi/2 issue is inherent to real32
if (abs(lat) <= epsilon(lat)) then
  radius_normal = E%SemimajorAxis
else
  radius_normal = E%SemimajorAxis**2 / sqrt( E%SemimajorAxis**2 * cos(lat)**2 + E%SemiminorAxis**2 * sin(lat)**2 )
endif

end function radius_normal

end submodule ecef
