module maptran

use, intrinsic:: ieee_arithmetic, only: ieee_quiet_nan, ieee_value

implicit none (type, external)
private
public :: pi, ecef2geodetic, geodetic2ecef, aer2enu, enu2aer, aer2ecef, ecef2aer, &
          enu2ecef, ecef2enu, ecef2enuv, aer2geodetic, geodetic2enu, enu2uvw,&
          geodetic2aer,enu2geodetic,degrees,radians, anglesep, &
          lookAtSpheroid, Ellipsoid, &
          haversine

type :: Ellipsoid
  real :: SemimajorAxis, Flattening, SemiminorAxis
end type

real, parameter :: pi = 4 * atan(1.)

type(Ellipsoid), parameter, public :: wgs84Ellipsoid = &
  Ellipsoid(SemimajorAxis=6378137., &
            SemiminorAxis=6378137. * (1 - 1. / 298.2572235630), &
            Flattening = 1. / 298.2572235630)

interface ! aer
module elemental subroutine ecef2aer(x, y, z, lat0, lon0, alt0, az, el, slantRange, spheroid, deg)
real, intent(in) :: x,y,z, lat0, lon0, alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real, intent(out) :: az,el, slantRange
end subroutine ecef2aer

module elemental subroutine aer2ecef(az, el, slantRange, lat0, lon0, alt0, x,y,z, spheroid, deg)
real, intent(in) :: az,el, slantRange, lat0, lon0, alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real, intent(out) :: x,y,z
end subroutine aer2ecef

module elemental subroutine geodetic2aer(lat, lon, alt, lat0, lon0, alt0, az, el, slantRange, spheroid, deg)
real, intent(in) :: lat,lon,alt, lat0, lon0, alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real, intent(out) :: az, el, slantRange
end subroutine geodetic2aer

module elemental subroutine aer2geodetic(az, el, slantRange, lat0, lon0, alt0, lat1, lon1, alt1, spheroid, deg)
real, intent(in) :: az, el, slantRange, lat0, lon0, alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real, intent(out) :: lat1,lon1,alt1
end subroutine aer2geodetic

end interface

interface ! ecef

module elemental subroutine ecef2geodetic(x, y, z, lat, lon, alt, spheroid, deg)
real, intent(in) :: x,y,z
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real, intent(out) :: lat, lon
real, intent(out), optional :: alt
end subroutine ecef2geodetic

module elemental subroutine geodetic2ecef(llat,llon,alt, x,y,z, spheroid, deg)
real, intent(in) :: llat,llon, alt
real, intent(out) :: x,y,z
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
end subroutine geodetic2ecef

module elemental subroutine enu2ecef(e, n, u, lat0, lon0, alt0, x, y, z, spheroid, deg)
real, intent(in) :: e,n,u,lat0,lon0,alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real, intent(out) :: x,y,z
end subroutine enu2ecef

module elemental subroutine ecef2enu(x, y, z, lat0, lon0, alt0, east, north, up, spheroid, deg)
real, intent(in) :: x,y,z,lat0,lon0,alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real, intent(out) :: east,north,up
end subroutine ecef2enu

module elemental subroutine ecef2enuv(u, v, w, llat0, llon0, east, north, up, deg)
real, intent(in) :: u,v,w
real, intent(in) :: llat0,llon0
logical, intent(in), optional :: deg
real, intent(out) :: east, north, up
end subroutine ecef2enuv

end interface

interface ! enu

module elemental subroutine geodetic2enu(lat, lon, alt, lat0, lon0, alt0, east, north, up, spheroid, deg)
real, intent(in) :: lat, lon, alt, lat0, lon0, alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real, intent(out) :: east, north, up
end subroutine geodetic2enu

module elemental subroutine enu2geodetic(east, north, up, lat0, lon0, alt0, lat, lon, alt, spheroid, deg)
real, intent(in) :: east, north, up, lat0, lon0, alt0
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real, intent(out) :: lat, lon, alt
end subroutine enu2geodetic

module elemental subroutine aer2enu(az, el, slantRange, east, north, up, deg)
real, intent(in) :: az,el  !< not value due to ifort segfault bug
real, intent(in) :: slantRange
logical, intent(in), optional :: deg
real,intent(out) :: east, north,up
end subroutine aer2enu

module elemental subroutine enu2aer(east, north, up, az, elev, slantRange, deg)
real,intent(in) :: east, north, up
logical, intent(in), optional :: deg
real, intent(out) :: az, elev, slantRange
end subroutine enu2aer


module elemental subroutine lookAtSpheroid(lat0, lon0, h0, az, tilt, lat, lon, srange, spheroid, deg)
real, intent(in) :: lat0, lon0, h0, az, tilt
type(Ellipsoid), intent(in), optional :: spheroid
logical, intent(in), optional :: deg
real, intent(out) :: lat, lon, srange
end subroutine lookAtSpheroid

end interface

interface ! utils.f90
module elemental subroutine enu2uvw(east,north,up, llat0,llon0, u,v,w, deg)
real, intent(in) :: east,north,up
real, intent(in) :: llat0,llon0
real, intent(out) :: u,v,w
logical, intent(in), optional :: deg
end subroutine enu2uvw

module elemental real function anglesep(lon0,lat0,lon1,lat1)
real, intent(in) :: lat0,lon0,lat1,lon1
end function anglesep
end interface

contains


elemental real function haversine(theta)
!! theta: angle in RADIANS
real, intent(in) :: theta

haversine =  (1 - cos(theta)) / 2

end function haversine


elemental real function degrees(rad)
real, intent(in) :: rad

degrees = 180 / pi * rad
end function degrees


elemental real function radians(deg)
real, intent(in) :: deg

radians = pi / 180 * deg
end function radians

end module maptran
