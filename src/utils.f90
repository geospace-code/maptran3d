submodule (maptran) utils
implicit none

contains

module procedure enu2uvw
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

real(wp) :: t, lat0, lon0
logical :: d

d=.true.
if (present(deg)) d = deg

lat0 = llat0
lon0 = llon0
if (d) then
  lat0 = radians(lat0)
  lon0 = radians(lon0)
endif


t = cos(lat0) * up - sin(lat0) * north
w = sin(lat0) * up + cos(lat0) * north

u = cos(lon0) * t - sin(lon0) * east
v = sin(lon0) * t + cos(lon0) * east

end procedure enu2uvw


module procedure anglesep
! angular separation between two points on sphere
! all input/output in DEGREES

real(wp) :: la0,lo0,la1,lo1

la0 = radians(lat0)
lo0 = radians(lon0)
la1 = radians(lat1)
lo1 = radians(lon1)

anglesep = degrees(2 * asin(sqrt(haversine(la1 - la0) + &
                   cos(la0) * cos(la1) * haversine(lo1 - lo0))))

end procedure anglesep

end submodule utils
