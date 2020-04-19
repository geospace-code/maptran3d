submodule (maptran) aer

implicit none (type, external)

contains


module procedure ecef2aer
!! ecef2aer  convert ECEF of target to azimuth, elevation, slant range from observer
!!
!! Inputs
!! ------
!! x,y,z: Earth Centered Earth Fixed (ECEF) coordinates of test point (meters)
!! lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
!! spheroid: Ellipsoid parameter struct
!! deg: .true.: degrees
!!
!! Outputs
!! -------
!! az, el, slantrange: look angles and distance to point under test (degrees, degrees, meters)
!! az: azimuth clockwise from local north
!! el: elevation angle above local horizon

real(wp) :: east, north, up

call ecef2enu(x, y, z, lat0, lon0, alt0, east, north, up, spheroid, deg)
call enu2aer(east, north, up, az, el, slantRange, deg)

end procedure ecef2aer


module procedure aer2ecef
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

end procedure aer2ecef


module procedure geodetic2aer
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

real(wp) :: east,north,up


call geodetic2enu(lat, lon, alt, lat0, lon0, alt0, east,north,up, spheroid, deg)
call enu2aer(east, north, up, az, el, slantRange, deg)

end procedure geodetic2aer


module procedure aer2geodetic
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

real(wp) :: x,y,z

call aer2ecef(az, el, slantRange, lat0, lon0, alt0, x, y, z, spheroid, deg)

call ecef2geodetic(x, y, z, lat1, lon1, alt1, spheroid, deg)
end procedure aer2geodetic

end submodule aer
