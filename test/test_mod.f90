program test_maptran
!! error tolerances are set for single precision, double precision is much more precise.
use, intrinsic:: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
use, intrinsic:: iso_fortran_env, only: stderr=>error_unit
use maptran
use vallado
use assert, only: assert_isclose

implicit none (type, external)

type(Ellipsoid), parameter :: spheroid = wgs84Ellipsoid

real, parameter :: &
lat = 42, lon= -82, alt = 200, &
az = 33, el=70, rng= 1e3, &
x0 = 660.6752518e3, y0 = -4700.9486832e3, z0 = 4245.7376622e3, &
xl = 6.609301927610815e+5, yl = -4.701424222957011e6, zl = 4.246579604632881e+06, & !< aer2ecef
east = 186.277521, north = 286.84222, up = 939.69262, & !< aer2enu
ha = 45.482789587392013

type(datetime), parameter :: t0 = datetime(2014,4,6,8,0,0) !< UTC

real :: ea, eb, atol_dist, atol_deg

select case (storage_size(0.))
  case (32)
    atol_dist = 1 !< 1 meter
    atol_deg = 0.1
  case (64)
    atol_dist = 0.001 !< 1 mm
    atol_deg = 0.01
  case (128)
    atol_dist = 0.001 !< 1 mm
    atol_deg = 1e-6
  case default
    error stop "unknown precision"
end select

print *,'BEGIN: Maptran ',storage_size(0.),' bits.'
print '(A,2ES12.3)', 'Abs_tol dist, deg:',atol_dist, atol_deg

ea = spheroid%SemimajorAxis
eb = spheroid%SemiminorAxis

!! scalar degrees

call test_enu2aer(east,north, up, az,el,rng)
print *, "OK: enu2aer"

call test_aer2ecef(az, el, rng, lat, lon, alt)
print *, "OK: aer2ecef"

call test_aer2enu(az, el, rng, east, north, up)
print *, "OK: aer2enu"

call test_geodetic2ecef(lat, lon, alt)
print *, "OK: geodetic2ecef"

call test_ecef2geodetic(x0,y0,z0)
print *, "OK: ecef2geodetic"

call test_enu_ecef(east, north, up, lat, lon, alt, xl, yl, zl)
print *, "OK enu2ecef, ecef2enu"

call test_ecef2aer(xl, yl, zl, lat, lon, alt)
print *, "OK: ecef2aer"

call test_lookatspheroid(lat, lon, alt)
print *, "OK: lookAtSpheroid"

call test_geodetic_aer(lat, lon, alt)
print *, "OK: geodetic2aer, aer2geodetic"

call test_geodetic_enu(lat, lon, alt)
print *, "OK: geodetic2enu, enu2geodetic"

!! ## array tests

call test_array()
print *, "OK: test_array"

!! ## Vallado Tests
call test_vallado(t0)
print *, "OK: vallado"

!! ## Meeus tests

call assert_isclose(anglesep(35.,23., 84.,20.), ha,&
  err_msg='angle_sep')
print *, "OK: Meeus"

print *,'OK: Maptran ',storage_size(0.),' bits'

contains


subroutine test_vallado(t0)

type(datetime), intent(in) :: t0

real, parameter :: azi = 180.1,  eli = 80, jd0 = 2456753.833333

real :: azrd,elrd,rae,dae,jd

jd = toJulian(t0)

!! http://aa.usno.navy.mil/jdconverter?ID=AA&year=2014&month=4&day=6&era=1&hr=8&min=0&sec=0.0

call assert_isclose(jd, jd0, err_msg='toJulian')

call assert_isclose(toGST(jd), 5.4896448816, &
                    rtol=0.1, err_msg='toGST')
call assert_isclose(toLST(radians(-148.), jd), 2.90658, &
                    rtol=0.2, err_msg='toLST')

call azel2radec(azi,eli,lat,lon, jd, rae, dae)
call radec2azel(rae,dae,lat,lon,jd,azrd,elrd)
call assert_isclose([azi,eli],[azrd,elrd],err_msg='azel<->radec')
end subroutine test_vallado


subroutine test_geodetic_enu(lat, lon, alt)
real, intent(in) :: lat, lon, alt
real :: lt, ln, at, e, n, u, ee, nn, uu

call geodetic2enu(lat, lon, alt-1, lat, lon, alt, e, n, u)
call assert_isclose([e,n,u], [0., 0., -1.], atol=atol_dist)

call enu2geodetic(e,n,u, lat,lon,alt, lt, ln, at)
call assert_isclose([lt,ln], [lat,lon])
call assert_isclose(at, alt-1, atol=atol_dist)

call geodetic2enu(radians(lt), radians(ln), at, &
  radians(lat),radians(lon),alt, ee,nn,uu, deg=.false.)
call assert_isclose([ee,nn,uu],[e,n,u], atol=atol_dist)

call enu2geodetic(0.,0., -1., radians(lat),radians(lon), 0., lt,ln,at, deg=.false.)
call assert_isclose([degrees(lt),degrees(ln)], [lat,lon])
call assert_isclose(at, -1., atol=atol_dist)
end subroutine test_geodetic_enu


subroutine test_array()

integer,parameter :: N = 3

real, dimension(N), parameter :: alat = [42,52,62], &
                                    deg0 = [15,30,45], &
                                    aaz = [33,43,53]
real, dimension(N) :: &
ax1, ay1, aaaz1, ax2, ay2,az2, aaaz2, ax3,ay3,aaaz3, &
ae1, an1, au1, ae2,an2,au2, &
alat2,alon2,aalt2, alat3,alon3,aalt3, alat4,alon4,aalt4, &
aaz2, ael2, arng2, aaz3,ael3,arng3, aaz4,ael4,arng4

call assert_isclose(degrees(radians(deg0)), deg0, err_msg='deg<->rad')
print *, "OK: degrees < = > radians"

call geodetic2ecef(alat,lon,alt,ax1,ay1,aaaz1)
call aer2enu(aaz, el, rng, ae1, an1, au1)
call aer2ecef(aaz, el, rng, lat,lon,alt, ax2, ay2, aaaz2)
call ecef2geodetic(ax1,ay1,aaaz1,alat2,alon2,aalt2)
call enu2aer(ae1,an1,au1, aaz2, ael2, arng2)
call ecef2aer(ax2,ay2,az2, lat,lon,alt, aaz3, ael3, arng3)
call aer2geodetic(aaz,el,rng,lat,lon,alt, alat3,alon3,aalt3)
call geodetic2enu(alat3, alon3, aalt3, lat, lon, alt, ae2, an2, au2)
call geodetic2aer(alat3,alon3, aalt3,lat,lon,alt, aaz4, ael4, arng4)
call enu2ecef(ae1,an1,au1,lat,lon,alt, ax3, ay3, aaaz3)
call enu2geodetic(ae2,an2,au2,lat,lon,alt,alat4, alon4, aalt4)

end subroutine test_array


subroutine test_geodetic_aer(lat, lon, alt)
real, intent(in) :: lat, lon, alt
real :: lt, ln, at, a, e, r

real, parameter :: lat1 = 42.0026, lon1 = -81.9978, alt1 = 1.1397e3

call aer2geodetic(0., -90., 1., lat, lon, alt, lt,ln,at)
call assert_isclose([lt,ln], [lat, lon])
call assert_isclose(at, alt-1, atol=atol_dist)

call geodetic2aer(0.,45.,-1000., 0., 45., 0., a, e, r)
call assert_isclose([a,e], [0.,-90.])
call assert_isclose(r, 1000., atol=atol_dist)

call aer2geodetic(radians(az),radians(el),rng, radians(lat),radians(lon),alt, &
  lt,ln,at, deg=.false.)
call assert_isclose([degrees(lt),degrees(ln)], [lat1,lon1])
call assert_isclose(at, alt1, atol=2*atol_dist)

call geodetic2aer(lt,ln, at, radians(lat),radians(lon),alt, a,e,r, deg=.false.)
call assert_isclose([degrees(a),degrees(e)],[az,el], atol=atol_deg)
call assert_isclose(r, rng, atol=atol_dist)

end subroutine test_geodetic_aer


subroutine test_ecef2aer(xl, yl, zl, lat, lon, alt)
real, intent(in) :: xl, yl, zl, lat, lon, alt
real :: a, e, r, x, y,z

call ecef2aer(xl,yl,zl, lat,lon,alt, a, e, r)
call assert_isclose([a,e], [az,el], atol=atol_deg)
call assert_isclose(r, rng, atol=atol_dist, rtol=0.001, &
  err_msg='ecef2aer-degrees1')

call ecef2aer(ea-1., 0., 0., 0., 0., 0., a, e, r)
call assert_isclose([a,e,r], [0., -90., 1.], &
  err_msg='ecef2aer-degrees2')

call ecef2aer(-ea+1., 0., 0., 0., 180., 0., a, e, r)
call assert_isclose([a,e,r], [0., -90., 1.], &
  err_msg='ecef2aer-degrees3')

call ecef2aer(0., ea-1., 0., 0., 90., 0., a, e, r)
call assert_isclose([a,e,r], [0., -90., 1.], &
  err_msg='ecef2aer-degrees4')

call ecef2aer(0., -ea+1., 0., 0., -90., 0., a, e, r)
call assert_isclose([a,e,r], [0., -90., 1.], &
  err_msg='ecef2aer-degrees5')

call ecef2aer(0., 0., eb-1000., 90., 0., 0., a, e, r)
call assert_isclose([a,e,r], [0., -90., 1000.], rtol=0.001)

call ecef2aer(0., 0., -eb+1000., -90., 0., 0., a, e, r)
call assert_isclose([a,e,r], [0., -90., 1000.], rtol=0.001)

call ecef2aer((ea-1000.)/sqrt(2.), (ea-1000.)/sqrt(2.), 0., &
              0., 45., 0., a, e, r)
call assert_isclose([a,e,r], [0., -90., 1000.], &
                      rtol=0.001, err_msg='ecef2aer-degrees')

call aer2ecef(radians(az),radians(el),rng, radians(lat),radians(lon),alt, x,y,z, deg=.false.)
call assert_isclose([x,y,z],[xl,yl,zl])

call ecef2aer(x,y,z, radians(lat),radians(lon),alt, a,e,r, deg=.false.)
call assert_isclose([degrees(a),degrees(e)],[az,el], atol=atol_deg)
call assert_isclose(r, rng, atol=atol_dist)

end subroutine test_ecef2aer


subroutine test_enu_ecef(east, north, up, lat, lon, alt, xl, yl, zl)
real, intent(in) :: east, north, up, lat, lon, alt, xl, yl, zl
real :: x, y, z, e, n, u

call enu2ecef(east, north, up, lat,lon,alt, x, y, z)
call assert_isclose([x,y,z], [xl,yl, zl], err_msg='enu2ecef-degrees')

call ecef2enu(x,y,z, lat,lon,alt, e,n,u)
call assert_isclose([e,n,u], [east, north, up], atol=atol_dist)

call enu2ecef(east, north, up, radians(lat),radians(lon),alt, x,y,z, deg=.false.)
call assert_isclose([x,y,z],[xl,yl,zl], err_msg='enu2ecef: rad')

call ecef2enu(x,y,z, radians(lat),radians(lon),alt, e,n,u, deg=.false.)
call assert_isclose([e,n,u],[east, north, up], atol=atol_dist)

end subroutine test_enu_ecef


subroutine test_ecef2geodetic(x0,y0,z0)
real, intent(in) :: x0, y0, z0
real :: lt, ln, at

call ecef2geodetic(x0,y0,z0,lt,ln,at)
call assert_isclose([lt,ln], [lat,lon])
call assert_isclose(at, alt, atol=atol_dist, &
  err_msg='ecef2geodetic-degrees1')

call ecef2geodetic(ea-1, 0., 0., lt, ln, at)
call assert_isclose([lt,ln,at], [0., 0., -1.], &
  err_msg='ecef2geodetic-degrees2')

call ecef2geodetic(0., ea-1, 0., lt, ln, at)
call assert_isclose([lt,ln,at], [0., 90., -1.], &
  err_msg='ecef2geodetic-degrees3')

call ecef2geodetic(0., -ea+1, 0., lt, ln, at)
call assert_isclose([lt,ln,at], [0., -90., -1.], &
  err_msg='ecef2geodetic-degrees4')

call ecef2geodetic(0., 0., eb-1, lt, ln, at)
call assert_isclose([lt,ln,at], [90., 0., -1.], &
  err_msg='ecef2geodetic-degrees5')

call ecef2geodetic(0., 0., -eb+1, lt, ln, at)
call assert_isclose([lt,ln,at], [-90., 0., -1.], &
  err_msg='ecef2geodetic-degrees6')

call ecef2geodetic(-ea+1, 0., 0., lt, ln, at)
call assert_isclose([lt,ln,at], [0., 180., -1.], &
  err_msg='ecef2geodetic-degrees7')

call ecef2geodetic((ea-1000)/sqrt(2.), &
                    (ea-1000)/sqrt(2.), 0., lt, ln, at)
call assert_isclose([lt,ln,at], [0., 45., -1000.], atol=atol_dist)

end subroutine test_ecef2geodetic


subroutine test_enu2aer(east, north, up, az, el, rng)
real, intent(in) :: east, north, up, az, el, rng
real :: a, e, r, ee, n, u

call enu2aer(east, north, up, a, e, r)
call assert_isclose([a, e, r], [az, el, rng], err_msg='enu2aer-degrees')

call enu2aer(1., 0., 0., a, e, r)
call assert_isclose([a, e, r], [90., 0., 1.])

call aer2enu(radians(az),radians(el),rng, ee,n,u, deg=.false.)
call assert_isclose([ee,n,u], [east, north, up])

call enu2aer(ee,n,u, a, e, r, deg=.false.)
call assert_isclose([degrees(a),degrees(e),r],[az,el,rng])
end subroutine test_enu2aer


subroutine test_aer2enu(az, el, rng, east, north, up)
real, intent(in) :: az, el, rng, east, north, up

real :: e,n,u

call aer2enu(az, el, rng, e, n, u)
call assert_isclose([e, n, u], [east, north, up])

call aer2enu(0., -90., 1., e, n, u)
call assert_isclose([e, n, u], [0., 0., -1.], atol=1e-6)

call aer2enu(0., 90., 1., e, n, u)
call assert_isclose([e, n, u], [0., 0., 1.], atol=1e-6)
end subroutine test_aer2enu


subroutine test_lookatspheroid(lat, lon, alt)
real, intent(in) :: lat, lon, alt
real, dimension(3) :: az, tilt, olat, olon, orng

real:: nan
nan = ieee_value(0., ieee_quiet_nan)
if(.not.ieee_is_nan(nan)) then
  write(stderr,*) "SKIP: lookatspheroid: compiler NaN not working"
  return
endif

az = [0., 10., 125.]
tilt = [30, 45, 90]

call lookAtSpheroid(lat, lon, alt, az, 0., olat, olon, orng)
call assert_isclose(olat, lat, err_msg='lookAtSpheroid:lat0')
call assert_isclose(olon, lon, err_msg='lookAtSpheroid:lon0')
call assert_isclose(orng, alt, rtol=0.01, err_msg='lookAtSpheroid:rng0')

call lookAtSpheroid(lat, lon, alt, az, tilt, olat, olon, orng)
call assert_isclose(olat, [42.00103959, 42.00177328, nan],  equal_nan=.true., err_msg='lookAtSpheroid:lat')
call assert_isclose(olon, [lon, -81.9995808, nan],  equal_nan=.true., err_msg='lookAtSpheroid:lon')
call assert_isclose(orng, [230.9413173, 282.84715651, nan], rtol=0.01, equal_nan=.true., err_msg='lookAtSpheroid:rng')

print *, "OK: lookatspheroid"
end subroutine test_lookatspheroid


subroutine test_geodetic2ecef(lat,lon, alt)
real, intent(in) :: lat, lon, alt
real :: x,y,z, lt, ln, at

call geodetic2ecef(lat,lon,alt, x, y, z)
call assert_isclose([x, y, z],[x0,y0,z0], err_msg='geodetic2ecef-degrees1')

call geodetic2ecef(0., 0., -1., x, y, z)
call assert_isclose([x, y, z], [ea-1, 0., 0.], err_msg='geodetic2ecef-degrees2')

call geodetic2ecef(0., 90., -1., x, y, z)
call assert_isclose([x, y, z], [0., ea-1, 0.], err_msg='geodetic2ecef-degrees3')

call geodetic2ecef(0., -90., -1., x, y, z)
call assert_isclose([x,y,z], [0., -ea+1, 0.], err_msg='geodetic2ecef-degrees4')

call geodetic2ecef(0., -180., -1., x, y, z)
call assert_isclose([x,y,z], [-ea+1, 0., 0.], atol=atol_dist, err_msg='geodetic2ecef-degrees5')
call geodetic2ecef(0., 180., -1., x, y, z)
call assert_isclose([x,y,z], [-ea+1, 0., 0.], atol=atol_dist, err_msg='geodetic2ecef-degrees6')

call geodetic2ecef(90., 0., -1., x, y, z)
call assert_isclose([x,y,z], [0., 0., eb-1], err_msg='geodetic2ecef-degrees7')

call geodetic2ecef(-90., 0., -1., x, y, z)
call assert_isclose([x,y,z], [0., 0., -eb+1], err_msg='geodetic2ecef-degrees8')

call geodetic2ecef(90., 15., -1., x, y, z)
call assert_isclose([x,y,z], [0., 0., eb-1], err_msg='geodetic2ecef-degrees9')

call geodetic2ecef(-45., 45., -1., x, y, z)
call assert_isclose([x,y,z], [3.194418645060574e+06, 3.194418645060574e+06, -4.487347701759138e+06], &
  err_msg='geodetic2ecef-degrees10')

call geodetic2ecef(45., -45., -1., x, y, z)
call assert_isclose([x,y,z], [3.194418645060574e+06, -3.194418645060574e+06, 4.487347701759138e+06], &
  err_msg='geodetic2ecef-degrees1')

call ecef2geodetic(x0,y0,z0, lt,ln,at, deg=.false.)
call assert_isclose([degrees(lt),degrees(ln)],[lat,lon], atol=atol_deg)
call assert_isclose(at, alt, atol=atol_dist)

call geodetic2ecef(radians(lat),radians(lon),alt, x,y,z, deg=.false.)
call assert_isclose([x,y,z],[x0,y0,z0])

end subroutine test_geodetic2ecef


subroutine test_aer2ecef(az, el, rng, lat, lon, alt)
real, intent(in) :: az, el, rng, lat, lon, alt
real :: x,y,z

call aer2ecef(az, el, rng, lat, lon, alt, x, y, z)
call assert_isclose([x, y, z], [xl,yl,zl], err_msg='aer2ecef1')

call aer2ecef(0., -90., 1., 0., 0., -1., x, y, z)
call assert_isclose([x, y, z], [ea-1., 0., 0.], atol=1e-6, err_msg='aer2ecef2')

call aer2ecef(0., -90., 1., 0., 90., -1., x, y, z)
call assert_isclose([x,y,z],[0., ea-1., 0.], atol=1e-6, err_msg='aer2ecef3')

call aer2ecef(0., -90., 1., 90., 0., -1., x, y, z)
call assert_isclose([x,y,z],[0., 0., eb-1.], atol=1e-6, err_msg='aer2ecef4')

call aer2ecef(0., -90., 1., 0., 90., -1., x, y, z)
call assert_isclose([x,y,z],[0., ea-1., 0.], atol=1e-6, err_msg='aer2ecef5')

call aer2ecef(0., -90., 1., 0., 45., -1., x, y, z)
call assert_isclose([x,y,z], &
        [(ea-1.)/sqrt(2.), (ea-1.)/sqrt(2.), 0.], atol=1e-6, err_msg='aer2ecef6')

end subroutine test_aer2ecef

end program
