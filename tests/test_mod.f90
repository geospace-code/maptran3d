!! error tolerances are set for single precision, double precision is much more precise.
use, intrinsic:: iso_fortran_env, only: compiler_version
use, intrinsic:: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
use maptran
use vallado
use assert, only: assert_allclose

implicit none

type(Ellipsoid), parameter :: spheroid = wgs84Ellipsoid

real(wp), parameter :: &
lat = 42, lon= -82, alt = 200, &
az = 33, el=70, rng= 1e3_wp, &
x0 = 660.6752518e3_wp, y0 = -4700.9486832e3_wp, z0 = 4245.7376622e3_wp, &
xl = 660.930e3_wp, yl = -4701.424e3_wp, zl = 4246.579e3_wp, & !< aer2ecef
er = 186.277521_wp, nr = 286.84222_wp, ur = 939.69262_wp, & !< aer2enu
lat1 = 42.0026_wp, lon1 = -81.9978_wp, alt1 = 1.1397e3_wp, & !< aer2geodetic
azi = 180.1_wp,  eli = 80._wp, &
ha = 45.482789587392013_wp
               
integer,parameter :: N = 3               
real(wp), dimension(N), parameter :: alat = [42,52,62], &
                                     deg0 = [15,30,45], &
                                     aaz = [33,43,53]
                                     
type(datetime), parameter :: t0 = datetime(2014,4,6,8,0,0) !< UTC 
real(wp), parameter ::  jd0 = 2456753.833333_wp


real(wp) :: azrd,elrd,rae,dae,jd, lat5(3), lon5(3), rng5(3), az5(3), tilt(3), ea, eb, atol_dist, atol_deg


real(wp), dimension(N) :: &
ax1, ay1, aaaz1, ax2, ay2,az2, aaaz2, ax3,ay3,aaaz3, &
ae1, an1, au1, ae2,an2,au2, &
alat2,alon2,aalt2, alat3,alon3,aalt3, alat4,alon4,aalt4, &
aaz2, ael2, arng2, aaz3,ael3,arng3, aaz4,ael4,arng4         

real(wp):: nan

select case (storage_size(nan))
  case (32)
    atol_dist = 1._wp !< 1 meter
    atol_deg = 0.1_wp
  case (64)
    atol_dist = 0.001_wp !< 1 mm
    atol_deg = 1e-3_wp
  case (128)
    atol_dist = 1e-6_wp
    atol_deg = 1e-8_wp
end select

nan = ieee_value(0._wp, ieee_quiet_nan)

ea = spheroid%SemimajorAxis
eb = spheroid%SemiminorAxis

print *, compiler_version()

!! ## scalar degrees
     
az5 = [0._wp, 10._wp, 125._wp];
tilt = [30, 45, 90];

call lookAtSpheroid(lat, lon, alt, az5, 0._wp, lat5, lon5, rng5)
call assert_allclose(lat5, lat) 
call assert_allclose(lon5, lon) 
call assert_allclose(rng5, alt, rtol=0.01_wp, err_msg='lookAtSpheroid')

call lookAtSpheroid(lat, lon, alt, az5, tilt, lat5, lon5, rng5)


call assert_allclose(lat5, [42.00103959_wp, 42.00177328_wp, nan],  equal_nan=.true.)
call assert_allclose(lon5, [lon, -81.9995808_wp, nan],  equal_nan=.true.)
call assert_allclose(rng5, [230.9413173_wp, 282.84715651_wp, nan], rtol=0.01_wp, equal_nan=.true.)


test_geodetic2ecef: block
  real(wp) :: x,y,z, lt, ln, at
  
  call geodetic2ecef(lat,lon,alt, x, y, z)
  call assert_allclose([x, y, z],[x0,y0,z0], &
                     err_msg='geodetic2ecef-degrees')
                     
  call geodetic2ecef(0._wp, 0._wp, -1._wp, x, y, z)
  call assert_allclose([x, y, z], [ea-1, 0._wp, 0._wp])

  call geodetic2ecef(0._wp, 90._wp, -1._wp, x, y, z)
  call assert_allclose([x, y, z], [0._wp, ea-1, 0._wp])

  call geodetic2ecef(0._wp, -90._wp, -1._wp, x, y, z)
  call assert_allclose([x,y,z], [0._wp, -ea+1, 0._wp])

  call geodetic2ecef(0._wp, -180._wp, -1._wp, x, y, z)
  call assert_allclose([x,y,z], [-ea+1, 0._wp, 0._wp], atol=atol_dist)
  call geodetic2ecef(0._wp, 180._wp, -1._wp, x, y, z)
  call assert_allclose([x,y,z], [-ea+1, 0._wp, 0._wp], atol=atol_dist)

  call geodetic2ecef(90._wp, 0._wp, -1._wp, x, y, z)
  call assert_allclose([x,y,z], [0._wp, 0._wp, eb-1])

  call geodetic2ecef(-90._wp, 0._wp, -1._wp, x, y, z)
  call assert_allclose([x,y,z], [0._wp, 0._wp, -eb+1])

  call geodetic2ecef(90._wp, 15._wp, -1._wp, x, y, z)
  call assert_allclose([x,y,z], [0._wp, 0._wp, eb-1])

  call geodetic2ecef(-45._wp, 45._wp, -1._wp, x, y, z)
  call assert_allclose([x,y,z], [3.194418645060574e+06_wp, 3.194418645060574e+06_wp, -4.487347701759138e+06_wp])

  call geodetic2ecef(45._wp, -45._wp, -1._wp, x, y, z)
  call assert_allclose([x,y,z], [3.194418645060574e+06_wp, -3.194418645060574e+06_wp, 4.487347701759138e+06_wp])
  
  call ecef2geodetic(x0,y0,z0, lt,ln,at, deg=.false.)
  call assert_allclose([degrees(lt),degrees(ln)],[lat,lon], atol=atol_deg)
  call assert_allclose(at, alt, atol=atol_dist)
  
  call geodetic2ecef(radians(lat),radians(lon),alt, x,y,z, deg=.false.)
  call assert_allclose([x,y,z],[x0,y0,z0])
end block test_geodetic2ecef


test_aer2enu: block
  real(wp) :: e,n,u

  call aer2enu(az, el, rng, e, n, u)
  call assert_allclose([e, n, u], [er, nr, ur])

  call aer2enu(0._wp, -90._wp, 1._wp, e, n, u)
  call assert_allclose([e, n, u], [0._wp, 0._wp, -1._wp], atol=1e-6_wp)

  call aer2enu(0._wp, 90._wp, 1._wp, e, n, u)
  call assert_allclose([e, n, u], [0._wp, 0._wp, 1._wp], atol=1e-6_wp)
end block test_aer2enu


test_aer2ecef: block
  real(wp) :: x,y,z

  call aer2ecef(az, el, rng, lat, lon, alt, x, y, z)
  call assert_allclose([x, y, z], [xl,yl,zl])

  call aer2ecef(0._wp, -90._wp, 1._wp, 0._wp, 0._wp, -1._wp, x, y, z)
  call assert_allclose([x, y, z], [ea-1._wp, 0._wp, 0._wp], atol=1e-6_wp)

  call aer2ecef(0._wp, -90._wp, 1._wp, 0._wp, 90._wp, -1._wp, x, y, z)
  call assert_allclose([x,y,z],[0._wp, ea-1._wp, 0._wp], atol=1e-6_wp)

  call aer2ecef(0._wp, -90._wp, 1._wp, 90._wp, 0._wp, -1._wp, x, y, z)
  call assert_allclose([x,y,z],[0._wp, 0._wp, eb-1._wp], atol=1e-6_wp)

  call aer2ecef(0._wp, -90._wp, 1._wp, 0._wp, 90._wp, -1._wp, x, y, z)
  call assert_allclose([x,y,z],[0._wp, ea-1._wp, 0._wp], atol=1e-6_wp)

  call aer2ecef(0._wp, -90._wp, 1._wp, 0._wp, 45._wp, -1._wp, x, y, z)
  call assert_allclose([x,y,z], &
         [(ea-1._wp)/sqrt(2._wp), (ea-1._wp)/sqrt(2._wp), 0._wp], atol=1e-6_wp)
end block test_aer2ecef


test_ecef2geodetic: block
  real(wp) :: lt, ln, at

  call ecef2geodetic(x0,y0,z0,lt,ln,at)
  call assert_allclose([lt,ln], [lat,lon])
  call assert_allclose(at, alt, atol=atol_dist)

  call ecef2geodetic(ea-1, 0._wp, 0._wp, lt, ln, at)
  call assert_allclose([lt,ln,at], [0._wp, 0._wp, -1._wp], &
                       err_msg='ecef2geodetic-degrees')

  call ecef2geodetic(0._wp, ea-1, 0._wp, lt, ln, at)
  call assert_allclose([lt,ln,at], [0._wp, 90._wp, -1._wp], &
                       err_msg='ecef2geodetic-degrees')

  call ecef2geodetic(0._wp, -ea+1, 0._wp, lt, ln, at)
  call assert_allclose([lt,ln,at], [0._wp, -90._wp, -1._wp], &
                       err_msg='ecef2geodetic-degrees')

  call ecef2geodetic(0._wp, 0._wp, eb-1, lt, ln, at)
  call assert_allclose([lt,ln,at], [90._wp, 0._wp, -1._wp], &
                      err_msg='ecef2geodetic-degrees')

  call ecef2geodetic(0._wp, 0._wp, -eb+1, lt, ln, at)
  call assert_allclose([lt,ln,at], [-90._wp, 0._wp, -1._wp], &
                      err_msg='ecef2geodetic-degrees')

  call ecef2geodetic(-ea+1, 0._wp, 0._wp, lt, ln, at)
  call assert_allclose([lt,ln,at], [0._wp, 180._wp, -1._wp], &
                      err_msg='ecef2geodetic-degrees')

  call ecef2geodetic((ea-1000)/sqrt(2._wp), &
                     (ea-1000)/sqrt(2._wp), 0._wp, lt, ln, at)
  call assert_allclose([lt,ln,at], [0._wp, 45._wp, -1000._wp], atol=atol_dist)
  
end block test_ecef2geodetic


test_enu2aer: block
  real(wp) :: a, e, r, ee, n, u

  call enu2aer(er,nr,ur, a, e, r)
  call assert_allclose([a, e, r], [az, el, rng], err_msg='enu2aer-degrees')

  call enu2aer(1._wp, 0._wp, 0._wp, a, e, r)
  call assert_allclose([a, e, r], [90._wp, 0._wp, 1._wp])
  
  call aer2enu(radians(az),radians(el),rng, ee,n,u, deg=.false.)
  call assert_allclose([ee,n,u], [er,nr,ur])

  call enu2aer(ee,n,u, a, e, r, deg=.false.)
  call assert_allclose([degrees(a),degrees(e),r],[az,el,rng])
end block test_enu2aer


test_ecef2aer: block
  real(wp) :: a, e, r, x, y,z

  call ecef2aer(xl,yl,zl, lat,lon,alt, a, e, r)
  call assert_allclose([a,e], [az,el], atol=atol_deg, rtol=0.001_wp)
  call assert_allclose(r, rng, atol=atol_dist, rtol=0.001_wp)

  call ecef2aer(ea-1._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, a, e, r)
  call assert_allclose([a,e,r], [0._wp, -90._wp, 1._wp], &
                       err_msg='ecef2aer-degrees')

  call ecef2aer(-ea+1._wp, 0._wp, 0._wp, 0._wp, 180._wp, 0._wp, a, e, r)
  call assert_allclose([a,e,r], [0._wp, -90._wp, 1._wp], &
                       err_msg='ecef2aer-degrees')

  call ecef2aer(0._wp, ea-1._wp, 0._wp, 0._wp, 90._wp, 0._wp, a, e, r)
  call assert_allclose([a,e,r], [0._wp, -90._wp, 1._wp], &
                       err_msg='ecef2aer-degrees')
                       
  call ecef2aer(0._wp, -ea+1._wp, 0._wp, 0._wp, -90._wp, 0._wp, a, e, r)
  call assert_allclose([a,e,r], [0._wp, -90._wp, 1._wp], &
                       err_msg='ecef2aer-degrees')

  call ecef2aer(0._wp, 0._wp, ea-1._wp, 90._wp, 0._wp, 0._wp, a, e, r)
  call assert_allclose([a,e,r], [0._wp, -90._wp, 1._wp], &
                       rtol=1e6_wp, err_msg='ecef2aer-degrees')
  call ecef2aer(0._wp, 0._wp, -ea+1._wp, -90._wp, 0._wp, 0._wp, a, e, r)
  call assert_allclose([a,e,r], [0._wp, -90._wp, 1._wp], &
                       rtol=1e6_wp, err_msg='ecef2aer-degrees')

  call ecef2aer((ea-1000._wp)/sqrt(2._wp), (ea-1000._wp)/sqrt(2._wp), 0._wp, &
                0._wp, 45._wp, 0._wp, a, e, r)
  call assert_allclose([a,e,r], [0._wp, -90._wp, 1000._wp], &
                       rtol=0.001_wp, err_msg='ecef2aer-degrees')
                       
 call aer2ecef(radians(az),radians(el),rng, radians(lat),radians(lon),alt, x,y,z, deg=.false.)
  call assert_allclose([x,y,z],[xl,yl,zl])

  call ecef2aer(x,y,z, radians(lat),radians(lon),alt, a,e,r, deg=.false.)
  call assert_allclose([degrees(a),degrees(e)],[az,el], atol=atol_deg)
  call assert_allclose(r, rng, atol=atol_dist)

end block test_ecef2aer


geodetic_aer: block
  real(wp) :: lt, ln, at, a, e, r

  call aer2geodetic(0._wp, -90._wp, 1._wp, lat, lon, alt, lt,ln,at)
  call assert_allclose([lt,ln], [lat, lon])
  call assert_allclose(at, alt-1, atol=atol_dist)
  
  
  call geodetic2aer(0._wp,45._wp,-1000._wp, 0._wp, 45._wp, 0._wp, a, e, r)
  call assert_allclose([a,e], [0._wp,-90._wp])
  call assert_allclose(r, 1000._wp, atol=atol_dist)
  
  call aer2geodetic(radians(az),radians(el),rng, radians(lat),radians(lon),alt, &
    lt,ln,at, deg=.false.)
  call assert_allclose([degrees(lt),degrees(ln)], [lat1,lon1])
  call assert_allclose(at, alt1, atol=atol_dist)
  
  call geodetic2aer(lt,ln, at, radians(lat),radians(lon),alt, a,e,r, deg=.false.)
  call assert_allclose([degrees(a),degrees(e)],[az,el], atol=atol_deg)
  call assert_allclose(r, rng, atol=atol_dist)

end block geodetic_aer


geodetic_enu: block
  real(wp) :: lt, ln, at, e, n, u, ee, nn, uu

  call geodetic2enu(lat, lon, alt-1, lat, lon, alt, e, n, u)
  call assert_allclose([e,n,u], [0._wp, 0._wp, -1._wp], atol=atol_dist)

  call enu2geodetic(e,n,u, lat,lon,alt, lt, ln, at)
  call assert_allclose([lt,ln], [lat,lon])
  call assert_allclose(at, alt-1, atol=atol_dist)
  
  call geodetic2enu(radians(lt), radians(ln), at, &
    radians(lat),radians(lon),alt, ee,nn,uu, deg=.false.)
  call assert_allclose([ee,nn,uu],[e,n,u], atol=atol_dist)
  
  call enu2geodetic(0._wp,0._wp, -1._wp, radians(lat),radians(lon), 0._wp, lt,ln,at, deg=.false.)
  call assert_allclose([degrees(lt),degrees(ln)], [lat,lon])
  call assert_allclose(at, -1._wp, atol=atol_dist)
end block geodetic_enu


enu_ecef: block
  real(wp) :: x, y, z, e, n, u

  call enu2ecef(er,nr,ur, lat,lon,alt, x, y, z)
  call assert_allclose([x,y,z], [xl,yl, zl], err_msg='enu2ecef-degrees')

  call ecef2enu(x,y,z, lat,lon,alt, e,n,u)
  call assert_allclose([e,n,u], [er,nr,ur], atol=atol_dist)
  
  call enu2ecef(er,nr,ur, radians(lat),radians(lon),alt, x,y,z, deg=.false.)
  call assert_allclose([x,y,z],[xl,yl,zl], err_msg='enu2ecef: rad')

  call ecef2enu(x,y,z, radians(lat),radians(lon),alt, e,n,u, deg=.false.)
  call assert_allclose([e,n,u],[er,nr,ur], atol=atol_dist)

end block enu_ecef


call assert_allclose(degrees(radians(deg0)), deg0, err_msg='deg<->rad')


!! ## array tests

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

!! ## Vallado Tests
jd = toJulian(t0)

!! http://aa.usno.navy.mil/jdconverter?ID=AA&year=2014&month=4&day=6&era=1&hr=8&min=0&sec=0.0

call assert_allclose(jd, jd0, err_msg='toJulian')

call assert_allclose(toGST(jd), 5.4896448816_wp, &
                    rtol=0.1_wp, err_msg='toGST')
call assert_allclose(toLST(radians(-148._wp), jd), 2.90658_wp, &
                    rtol=0.2_wp, err_msg='toLST')

call azel2radec(azi,eli,lat,lon, jd, rae, dae)
call radec2azel(rae,dae,lat,lon,jd,azrd,elrd)
call assert_allclose([azi,eli],[azrd,elrd],err_msg='azel<->radec')

!! ## Meeus tests

call assert_allclose(anglesep(35._wp,23._wp, 84._wp,20._wp), ha,& 
                   err_msg='angle_sep')


print *,'OK: Maptran ',storage_size(0._wp),' bits'

end program
