module vallado
!! based on http://www.smad.com/vallado/

implicit none (type, external)
private
public:: toGST, toJulian, toLST, radec2azel, azel2radec

real, parameter :: pi = 4 * atan(1.0)

type,public :: datetime
integer :: year, month, day, hour, minute
real :: second
end type


contains

elemental subroutine azel2radec(azz,ell,llat,llon,jd, ra,decl)
! convert azimuth, elevation to right ascension, declination
!
! inputs
! ------
! az, el: azimuth, elevation (degrees)
! lat, lon: geodetic latitude, longitude (degrees)
! jd: julian date (decimal)
!
! outputs
! -------
! ra, decl: right ascension, declination (degrees)

real, intent(in) :: azz,ell,llat,llon
real, intent(in) :: jd ! Julian Date
real, intent(out) :: ra, decl

real :: lst, lha, sinv, cosv, az, el, lat, lon

az = radians(azz)
el = radians(ell)
lat = radians(llat)
lon = radians(llon)

Decl = ASIN(SIN(El)*SIN(lat) + &
            COS(el)*COS(lat)*COS(Az) )

Sinv = -(SIN(az)*COS(el)*COS(lat)) / &
        (COS(lat)*COS(Decl))

Cosv = (SIN(el) - SIN(lat)*SIN(decl)) / &
       (COS(lat)*COS(Decl))

LHA  = ATAN2(Sinv,Cosv)
lst = toLST(Lon, JD)

ra = modulo(degrees(LST - LHA), 360.)
decl = degrees(decl)

end subroutine azel2radec


elemental SUBROUTINE radec2azel(rra,dDecl,llat,llon,jd, Az,El)
! convert right ascension, declination to azimuth, elevation
!
! inputs
! ------
! ra, decl: right ascension, declination (degrees)
! lat, lon: geodetic latitude, longitude (degrees)
! jd: julian date (decimal)
!
! outputs
! -------
! az, el: azimuth, elevation (degrees)

REAL, intent(in) :: rra,dDecl, llat, llon
real, intent(in) :: jd
real, intent(out) :: Az,El

REAL :: Sinv, Cosv, LHA, lat, lon, ra, decl

lat = radians(llat)
lon = radians(llon)
ra = radians(rra)
decl = radians(ddecl)

LHA = toLST(Lon, JD) - ra

El = ASIN( SIN(Decl)*SIN(lat) + &
            COS(Decl)*COS(lat)*COS(LHA) )

Sinv = -SIN(LHA)*COS(Decl)*COS(lat) / &
      (COS(el)*COS(lat))

Cosv = ( SIN(Decl)-SIN(el)*SIN(lat) ) / &
       (COS(el)*COS(lat))

Az = modulo(degrees(ATAN2(Sinv,Cosv)), 360.)
el = degrees(el)


END subroutine radec2azel


elemental real function toLST(Lon, JD) result(LST)
!! Julian Date => local sidereal time
!!
!! ## inputs
!!
!! *  lon: geodetic longitude (radians)
!! *  jd: Julian Date (decimal)

REAL, intent(in) :: Lon, JD

LST = Lon + toGST(jd)

LST = modulo(LST, 2*pi )

END function toLST


elemental real function toJulian(t) result(jd)
! Gregorian date, time => Julian Date
!
! inputs
! ------
! time: Gregorian user-defined type
!
! output
! -----
! JD: Julian Date

type(datetime), intent(in) :: t
real :: B, y, m

y = t%year
m = t%month

IF ( M <= 2 ) THEN
  Y = y - 1
  M = m + 12
ENDIF

B = 2 - INT(Y*0.01) + INT(INT(Y*0.01)*0.25)

JD= INT( 365.25*(Y + 4716) ) + &
    INT( 30.6001*(M+1) ) + &
    t%day + B - 1524.5 + &
    ( (t%second/60.0 + t%minute ) / 60.0 + t%hour ) / 24.0

END function toJulian


elemental real FUNCTION toGST(JD) result(GST)
! Julian Date => to Greenwich Sidereal Time
!
! inputs
! ------
! JD: Julian Date (decimal)
!
! output
! ------
! GST: Greenwich Sidereal Time (decimal)

real, intent(in) :: JD
real :: TUT1


TUT1= ( JD - 2451545. ) / 36525.
gst = -6.2e-6*TUT1**3 + &
            0.093104*TUT1**2 + &
           (876600.*3600. + 8640184.812866)*TUT1 + &
            67310.54841

gst = modulo(radians(gst) / 240., 2*pi) ! 360/86400 = 1/240, to deg, to rad

end function toGST


elemental real function degrees(rad)
real, intent(in) :: rad

degrees = 180 / pi * rad
end function degrees


elemental real function radians(deg)
real, intent(in) :: deg

radians = pi / 180 * deg
end function radians


end module vallado
