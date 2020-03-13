program benchmark_maptran
!! Benchmark speed of array
use, intrinsic:: iso_fortran_env, only: int64, real64
use maptran, only: wp, geodetic2ecef

implicit none

integer(int64) :: tic, toc, rate
real(real64) :: time
integer :: N, i
real(wp), allocatable, dimension(:,:,:) :: lat, lon, alt, x, y, z
character(8) :: buf

N = 100
call get_command_argument(1, buf, status=i)
if (i==0) read(buf,*) N
allocate(lat(N,N,N))
allocate(lon, alt, x, y, z, mold=lat)

call random_number(lat)
call random_number(lon)
call random_number(alt)

call system_clock(count_rate=rate) !< set timer to highest precision

call system_clock(tic)
call geodetic2ecef(lat, lon, alt, x, y, z)
call system_clock(toc)

time = (toc-tic) / real(rate, wp)

print '(F7.1,A,I8,A,I4)',time*1000,' ms for ',N,'^3 array geodetic2ecef,  real bits:',storage_size(z)

end program
