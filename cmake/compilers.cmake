if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  if(NOT WIN32)
    string(APPEND CMAKE_Fortran_FLAGS " -march=native -stand f18 -traceback -warn -heap-arrays")
  else()
    string(APPEND CMAKE_Fortran_FLAGS " /arch:native /stand:f18 /traceback /warn /heap-arrays")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 8)
    string(APPEND CMAKE_Fortran_FLAGS " -std=f2018")
  endif()
  # not mtune=native for CPU compatibility
  string(APPEND CMAKE_Fortran_FLAGS " -mtune=native -fimplicit-none -Wall -Wextra")
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -ffpe-trap=overflow -Warray-bounds")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
  string(APPEND CMAKE_Fortran_FLAGS " -C -Mdclchk")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL NAG)
  string(APPEND CMAKE_Fortran_FLAGS " -f2018 -C -colour -gline -nan -info -u")
endif()

include(CheckFortranSourceRuns)
check_fortran_source_runs("use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
real :: r
r = ieee_value(1., ieee_quiet_nan)
if (.not.ieee_is_nan(r)) error stop
end program" f03nan)