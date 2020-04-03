set(CMAKE_CONFIGURATION_TYPES "Release;RelWithDebInfo;Debug" CACHE STRING "Build type selections" FORCE)

if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  if(WIN32)
  string(APPEND CMAKE_Fortran_FLAGS " /arch:native /stand:f18 /traceback /warn /heap-arrays")
else()
  string(APPEND CMAKE_Fortran_FLAGS " -march=native -stand f18 -traceback -warn -heap-arrays")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  string(APPEND CMAKE_Fortran_FLAGS " -Wall -Wextra -fimplicit-none")
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -fcheck=all -Werror=array-bounds")
  # -march=native is not for all CPU arches with GCC.
  add_compile_options(-mtune=native)

  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 8)
    string(APPEND CMAKE_Fortran_FLAGS " -std=f2018")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
  string(APPEND CMAKE_Fortran_FLAGS " -C -Mdclchk")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Flang)
  string(APPEND CMAKE_Fortran_FLAGS " -W")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL NAG)
  string(APPEND CMAKE_Fortran_FLAGS " -f2018 -C -colour -gline -nan -info -u")
endif()

include(CheckFortranSourceCompiles)
check_fortran_source_compiles("implicit none (external); end" f2018impnone SRC_EXT f90)
if(NOT f2018impnone)
  message(FATAL_ERROR "Compiler does not support Fortran 2018 IMPLICIT NONE (EXTERNAL): ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
endif()

include(CheckFortranSourceRuns)
check_fortran_source_runs("use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
real :: r
r = ieee_value(1., ieee_quiet_nan)
if (.not.ieee_is_nan(r)) error stop
end program" f03nan)
