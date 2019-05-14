
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)

if(NOT WIN32)
  set(FFLAGS -stand f18 -implicitnone -traceback -warn -heap-arrays)
else()
  set(FFLAGS /stand:f18 /4Yd /traceback /warn /heap-arrays)
  # Note: -g is /debug:full for ifort Windows
endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 8)
    list(APPEND FFLAGS -std=f2018)
  endif()

  list(APPEND FFLAGS -fimplicit-none -march=native -Wall -Wextra -Wpedantic -Warray-bounds)# -Wfatal-errors)
  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    list(APPEND FFLAGS -ffpe-trap=overflow)
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Flang)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL NAG)
  list(APPEND FFLAGS -f2008 -C -colour -gline -nan -info -u)
endif()

include(CheckFortranSourceCompiles)

check_fortran_source_compiles("use, intrinsic:: ieee_arithmetic; end"
  f08ieee SRC_EXT f90)
if(NOT f08ieee)
  message(FATAL_ERROR "IEEE_arithmetic not supported by " ${CMAKE_Fortran_COMPILER_ID})
endif()
