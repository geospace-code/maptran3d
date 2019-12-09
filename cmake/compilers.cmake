
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)

if(NOT WIN32)
  set(FFLAGS -stand f18 -warn declarations -traceback -warn -heap-arrays)
else()
  set(FFLAGS /stand:f18 /warn:declarations /traceback /warn /heap-arrays)
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
  set(FFLAGS -C -Mdclchk)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Flang)
  set(CFLAGS -W)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL NAG)
  list(APPEND FFLAGS -f2008 -C -colour -gline -nan -info -u)
endif()

include(CheckFortranSourceCompiles)
