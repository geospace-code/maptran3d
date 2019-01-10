if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_compile_options(-g -O0)
else()
  add_compile_options(-O3)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)

  list(APPEND FFLAGS -traceback -warn -implicitnone)
  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    list(APPEND FFLAGS -check all -debug extended)
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
  list(APPEND FFLAGS -Mallocatable=03)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL NAG)
  list(APPEND FFLAGS -f2008 -C -colour -gline -nan -info -u)
endif()

include(CheckFortranSourceCompiles)

check_fortran_source_compiles("block; end block; end" 
  f08block SRC_EXT f90)
if(NOT f08block)
  message(WARNING "f2008 BLOCK not supported by " ${CMAKE_Fortran_COMPILER_ID})
endif()

check_fortran_source_compiles("use, intrinsic:: ieee_arithmetic; end" 
  f08ieee SRC_EXT f90)
if(NOT f08ieee)
  message(FATAL_ERROR "IEEE_arithmetic not supported by " ${CMAKE_Fortran_COMPILER_ID})
endif()
                       
check_fortran_source_compiles("error stop; end" 
  f08errorstop SRC_EXT f90)
if(NOT f08errorstop)
  set(f08errorstop 0)
endif()
