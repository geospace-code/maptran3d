cmake_minimum_required(VERSION 3.7)
project(MapTran Fortran)
enable_testing()

set(CTEST_OUTPUT_ON_FAILURE ON)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_compile_options(-g -O0)
else()
  add_compile_options(-O3)
endif()

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  list(APPEND FFLAGS -check all -traceback -warn -debug extended)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 8)
    list(APPEND FFLAGS -std=f2018)
  endif()

  list(APPEND FFLAGS -march=native -ffpe-trap=zero,overflow -fbacktrace 
                    -Wall -Wextra -Wpedantic -Warray-bounds)# -Wfatal-errors)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)
  list(APPEND FFLAGS -Mallocatable=03)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Flang)
  list(APPEND FFLAGS -Mallocatable=03)
  link_libraries(-static-flang-libs)
endif()
    

# --- assert library
add_library(assert src/assert.f90)
target_compile_options(assert PRIVATE ${FFLAGS})

set(f08 GNU Intel)
if(CMAKE_Fortran_COMPILER_ID IN_LIST f08)
 target_sources(assert PRIVATE src/error2008.f90)
else()  # non-Fortran 2008 compilers
  target_sources(assert PRIVATE src/error2003.f90)
endif()
# --- other libraries
add_library(maptran src/maptran.f90 src/vallado.f90)
target_link_libraries(maptran PRIVATE assert)
target_compile_options(maptran PRIVATE ${FFLAGS})
# --- testing
add_executable(testmaptran tests/test_mod.f90)
target_link_libraries(testmaptran PRIVATE assert maptran)
target_compile_options(testmaptran PRIVATE ${FFLAGS})
add_test(NAME Maptran COMMAND testmaptran)
set_tests_properties(Maptran PROPERTIES TIMEOUT 30)

