cmake_minimum_required(VERSION 3.7)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release")
endif()
project(MapTran Fortran)
enable_testing()

include(cmake/compilers.cmake)

if(NOT realbits)
  set(realbits 64)
endif()

# --- Maptran library
add_library(maptran src/maptran.F90 src/vallado.F90)
target_compile_options(maptran PRIVATE ${FFLAGS})
target_compile_definitions(maptran PUBLIC REALBITS=${realbits})


# --- testing
check_fortran_source_compiles("block; end block; end" f08block SRC_EXT f90)

if(f08block)
add_executable(testmaptran tests/test_mod.f90 src/assert.F90)
target_link_libraries(testmaptran PRIVATE maptran ${FLIBS})
target_compile_options(testmaptran PRIVATE ${FFLAGS})
add_test(NAME Maptran COMMAND testmaptran)
set_tests_properties(Maptran PROPERTIES TIMEOUT 30)
endif()
