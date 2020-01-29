cmake_minimum_required(VERSION 3.12)  # 3.12 for Fortran submodule
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release")
endif()
project(MapTran
  LANGUAGES Fortran
  VERSION 1.1.0)
enable_testing()

if(NOT realbits)
  set(realbits 64)
endif()

include(cmake/compilers.cmake)

# --- Maptran library
add_library(maptran src/maptran.F90 src/aer.f90 src/ecef.f90 src/enu.f90 src/vallado.F90)
target_compile_definitions(maptran PRIVATE REALBITS=${realbits})
target_include_directories(maptran INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(maptran PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# --- testing

add_executable(testmaptran tests/test_mod.f90 src/assert.F90)
target_link_libraries(testmaptran maptran)
add_test(NAME Maptran COMMAND testmaptran)
set_tests_properties(Maptran PROPERTIES TIMEOUT 30)
