cmake_minimum_required(VERSION 3.7)
project(MapTran Fortran)
enable_testing()

set(CTEST_OUTPUT_ON_FAILURE ON)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/)

include(cmake/compilers.cmake)
    
if(NOT realbits)
  set(realbits 64)
endif()

if(NOT WIN32)
  option(BUILD_SHARED_LIBS "Shared libraries .so. Feel free to turn off." ON)
endif()

# --- Maptran library
add_library(maptran src/maptran.F90 src/vallado.F90)
target_compile_options(maptran PRIVATE ${FFLAGS})
target_compile_definitions(maptran PUBLIC REALBITS=${realbits})
# --- testing
add_executable(testmaptran tests/test_mod.f90 src/assert.F90)
target_link_libraries(testmaptran PRIVATE maptran ${FLIBS})
target_compile_options(testmaptran PRIVATE ${FFLAGS})
target_compile_definitions(testmaptran PRIVATE F08=${f08errorstop})
add_test(NAME Maptran COMMAND testmaptran)
set_tests_properties(Maptran PROPERTIES TIMEOUT 30)


