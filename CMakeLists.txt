cmake_minimum_required(VERSION 3.7)
project(MapTran Fortran)
enable_testing()

set(CTEST_OUTPUT_ON_FAILURE ON)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/)

include(compilers.cmake)
    

# --- assert library
add_library(assert src/assert.F90)
target_compile_options(assert PRIVATE ${FFLAGS})
target_link_libraries(assert PRIVATE ${FLIBS})
if(f08errorstop)
  set_target_properties(assert PROPERTIES 
                        COMPILE_DEFINITIONS F08)
endif()

# --- other libraries
add_library(maptran src/maptran.f90 src/vallado.f90)
target_link_libraries(maptran PRIVATE assert ${FLIBS})
target_compile_options(maptran PRIVATE ${FFLAGS})
# --- testing
add_executable(testmaptran tests/test_mod.f90)
target_link_libraries(testmaptran PRIVATE assert maptran ${FLIBS})
target_compile_options(testmaptran PRIVATE ${FFLAGS})
add_test(NAME Maptran COMMAND testmaptran)
set_tests_properties(Maptran PROPERTIES TIMEOUT 30)


