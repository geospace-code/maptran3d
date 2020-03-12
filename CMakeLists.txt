cmake_minimum_required(VERSION 3.14)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release")
endif()
project(MapTran
  LANGUAGES Fortran
  VERSION 1.1.1)
enable_testing()

if(NOT realbits)
  set(realbits 64)
endif()

include(cmake/compilers.cmake)

# --- Maptran library
add_library(maptran)
target_compile_definitions(maptran PRIVATE REALBITS=${realbits})
target_include_directories(maptran INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)
set_target_properties(maptran PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
add_library(maptran::maptran ALIAS maptran)
install(TARGETS maptran
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

add_subdirectory(src)

# --- testing

add_executable(test_maptran tests/test_mod.f90 src/assert.F90)
target_link_libraries(test_maptran maptran::maptran)
target_compile_definitions(test_maptran PRIVATE REALBITS=${realbits})
add_test(NAME Maptran COMMAND $<TARGET_FILE:test_maptran>)
set_tests_properties(Maptran PROPERTIES TIMEOUT 30)
