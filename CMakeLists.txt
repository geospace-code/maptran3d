cmake_minimum_required(VERSION 3.14)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "default build type")
endif()

project(MapTran
  LANGUAGES Fortran
  VERSION 1.1.2
  HOMEPAGE_URL https://github.com/geospace-code/maptran3d)

include(CTest)

include(cmake/compilers.cmake)

if(realbits EQUAL 32)
  set(wp_real "wp=>real32")
else()
  set(wp_real "wp=>real64")
endif()

# --- Maptran library
add_subdirectory(src)

if(BUILD_TESTING)
  add_subdirectory(src/tests)
endif()
