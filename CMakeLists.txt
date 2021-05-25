cmake_minimum_required(VERSION 3.14...3.21)

file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION PROJECT_VERSION
  REGEX "^([0-9]+\.[0-9]+\.[0-9]+)" LIMIT_INPUT 16 LENGTH_MAXIMUM 16 LIMIT_COUNT 1)

project(maptran
  LANGUAGES Fortran
  VERSION ${PROJECT_VERSION}
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

include(cmake/install.cmake)
