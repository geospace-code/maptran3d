cmake_minimum_required(VERSION 3.14)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "default build type")
endif()

project(MapTran
  LANGUAGES Fortran
  VERSION 1.1.2
  HOMEPAGE_URL https://github.com/geospace-code/maptran3d)

enable_testing()
include(CTest)
include(FeatureSummary)

include(cmake/compilers.cmake)

if(realbits EQUAL 32)
  set(wp_real "wp=>real32")
else()
  set(wp_real "wp=>real64")
endif()

# OPTIONAL link-time optimization
if(PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  if(CMAKE_BUILD_TYPE STREQUAL Release)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT lto_ok OUTPUT _err)

    if(lto_ok)
      set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
  endif()
endif()

set(CTEST_TEST_TIMEOUT 30)

# --- Maptran library
add_subdirectory(src)

if(PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  add_subdirectory(src/tests)

  add_feature_info(ipo lto_ok "Interprocess / Link optimization ${_err}")
  feature_summary(WHAT ALL)
endif()
