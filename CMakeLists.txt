cmake_minimum_required(VERSION 3.14)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release")
endif()
project(MapTran
  LANGUAGES Fortran
  VERSION 1.1.1
  HOMEPAGE_URL https://github.com/geospace-code/maptran3d)
enable_testing()

if(NOT realbits)
  set(realbits 64)
endif()

include(cmake/compilers.cmake)

# OPTIONAL link-time optimization
if(CMAKE_BUILD_TYPE STREQUAL Release)
	include(CheckIPOSupported)
	check_ipo_supported(RESULT lto_ok OUTPUT _err)

	if(lto_ok)
		message(STATUS "IPO / LTO enabled")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(STATUS "IPO / LTO disabled: ${_err}")
	endif()
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set(CTEST_TEST_TIMEOUT 30)

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
add_subdirectory(src/tests)
