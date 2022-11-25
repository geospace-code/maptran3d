cmake_minimum_required(VERSION 3.14...3.25)

project(maptran
LANGUAGES Fortran
VERSION 1.2.0
HOMEPAGE_URL https://github.com/geospace-code/maptran3d
)

enable_testing()

option(BUILD_TESTING "self tests" on)

if(NOT realbits)
  set(realbits 32)
endif()

include(cmake/compilers.cmake)

# --- Maptran library
add_library(maptran)
target_include_directories(maptran PUBLIC
$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
$<INSTALL_INTERFACE:include>
)
set_target_properties(maptran PROPERTIES
Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/maptran.mod TYPE INCLUDE)
install(TARGETS maptran EXPORT ${PROJECT_NAME}-targets)

add_subdirectory(src)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()

include(cmake/install.cmake)
