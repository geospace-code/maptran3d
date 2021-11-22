cmake_minimum_required(VERSION 3.19...3.22)

file(READ ${CMAKE_CURRENT_LIST_DIR}/codemeta.json _j)
string(JSON PROJECT_VERSION GET ${_j} version)

project(maptran
LANGUAGES Fortran
VERSION ${PROJECT_VERSION}
HOMEPAGE_URL https://github.com/geospace-code/maptran3d
)

include(CTest)

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
  add_subdirectory(src/tests)
endif()

include(cmake/install.cmake)
