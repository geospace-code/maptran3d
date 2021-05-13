configure_file(maptran.in.f90 maptran.f90 @ONLY)
configure_file(vallado.in.f90 vallado.f90 @ONLY)

add_library(maptran)
target_sources(maptran PRIVATE
${CMAKE_CURRENT_BINARY_DIR}/maptran.f90 ${CMAKE_CURRENT_BINARY_DIR}/vallado.f90
aer.f90 ecef.f90 enu.f90 sphere.f90 utils.f90)

target_include_directories(maptran
  PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)

set_target_properties(maptran PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/maptran.mod
  DESTINATION include)

install(TARGETS maptran EXPORT ${PROJECT_NAME}Targets)
