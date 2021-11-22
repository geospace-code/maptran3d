configure_file(maptran.in.f90 maptran.f90 @ONLY)
configure_file(vallado.in.f90 vallado.f90 @ONLY)

target_sources(maptran PRIVATE
${CMAKE_CURRENT_BINARY_DIR}/maptran.f90
${CMAKE_CURRENT_BINARY_DIR}/vallado.f90
aer.f90 ecef.f90 enu.f90 sphere.f90 utils.f90
)
