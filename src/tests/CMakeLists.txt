add_executable(test_maptran test_mod.f90 assert.F90)
target_link_libraries(test_maptran maptran::maptran)
target_compile_definitions(test_maptran PRIVATE REALBITS=${realbits})
add_test(NAME unit:maptran COMMAND $<TARGET_FILE:test_maptran>)

if(NOT PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  return()
endif()

add_executable(benchmark_maptran benchmark.f90)
target_link_libraries(benchmark_maptran maptran::maptran)
