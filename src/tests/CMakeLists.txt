add_executable(test_maptran)
target_sources(test_maptran PRIVATE test_mod.f90 assert.f90)
target_link_libraries(test_maptran maptran::maptran)
add_test(NAME unit:maptran COMMAND $<TARGET_FILE:test_maptran>)
set_target_properties(test_maptran PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

add_executable(benchmark_maptran benchmark.f90)
target_link_libraries(benchmark_maptran maptran::maptran)
set_target_properties(benchmark_maptran PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
