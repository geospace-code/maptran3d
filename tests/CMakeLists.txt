add_executable(test_maptran test_mod.f90 assert.f90)
target_link_libraries(test_maptran PRIVATE maptran)
add_test(NAME unit:maptran COMMAND test_maptran WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

add_executable(benchmark_maptran benchmark.f90)
target_link_libraries(benchmark_maptran PRIVATE maptran)
