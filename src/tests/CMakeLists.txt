add_executable(test_maptran test_mod.f90 assert.f90)
target_link_libraries(test_maptran PRIVATE maptran)
add_test(NAME unit:maptran COMMAND $<TARGET_FILE:test_maptran>)

add_executable(benchmark_maptran benchmark.f90)
target_link_libraries(benchmark_maptran PRIVATE maptran)
