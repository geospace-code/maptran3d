configure_file(assert.in.f90 assert.f90)

add_executable(test_maptran test_mod.f90 ${CMAKE_CURRENT_BINARY_DIR}/assert.f90)
target_link_libraries(test_maptran maptran::maptran)
target_compile_definitions(test_maptran PRIVATE REALBITS=${realbits})
add_test(NAME unit:maptran COMMAND $<TARGET_FILE:test_maptran>)

add_executable(benchmark_maptran benchmark.f90)
target_link_libraries(benchmark_maptran maptran::maptran)
