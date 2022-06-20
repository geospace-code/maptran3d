add_executable(test_maptran test_mod.f90 assert.f90)
target_link_libraries(test_maptran PRIVATE maptran)
add_test(NAME unit:maptran COMMAND test_maptran)

add_executable(benchmark_maptran benchmark.f90)
target_link_libraries(benchmark_maptran PRIVATE maptran)

get_property(test_names DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)
set_tests_properties(${test_names} PROPERTIES TIMEOUT 60)

if(BUILD_SHARED_LIBS AND WIN32)
  set_tests_properties(${test_names} PROPERTIES
  ENVIRONMENT_MODIFICATION "PATH=path_list_append:$<TARGET_FILE_DIR:maptran>")
endif()
