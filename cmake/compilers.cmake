if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  add_compile_options(
  -warn -traceback
  "$<$<CONFIG:Debug,RelWithDebInfo>:-check>"
  "$<IF:$<EQUAL:${realbits},32>,,-r8>"
  )
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  add_compile_options(
  -Wall -fimplicit-none
  "$<$<CONFIG:Debug,RelWithDebInfo>:-fcheck=all;-Werror=array-bounds>"
  "$<IF:$<EQUAL:${realbits},32>,,-fdefault-real-8>"
  )
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL NAG)
  add_compile_options(-f2018 -C -colour -gline -nan -info -u)
endif()
