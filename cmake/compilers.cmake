set(CMAKE_CONFIGURATION_TYPES "Release;RelWithDebInfo;Debug" CACHE STRING "Build type selections" FORCE)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  add_compile_options(
  $<IF:$<BOOL:${WIN32}>,/QxHost,-xHost>
  "$<$<COMPILE_LANGUAGE:Fortran>:-warn;-traceback;-heap-arrays>"
  )
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  add_compile_options(-mtune=native -Wall -Wextra
  "$<$<COMPILE_LANGUAGE:Fortran>:-fimplicit-none>"
  "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fcheck=all;-Werror=array-bounds>"
  )
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL NAG)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-f2018;-C;-colour;-gline;-nan;-info;-u>")
endif()
