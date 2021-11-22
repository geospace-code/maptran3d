project('MapTran', 'fortran',
  meson_version : '>=0.57.0',
  default_options : ['buildtype=release'])

fc = meson.get_compiler('fortran')
if fc.get_id() == 'intel'
  add_project_arguments('-heap-arrays', language : 'fortran')
elif fc.get_id() == 'intel-cl'
  add_project_arguments('/fpp', '/heap-arrays', language : 'fortran')
endif

subdir('src')
# --- Maptran library
maptran = library('maptran',
sources: srcs,
install: true)

subdir('src/tests')
