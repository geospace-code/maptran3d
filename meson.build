project('MapTran', 'fortran',
  version: '1.1.0',
  meson_version : '>=0.51.2',
  default_options : ['default_library=static', 'buildtype=release', 'warning_level=3', 'fortran_std=f2018'])

realbits = '-DREALBITS=' + get_option('realbits')

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
  fortran_args : realbits)

# --- testing
if fc.links('block; end block; end', name: 'F2008 Block syntax')
  mtexe = executable('testmaptran', 'tests/test_mod.f90', 'src/assert.F90',
    link_with : maptran,
    fortran_args : realbits)

  test('Maptran'+realbits, mtexe,
    timeout: 30)
endif