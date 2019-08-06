project('MapTran', 'fortran',
  meson_version : '>=0.51.2',
  default_options : ['default_library=static', 'buildtype=release'])

realbits = '-DREALBITS='+get_option('realbits')

fc = meson.get_compiler('fortran')
if fc.get_id() == 'gcc'
  add_project_arguments('-fimplicit-none', '-Wall', '-Wextra', '-Wpedantic', language : 'fortran')
  if fc.version().version_compare('>=8')
    add_project_arguments('-std=f2018', language : 'fortran')
  endif
elif fcid == 'intel'
  add_project_arguments('-warn', '-heap-arrays', '-stand f18', language : 'fortran')
elif fcid == 'intel-cl'
  add_project_arguments('/warn', '/heap-arrays', '/stand:f18', language : 'fortran')
elif fcid == 'pgi'
  add_project_arguments('-C', '-Mdclchk', language : 'fortran')
endif

# --- Maptran library
maptran = library('maptran', 'src/maptran.F90', 'src/vallado.F90', 'src/aer.f90', 'src/ecef.f90', 'src/enu.f90',
  fortran_args : realbits)

# --- testing
mtexe = executable('testmaptran', 'tests/test_mod.f90', 'src/assert.F90',
  link_with : maptran,
  fortran_args : realbits)

test('Maptran'+realbits, mtexe,
  timeout: 30)
