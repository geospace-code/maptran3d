project('MapTran', 'fortran',
  meson_version : '>=0.51.2',
  default_options : ['default_library=static', 'buildtype=release', 'warning_level=3'])

realbits = '-DREALBITS='+get_option('realbits')

fc = meson.get_compiler('fortran')
f18flag = fc.first_supported_argument(['-std=f2018', '-stand f18', '/stand:f18'])
if fc.get_id() == 'gcc'
  add_project_arguments('-fimplicit-none', f18flag, language : 'fortran')
elif fc.get_id() == 'intel'
  add_project_arguments('-warn', '-heap-arrays', f18flag, language : 'fortran')
elif fc.get_id() == 'intel-cl'
  add_project_arguments('/fpp', '/warn', '/heap-arrays', f18flag, language : 'fortran')
elif fc.get_id() == 'pgi'
  add_project_arguments('-C', '-Mdclchk', language : 'fortran')
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