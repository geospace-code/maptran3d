project('MapTran', 'fortran', default_options : ['default_library=static'])

REALBITS = '-DREALBITS='+get_option('realbits')

fc = meson.get_compiler('fortran')
if fc.get_id() == 'gcc'
  add_global_arguments('-fimplicit-none', '-Wall', '-Wextra', '-Wpedantic', language : 'fortran')
endif

# --- Maptran library
maptran = library('maptran', 'src/maptran.F90', 'src/vallado.F90',
  fortran_args : REALBITS)
  
# --- testing
mtexe = executable('testmaptran', 'tests/test_mod.f90', 'src/assert.F90',
  link_with : maptran,
  fortran_args : [REALBITS, '-DF08=1'])

test('Maptran', mtexe)
