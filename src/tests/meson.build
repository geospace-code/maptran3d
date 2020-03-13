mtexe = executable('testmaptran',
  sources: ['test_mod.f90', 'assert.F90'],
  link_with : maptran,
  fortran_args : realbits)

test('Maptran'+realbits, mtexe,
  timeout: 30)

executable('benchmark_maptran', 'benchmark.f90',
  link_with: maptran)
