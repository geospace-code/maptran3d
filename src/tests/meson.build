mtexe = executable('testmaptran',
sources: ['test_mod.f90', 'assert.f90'],
link_with : maptran)

test('Maptran' + get_option('realbits'), mtexe, timeout: 30)

executable('benchmark_maptran', 'benchmark.f90', link_with: maptran)
