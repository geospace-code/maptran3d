configure_file(
input : 'assert.in.f90',
output : 'assert.f90',
configuration : wp_conf
)

mtexe = executable('testmaptran',
sources: ['test_mod.f90', meson.current_build_dir() / 'assert.f90'],
link_with : maptran)

test('Maptran' + get_option('realbits'),
mtexe,
timeout: 30)

executable('benchmark_maptran', 'benchmark.f90',
link_with: maptran)
