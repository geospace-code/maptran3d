configure_file(
input : 'maptran.in.f90',
output : 'maptran.f90',
configuration : wp_conf)

configure_file(
input : 'vallado.in.f90',
output : 'vallado.f90',
configuration : wp_conf
)

srcs = files('aer.f90', 'ecef.f90', 'enu.f90', 'sphere.f90', 'utils.f90',
meson.current_build_dir() / 'maptran.f90',
meson.current_build_dir() / 'vallado.f90')
