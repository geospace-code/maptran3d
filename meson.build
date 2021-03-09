project('MapTran', 'fortran',
  version: '1.1.4',
  meson_version : '>=0.51.2',
  default_options : ['default_library=static', 'buildtype=release', 'warning_level=2'])

fc = meson.get_compiler('fortran')
if fc.get_id() == 'intel'
  add_project_arguments('-heap-arrays', language : 'fortran')
elif fc.get_id() == 'intel-cl'
  add_project_arguments('/fpp', '/heap-arrays', language : 'fortran')
endif

wp_real = get_option('realbits')=='32' ? 'wp=>real32' : 'wp=>real64'

wp_conf = configuration_data()
wp_conf.set('wp_real', wp_real)

subdir('src')
# --- Maptran library
maptran = library('maptran',
sources: srcs,
install: true)

subdir('src/tests')
