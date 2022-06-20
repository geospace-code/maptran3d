project('MapTran', 'fortran',
meson_version : '>=0.57.0'
)

subdir('src')
# --- Maptran library
maptran = library('maptran',
sources: srcs,
install: true)

subdir('test')
