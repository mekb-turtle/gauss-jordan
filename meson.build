project(
    'gauss-jordan',
    'c',
    version: '0.1',
    license: 'MPL-2.0',
    default_options: ['warning_level=3'],
)

main_name = 'gauss-jordan'

# define source files
src_main = files('src' / 'main.c')
src = files('src'/'gauss-jordan.c', 'src'/'gauss-jordan.h', 'src'/'rational.c', 'src'/'rational.h')

# define project metadata
url = 'https://github.com/mekb-turtle/gauss-jordan'
name = meson.project_name()
version = meson.project_version()

cc = meson.get_compiler('c')
libm = cc.find_library('m', required : false)

dependencies = [libm]

# used if you want to install a library
libdir = get_option('prefix') / get_option('libdir')

# begin meson configuration

add_project_arguments(
    f'-DPROJECT_NAME="@name@"',
    f'-DPROJECT_VERSION="@version@"',
    f'-DPROJECT_URL="@url@"',
    language: 'c',
)

fs = import('fs')

# main executable
exe = executable(
    main_name,
    sources: src + src_main,
    install: true,
    dependencies: dependencies,
    install_rpath: libdir,
)
