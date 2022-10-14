# This file is based on https://github.com/pybind/python_example
# which is under a BSD-like license:
# https://github.com/pybind/python_example/blob/master/LICENSE

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from setuptools import distutils
import os
import sys

USE_SYSTEM_ZLIB = False
MIN_PYBIND_VER = '2.6.2'

def read_version_from_header():
    with open('include/gemmi/version.hpp') as f:
        for line in f:
            if line.startswith('#define GEMMI_VERSION '):
                return line.split()[2].strip('"dev')

__version__ = read_version_from_header()

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __str__(self):
        import pybind11
        return pybind11.get_include()

if USE_SYSTEM_ZLIB:
    zlib_library = 'z'
    zlib_include_dirs = []
    build_libs = []
else:
    zlib_library = 'gemmi_zlib'
    zlib_include_dirs = ['third_party/zlib']
    zlib_files = ['third_party/zlib/%s.c' % name for name in
                  ['adler32', 'crc32', 'gzlib', 'gzread', 'inflate',
                   'inftrees', 'inffast', 'zutil']]
    zlib_macros = [('NO_GZCOMPRESS', '1')]
    if os.name != 'nt':
        zlib_macros += [('Z_HAVE_UNISTD_H', '1')]
    build_libs = [('gemmi_zlib', {'sources': zlib_files,
                                  'macros': zlib_macros})]

ext_modules = [
    Extension('gemmi',
              ['python/%s.cpp' % name for name in
                  ['gemmi', 'align', 'ccp4', 'chemcomp', 'cif', 'custom',
                   'elem', 'hkl', 'grid', 'meta', 'mol', 'monlib', 'mtz',
                   'read', 'recgrid', 'scaling', 'search', 'sf', 'sym',
                   'topo', 'unitcell', 'write']],
              include_dirs=zlib_include_dirs + [
                  'include',
                  'third_party',
                  # Path to pybind11 headers
                  get_pybind_include(),
              ],
              libraries=[zlib_library],
              language='c++'),
]


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp', delete=False) as f:
        # Don't trigger -Wunused-parameter.
        f.write('int main (int, char **) { return 0; }')
        fname = f.name
    try:
        compiler.compile([fname], extra_postargs=[flagname])
    except distutils.errors.CompileError:
        return False
    finally:
        try:
            os.remove(fname)
        except OSError:
            pass
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14/17] compiler flag.

    The newer version is prefered over c++11 (when it is available).
    """
    flags = ['-std=c++20', '-std=c++17', '-std=c++14', '-std=c++11']

    # C++17 on Mac requires higher -mmacosx-version-min, skip it for now
    if sys.platform == 'darwin':
        flags = flags[2:]

    for flag in flags:
        if has_flag(compiler, flag):
            return flag

    raise RuntimeError('Unsupported compiler -- at least C++11 support '
                       'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc', '/D_CRT_SECURE_NO_WARNINGS'],
        'unix': [],
    }
    l_opts = {
        'msvc': [],
        'unix': [],
    }

    if sys.platform == 'win32':
        if sys.version_info[0] == 2:
            # without these variables distutils insist on using VS 2008
            os.environ['DISTUTILS_USE_SDK'] = '1'
            os.environ['MSSdk'] = '1'
        if sys.version_info[0] >= 3:
            c_opts['msvc'].append('/D_UNICODE')

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])

        if sys.platform == 'darwin':
            darwin_opts = []
            if 'MACOSX_DEPLOYMENT_TARGET' not in os.environ:
                import platform
                mac_ver = platform.mac_ver()
                current_macos = tuple(int(x) for x in mac_ver[0].split(".")[:2])
                if current_macos > (10, 9):
                    darwin_opts.append('-mmacosx-version-min=10.9')
            if has_flag(self.compiler, '-stdlib=libc++'):
                darwin_opts.append('-stdlib=libc++')
            opts += darwin_opts
            link_opts += darwin_opts

        if ct == 'unix':
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
            if has_flag(self.compiler, '-g0'):
                opts.append('-g0')
            if has_flag(self.compiler, '-Wl,-s'):
                link_opts.append('-Wl,-s')
        elif ct.startswith('mingw'):
            #opts.append('-std=c++14')
            opts.append(cpp_flag(self.compiler))
            opts.append('-fvisibility=hidden')
            opts.append('-g0')
            link_opts.append('-Wl,-s')
        for ext in self.extensions:
            ext.define_macros = [('VERSION_INFO',
                                  '"%s"' % self.distribution.get_version())]
            ext.extra_compile_args = opts
            ext.extra_link_args = link_opts
        build_ext.build_extensions(self)

setup(
    name='gemmi',
    version=__version__,
    author='Marcin Wojdyr',
    author_email='wojdyr@gmail.com',
    url='https://project-gemmi.github.io/',
    description='library for structural biology',
    long_description='''\
    Library for macromolecular crystallography and structural bioinformatics.
    For working with coordinate files (mmCIF, PDB, mmJSON),
    refinement restraints (monomer library), electron density maps (CCP4),
    and crystallographic reflection data (MTZ, SF-mmCIF). It understands
    crystallographic symmetries, it knows how to switch between the real
    and reciprocal space and it can do a few other things.

    The setup.py script builds only Python extension.
    Use cmake to build also a command-line program.
    ''',
    long_description_content_type='text/plain',
    libraries=build_libs,
    ext_modules=ext_modules,
    packages=['gemmi-examples'],
    package_dir={'gemmi-examples': 'examples'},
    install_requires=[],
    setup_requires=['pybind11>=' + MIN_PYBIND_VER],
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,
    license='MPL-2.0',  # or, at your option, LGPL-3.0
    keywords=('structural bioinformatics, structural biology, crystallography,'
              ' CIF, mmCIF, PDB, CCP4, MTZ'),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: C++',
        'Programming Language :: Python',
    ],
)
