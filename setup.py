# This file is based on https://github.com/pybind/python_example
# which is under a BSD-like license:
# https://github.com/pybind/python_example/blob/master/LICENSE

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from setuptools import distutils
import os
import sys

USE_SYSTEM_ZLIB = False
MIN_PYBIND_VER = '2.4.3'

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

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        # We should have pybind11 installed by now, but old pip
        # before https://github.com/pypa/pip/pull/2616
        # would not handle dependencies properly.
        # So in such case at least hint a workaround to the user.
        try:
            import pybind11
        except ImportError:
            print('\n' + 50*'*')
            print('*****  Please try to install pybind11 first  *****')
            print(50*'*' + '\n')
            sys.exit(1)
        if pybind11.__version__ < MIN_PYBIND_VER:
            print('\n' + 50*'*')
            print('Use pybind11 >= %s.' % MIN_PYBIND_VER)
            print('You have pybind11 %s in %s'
                  % (pybind11.__version__, os.path.dirname(pybind11.__file__)))
            sys.exit(1)
        return pybind11.get_include(self.user)

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
                  ['gemmi', 'cif', 'hkl', 'grid', 'mol', 'read',
                   'monlib', 'sym', 'unitcell', 'write']],
              include_dirs=zlib_include_dirs + [
                  'include',
                  'third_party',
                  # Path to pybind11 headers
                  get_pybind_include(),
                  get_pybind_include(user=True)
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
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.

    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc', '/D_CRT_SECURE_NO_WARNINGS'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    elif sys.platform == 'win32' and sys.version_info[0] == 2:
        # without these variables distutils insist on using VS 2008
        os.environ['DISTUTILS_USE_SDK'] = '1'
        os.environ['MSSdk'] = '1'

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' %
                        self.distribution.get_version())
        elif ct.startswith('mingw'):
            # avoid has_flag() - it didn't work for us in msys2/mingw setup
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append('-std=c++14')
            opts.append('-fvisibility=hidden')
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)

setup(
    name='gemmi',
    version=__version__,
    author='Marcin Wojdyr',
    author_email='wojdyr@gmail.com',
    url='https://project-gemmi.github.io/',
    description='General MacroMolecular I/O',
    long_description='''\
    Library for macromolecular crystallography and structural bioinformatics.
    For working with coordinate files (mmCIF, PDB, mmJSON),
    refinement restraints (monomer library), electron density maps
    and crystallographic reflections.''',
    libraries=build_libs,
    ext_modules=ext_modules,
    packages=['gemmi-examples'],
    package_dir={'gemmi-examples': 'examples'},
    install_requires=['pybind11>=' + MIN_PYBIND_VER],
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,
    license='MPL-2.0',
    keywords=('structural bioinformatics, structural biology, crystallography,'
              ' CIF, mmCIF, PDB, CCP4'),
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
