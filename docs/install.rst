
Installation
============

.. highlight:: none

C++ library
-----------

It is a header-only library. You need to ensure that
the ``include`` directory is in your include path
when compiling your program. For example::

    git clone https://github.com/project-gemmi/gemmi.git
    c++ -std=c++11 -Igemmi/include -O2 my_program.cpp

If you want Gemmi to uncompress gzipped files on the fly
(i.e. if you ``#include <gemmi/gz.hpp>``)
you will also need to link your program with the zlib library.

If a file name is passed to Gemmi (through ``std::string``)
it is assumed to be in ASCII or UTF-8.

.. _install_py:

Python 2.7/3.x module
---------------------

From source
~~~~~~~~~~~

To install the gemmi module you need pip, git and not too old
C++ compiler (GCC 4.8+, Clang 3.4+, MSVC 2015+, ICC 17+)::

    pip install gemmi

(We have binary wheels for some Python versions on Windows, so the command
above may actually download binaries).

Alternatively, to install the latest version directly from the repository::

    pip install git+https://github.com/project-gemmi/gemmi.git

or clone the `project <https://github.com/project-gemmi/gemmi/>`_
(or download a zip file) and from the top-level directory do::

    pip install .

(Setuptools compile only one unit at a time and the whole process
will take several minutes. To make it faster, use
``cmake -D USE_PYTHON=1 .``).

If gemmi is already installed, uninstall the old version first
(``pip uninstall``) or add option ``--upgrade``.

On Windows Python 3.5+ should automatically find an appropriate compiler
(MSVC 2015+) . If the compiler is not installed, pip shows a message
with a download link.
For Python 2.7 pip prefers MSVC 2008, which is too old to compile gemmi.
You may still use MSVC 2015, 2017 or 2019, but before invoking pip you need to
set the compiler environment with one of these commands::

    "%VS140COMNTOOLS%\..\..\VC\vcvarsall.bat" x64
    "%VS140COMNTOOLS%\..\..\VC\vcvarsall.bat"

If you'd like to use PyPy instead of CPython -- PyPy2.7 >= 5.7 is supported
(although only occasionally tested -- open an issue if it doesn't work).

Binaries
~~~~~~~~

If you use the `CCP4 suite <https://www.ccp4.ac.uk/>`_,
you can find gemmi there.

If you use Anaconda Python, you can install
`package conda <https://github.com/conda-forge/gemmi-feedstock>`_
from conda-forge::

    conda install -c conda-forge gemmi

These distribution channels may have a previous version of gemmi.

Fortran and C bindings
----------------------

The Fortran bindings are in early stage and are not documented yet.
They use the ISO_C_BINDING module introduced in Fortran 2003.
You may see the ``fortran/`` directory to know what to expect.
The bindings and usage examples can be compiled with CMake::

    cmake -D USE_FORTRAN=1 .
    make

The C bindings are used only for making Fortran bindings,
but they should be usable on their own.
If you use cmake to build the project
you get a static library ``libcgemmi.a`` that can be used from C,
together with the :file:`fortran/*.h` headers.

Gemmi program
-------------

The library comes with a command-line program. To build it from source,
first make sure you have git, cmake and C++ compiler installed
(on Ubuntu: ``sudo apt install git cmake make g++``), then::

    git clone https://github.com/project-gemmi/gemmi.git
    cd gemmi
    cmake .
    make

Binaries are distributed with the CCP4 suite and with Global Phasing software.
They are also included in
`conda-forge packages <https://anaconda.org/conda-forge/gemmi/files>`_.
Additionally, the very latest Windows builds (as well as older ones)
can be downloaded from
`AppVeyor CI <https://ci.appveyor.com/project/wojdyr/gemmi>`_: click
the first (green) job, then in the Artifacts tab you should find gemmi.exe.

Credits
-------

This project is using code from a number of third-party open-source projects.

Projects used in the C++ library and included under
``include/gemmi/third_party/``:

* `PEGTL <https://github.com/taocpp/PEGTL/>`_ -- library for creating PEG
  parsers. License: MIT.
* `sajson <https://github.com/chadaustin/sajson>`_ -- high-performance
  JSON parser. License: MIT.
* `PocketFFT <https://gitlab.mpcdf.mpg.de/mtr/pocketfft>`_ -- FFT library.
  License: 3-clause BSD.
* `stb_sprintf <https://github.com/nothings/stb>`_ -- locale-independent
  snprintf() implementation. License: Public Domain.
* `fast_float <https://github.com/fastfloat/fast_float>`_ -- locale-independent
  number parsing. License: Apache 2.0.
* `tinydir <https://github.com/cxong/tinydir>`_ -- directory (filesystem)
  reader. License: 2-clause BSD.

Code derived from the following projects is used in the library:

* `ksw2 <https://github.com/lh3/ksw2>`_ -- sequence alignment in
  ``seqalign.hpp`` is based on the ksw_gg function from ksw2. License: MIT.
* `QCProt <https://theobald.brandeis.edu/qcp/>`_ -- superposition method
  in ``qcp.hpp`` is taken from QCProt and adapted to our project. License: BSD.
* `Larch <https://github.com/xraypy/xraylarch>`_ -- calculation of f' and f"
  in ``fprime.hpp`` is based on CromerLiberman code from Larch.
  License: 2-clause BSD.

Projects included under ``third_party/``, not used in the library itself,
but used in command-line utilities, python bindings or tests:

* `The Lean Mean C++ Option Parser <http://optionparser.sourceforge.net/>`_ --
  command-line option parser. License: MIT.
* `doctest <https://github.com/onqtam/doctest>`_ -- testing framework.
  License: MIT.
* `linalg.h <http://github.com/sgorsten/linalg/>`_ -- linear algebra library.
  License: Public Domain.
* `zlib <https://github.com/madler/zlib>`_ -- a subset of the zlib library
  for uncompressing gz files, used as a fallback when the zlib library
  is not found in the system. License: zlib.

Not distributed with Gemmi:

* `pybind11 <https://github.com/pybind/pybind11>`_ -- used for creating
  Python bindings. License: 3-clause BSD.
* `cctbx <https://github.com/cctbx/cctbx_project>`_ -- used in tests and
  in scripts that generated space group data. License: 3-clause BSD.

Email me if I forgot about something.

List of C++ headers
-------------------

Here is a list of C++ headers in ``gemmi/include/``.
This list also gives an overview of the library.

.. include:: headers.rst
