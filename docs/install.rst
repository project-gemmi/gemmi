
Installation
============

.. highlight:: none

C++ library
-----------

Before version 0.6 gemmi was a header-only library.
Many functions are still in headers. If you use only such function,
you only need to ensure that the ``include`` directory is in your
include path when compiling your program. For example::

    git clone https://github.com/project-gemmi/gemmi.git
    c++ -Igemmi/include -O2 my_program.cpp

Otherwise, you either need to build gemmi_cpp library,
or add (selected) files from src/ to your project.

If you use **CMake**, you may

* use find_package for installed gemmi::

    find_package(gemmi 0.6.4 CONFIG REQUIRED)

* or add gemmi as a git submodule and use add_subdirectory::

    add_subdirectory(gemmi EXCLUDE_FROM_ALL)

* or use FetchContent::

    add_subdirectory(gemmi EXCLUDE_FROM_ALL)
    include(FetchContent)
    FetchContent_Declare(
      gemmi
      GIT_REPOSITORY https://github.com/project-gemmi/gemmi.git
      GIT_TAG        ...
    )
    FetchContent_GetProperties(gemmi)
    if (NOT gemmi_POPULATED)
      FetchContent_Populate(gemmi)
      add_subdirectory(${gemmi_SOURCE_DIR} ${gemmi_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif()

Then, to find headers and link your target with the library, use::

    target_link_libraries(example PRIVATE gemmi::gemmi_cpp)

If only headers are needed, do::

    target_link_libraries(example PRIVATE gemmi::headers)

The gemmi::headers interface, which is also included in gemmi::gemmi_cpp,
adds two things: include dictory and *compile feature* cxx_std_11 (a minimal
requirement for the compilation).

----

Note on Unicode: if a file name is passed to Gemmi (through ``std::string``)
it is assumed to be in ASCII or UTF-8.

.. _install_py:

Python module
-------------

From PyPI
~~~~~~~~~

To install the gemmi module do::

    pip install gemmi

We have binary wheels for several Python versions (for all supported CPython
versions and one PyPy version), so the command usually downloads binaries.
If a matching wheel is not available,
the module is compiled from source -- it takes several minutes
and requires a C++ compiler.

Other binaries
~~~~~~~~~~~~~~

If you use the `CCP4 suite <https://www.ccp4.ac.uk/>`_,
you can find gemmi there.

If you use Anaconda Python, you can install
`package conda <https://github.com/conda-forge/gemmi-feedstock>`_
from conda-forge::

    conda install -c conda-forge gemmi

These distribution channels may have an older version of gemmi.

From git
~~~~~~~~

The latest version can be installed directly from the repository.
Either use::

    pip install git+https://github.com/project-gemmi/gemmi.git

or clone the `project <https://github.com/project-gemmi/gemmi/>`_
(or download a zip file) and from the top-level directory do::

    pip install .

On Windows Python should automatically find an appropriate compiler (MSVC).
If the compiler is not installed, pip shows a message with a download link.

If gemmi is already installed, uninstall the old version first
(``pip uninstall``) or add option ``--upgrade``.

Alternatively, you can build a cloned project with CMake::

    cmake -D USE_PYTHON=1 .
    make -j4 py

Fortran and C bindings
----------------------

The Fortran bindings are in early stage and are not documented yet.
They use the ISO_C_BINDING module introduced in Fortran 2003
and `shroud <https://github.com/LLNL/shroud>`_.
You may see the ``fortran/`` directory to know what to expect.
This directory contains Makefile -- run make to built the bindings.
(They are currently not integrated with the cmake build.)

..
 The bindings and usage examples can be compiled with CMake::

    cmake -D USE_FORTRAN=1 .
    make

The C bindings are used only for making Fortran bindings,
but they should be usable on their own.

..
 If you use cmake to build the project
 you get a static library ``libcgemmi.a`` that can be used from C,
 together with the :file:`fortran/*.h` headers.

Program
-------

The library comes with a command-line program also named ``gemmi``.

Binaries
~~~~~~~~

Binaries are distributed with the CCP4 suite and with Global Phasing software.
They are also in `PyPI <https://pypi.org/project/gemmi-program/>`_
(``pip install gemmi-program``) and
`conda-forge packages <https://anaconda.org/conda-forge/gemmi/files>`_.

The very latest builds (as well as a little older ones)
can be downloaded from CI jobs:

- for Windows --
  click the first (green) job in
  `AppVeyor CI <https://ci.appveyor.com/project/wojdyr/gemmi>`_
  and find gemmi.exe in the Artifacts tab,
- for Linux and Mac -- sign in to GitHub (no special permissions are needed,
  but GitHub requires sign-in for artifacts),
  click the first job (with ✅) in
  `GitHub Actions <https://github.com/project-gemmi/gemmi/actions/workflows/ci.yml>`_
  and download a zip file from the Artifacts section.

From source
~~~~~~~~~~~

To build it from source, first make sure you have git, cmake and C++ compiler
installed (on Ubuntu: ``sudo apt install git cmake make g++``), then::

    git clone https://github.com/project-gemmi/gemmi.git
    cd gemmi
    cmake .
    make

Testing
-------

The main automated tests are in Python::

    python3 -m unittest discover -v tests/

We also have doctest tests in the documentation, and some others.
All of them can be run from the ``run-tests.sh`` script in the repository.

Credits
-------

This project is using code from a number of third-party open-source projects.

Projects used in the C++ library, included under
``include/gemmi/third_party/`` (if used in headers) or ``third_party/``:

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

Projects included under ``third_party/`` that are not used in the library
itself, but are used in command-line utilities, python bindings or tests:

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
* `cctbx <https://github.com/cctbx/cctbx_project>`_ -- used in tests
  (if cctbx is not present, these tests are skipped) and
  in scripts that generated space group data and 2-fold twinning operations.
  License: 3-clause BSD.

Email me if I forgot about something.

List of C++ headers
-------------------

Here is a list of C++ headers in ``gemmi/include/``.
This list also gives an overview of the library.

.. include:: headers.rst
