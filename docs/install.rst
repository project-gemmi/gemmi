
Installation
============

.. highlight:: none

C++ library
-----------

Gemmi used to be a header-only library (until ver. 0.6.0).
Parts of the library (finding symmetry operations, parsing CIF grammar)
are still header-only; if you happen to use only these parts,
just ensure that Gemmi's `include` directory is in
your project's include path. For example::

    git clone https://github.com/project-gemmi/gemmi.git
    c++ -Igemmi/include -O2 my_program.cpp

However, in most cases, you need to build a library called gemmi_cpp
and link your project against it.

If you use **CMake**, you may:

* first install gemmi and then use find_package::

    find_package(gemmi 0.7.0 CONFIG REQUIRED)

* or add gemmi as a git submodule and use add_subdirectory::

    add_subdirectory(gemmi EXCLUDE_FROM_ALL)

* or use FetchContent::

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

Then link your target with the library (this also takes care of includes)::

    target_link_libraries(example PRIVATE gemmi::gemmi_cpp)

If a target only needs gemmi headers, do this instead::

    target_link_libraries(example PRIVATE gemmi::headers)

The gemmi::headers interface, which is also included in gemmi::gemmi_cpp,
adds two things: the include directory and the *compile feature* cxx_std_14
(a minimal requirement for compilation).

Gemmi can be compiled with either zlib or zlib-ng.
The only difference is that zlib-ng is faster.
Here are the relevant CMake options:

* FETCH_ZLIB_NG -- download, build statically, and use zlib-ng.
* USE_ZLIB_NG -- find zlib-ng installed on the system.
* INTERNAL_ZLIB -- compile third_party/zlib (a subset of zlib distributed
  with gemmi).
* None of the above -- find zlib installed on the system;
  if not found, use third_party/zlib.

On Windows, when a program or library is linked with a zlib(-ng) DLL,
it may require the DLL to be in the same directory.
It is simpler to build zlib-ng statically or use `-D FETCH_ZLIB_NG=ON`.

----

Note on Unicode: if a file name is passed to Gemmi (through `std::string`),
it is assumed to be in ASCII or UTF-8.

.. _install_py:

Python module
-------------

From PyPI
~~~~~~~~~

To install the gemmi module, run::

    pip install gemmi

We have binary wheels for several Python versions (for all supported CPython
versions and one PyPy version), so the command usually downloads binaries.
If a matching wheel is not available,
the module is compiled from source -- it takes a few minutes
and requires a C++ compiler that supports C++17.

Gemmi 0.7+ supports only Python 3.8+.

Other binaries
~~~~~~~~~~~~~~

You can find gemmi:

If you use the `CCP4 suite <https://www.ccp4.ac.uk/>`_,
you can find gemmi there.

If you use conda,
the `gemmi package <https://github.com/conda-forge/gemmi-feedstock>`_,
which includes also a command-line program and C++ dev files,
can be installed from conda-forge::

    conda install -c conda-forge gemmi

These distribution channels may have an older version of gemmi.

From git
~~~~~~~~

The latest version can be installed directly from the repository.
Either use::

    pip install git+https://github.com/project-gemmi/gemmi.git

or clone the `project <https://github.com/project-gemmi/gemmi/>`_
(or download a zip file) and from the top-level directory run::

    pip install .

Building with pip uses scikit-build-core and CMake underneath.
You can pass options to CMake either using the `--config-settings` option
in recent pip versions::

  pip install . --config-settings="cmake.args=-DFETCH_ZLIB_NG=ON"

or by using environment variables such as `CMAKE_ARGS`. See
`scikit-build-core docs <https://scikit-build-core.readthedocs.io/en/latest/configuration.html#configuring-cmake-arguments-and-defines>`_
for details.

If gemmi is already installed, uninstall the old version first
(`pip uninstall`) or add the `--upgrade` option.

Alternatively, you can manually install nanobind and cmake (using pip)
and build a cloned project directly with CMake::

    cmake -D USE_PYTHON=1 .
    make -j4 gemmi_py

Fortran and C bindings
----------------------

The Fortran bindings are in an early stage and are not documented yet.
They use the ISO_C_BINDING module introduced in Fortran 2003
and `shroud <https://github.com/LLNL/shroud>`_.
You can check the `fortran/` directory to see what to expect.
This directory contains a Makefile -- run make to build the bindings.
(They are currently not integrated with the CMake build.)

..
 The bindings and usage examples can be compiled with CMake::

    cmake -D USE_FORTRAN=1 .
    make

The C bindings are used only for making Fortran bindings,
but they should be usable on their own.

..
 If you use cmake to build the project
 you get a static library `libcgemmi.a` that can be used from C,
 together with the :file:`fortran/*.h` headers.

Program
-------

The library comes with a command-line program also named `gemmi`.

Binaries
~~~~~~~~

Binaries are distributed with the CCP4 suite and with Global Phasing software.
They are also in `PyPI <https://pypi.org/project/gemmi-program/>`_
(`pip install gemmi-program`),
`conda-forge packages <https://anaconda.org/conda-forge/gemmi/files>`_,
and a few Linux (and FreeBSD)
`distros <https://repology.org/project/gemmi/versions>`_.

The very latest builds (as well as a little older ones)
can be downloaded from CI jobs:

- For Windows --
  click the first (green) job in
  `AppVeyor CI <https://ci.appveyor.com/project/wojdyr/gemmi>`_
  and find gemmi.exe in the Artifacts tab (if there is also a dll file there,
  it's a dynamically linked build and both files are needed).
- For Linux and Mac -- sign in to GitHub (no special permissions are needed,
  but GitHub requires sign-in for artifacts), go to gemmi's
  `gemmi's CI workflow <https://github.com/project-gemmi/gemmi/actions/workflows/ci.yml>`_,
  click the latest job with ✅, scroll to the bottom of the page,
  and download one of the zip files from the Artifacts section.

From source
~~~~~~~~~~~

To build it from source, first make sure you have git, cmake and C++ compiler
installed (on Ubuntu: `sudo apt install git cmake make g++`), then::

    git clone https://github.com/project-gemmi/gemmi.git
    cd gemmi
    cmake .
    make

Alternatively, you can use `pip install git+https://...`, which installs
both the Python module and the program. If you are not using the Python module,
you can use pip to build only the program::

    pip install git+https://github.com/project-gemmi/gemmi.git --config-settings=cmake.args=-DONLY_PROGRAM=ON

Testing
-------

The main automated tests are in Python::

    python3 -m unittest discover -v tests/

We also have *Python doctest* tests in the documentation,
and a few other test routines.
All the commands used for testing are listed in the `run-tests.sh`
script in the repository.

Credits
-------

This project is using code from a number of third-party open-source projects.

Projects used in the C++ library, included under
`include/gemmi/third_party/` (if used in headers) or `third_party/`:

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
  `seqalign.hpp` is based on the ksw_gg function from ksw2. License: MIT.
* `QCProt <https://theobald.brandeis.edu/qcp/>`_ -- superposition method
  in `qcp.hpp` is taken from QCProt and adapted to our project. License: BSD.
* `Larch <https://github.com/xraypy/xraylarch>`_ -- calculation of f' and f"
  in `fprime.cpp` is based on CromerLiberman code from Larch.
  License: 2-clause BSD.

Projects included under `third_party/` that are not used in the library
itself, but are used in command-line utilities, python bindings or tests:

* `zpp serializer <https://github.com/eyalz800/serializer>`_ --
  serialization framework. License: MIT.
* `The Lean Mean C++ Option Parser <http://optionparser.sourceforge.net/>`_ --
  command-line option parser. License: MIT.
* `doctest <https://github.com/onqtam/doctest>`_ -- testing framework.
  License: MIT.
* `linalg.h <http://github.com/sgorsten/linalg/>`_ -- linear algebra library.
  License: Public Domain.
* `zlib <https://github.com/madler/zlib>`_ -- a subset of the zlib library
  for decompressing gz files, used as a fallback when the zlib library
  is not found in the system. License: zlib.

Not distributed with Gemmi:

* `nanobind <https://github.com/wjakob/nanobind>`_ -- used for creating
  Python bindings. License: 3-clause BSD.
* `zlib-ng <https://github.com/zlib-ng/zlib-ng>`_ -- optional, can be used
  instead of zlib for faster reading of gzipped files.
* `cctbx <https://github.com/cctbx/cctbx_project>`_ -- used in tests
  (if cctbx is not present, these tests are skipped) and
  in scripts that generated space group data and 2-fold twinning operations.
  License: 3-clause BSD.

Mentions:

* `NLOpt <https://github.com/stevengj/nlopt>`_
  was used to try out various optimization methods for class Scaling.
  License: MIT.

Email me if I forgot about something.

List of C++ headers
-------------------

Here is a list of C++ headers in `gemmi/include/`.
This list also provides an overview of the library.

.. include:: headers.rst
