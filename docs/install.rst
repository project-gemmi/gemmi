
Installation
============

.. highlight:: none

C++ library
-----------

Gemmi was a header-only library until ver. 0.6.0.
Some parts, such as CIF syntax parsing, remain header-only.
If you happen to use only these parts,
just ensure that Gemmi's `include` directory is in
your project's include path. For example::

    git clone https://github.com/project-gemmi/gemmi.git
    c++ -Igemmi/include -O2 my_program.cpp

However, in most cases, you need to build the gemmi_cpp library
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

.. _webassembly:

WebAssembly
-----------

The Gemmi library can be compiled with Emscripten to WebAssembly.
Since compiling the entire library is unlikely to be necessary,
we'll show how to compile a subset needed for a particular project,
adding bindings for JavaScript. We present two approaches.

With Embind
~~~~~~~~~~~

The `wasm/` subdirectory contains bindings that use Embind to expose
C++ classes to JavaScript. Currently, they consist of two parts:

* Minimal bindings to the macromolecular `Structure`
  that allow reading a PDB or mmCIF file and iterating over models, chains,
  residues and atoms (see `mol.test.js`). This serves as an example
  and a starting point for further work (which can be carried on either
  as part of gemmi or in the user's own project). Feel free to reach out
  if you have questions.
* Bindings to class `Mtz` that enable map calculation (via FFT)
  from map coefficients. This part was previously provided in the separate
  `mtz module <https://www.npmjs.com/package/mtz>`_,
  the first library to enable the use of MTZ files in molecular graphics apps.

The files from the Gemmi library used for building the wasm module
are listed as `GEMMI_OBJS` in the `Makefile`.

With C API
~~~~~~~~~~

As part of the Gemmi project, we maintain a set of
`web tools <https://project-gemmi.github.io/wasm/>`_ (mostly file converters),
which are single-page applications powered by Gemmi functions in WASM.
The source code of these tools is in the
`wasm repository <https://github.com/project-gemmi/wasm>`_ (not to be confused
with the wasm subdirectory of the gemmi repo -- one of them should be renamed).

These tools don't use the bindings described above. They demonstrate
an alternative approach. For each page we wrote a dedicated C++ function,
with a C API, that performs the bulk of the work. The bindings to these
functions are generated using Emscripten (without Embind).
This approach -- writing part of the web app in C++ -- is more performant,
as it keeps all computations on the WebAssembly side and minimizes
the number of calls across the JS/WASM boundary.

Check the Makefiles in subdirectories to see how the wasm modules are built.


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
