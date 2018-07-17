
Installation
============

.. highlight:: none

C++ library
-----------

It is a header-only library. You need to ensure that
the ``include`` and ``third_party`` directories are in your include path
when compiling your program. For example::

    git clone https://github.com/project-gemmi/gemmi.git
    c++ -std=c++11 -Igemmi/include -Igemmi/third_party -O2 my_program.cpp

If you want Gemmi to uncompress gzipped files on the fly
(i.e. if you ``#include <gemmi/gz.hpp>``)
you will also need to link your program with the zlib library.

.. _install_py:

Python 2.7/3.x module
---------------------

To install the gemmi module you need pip, git and not too old
C++ compiler (GCC 4.8+, Clang 3.4+, MSVC 2015+, ICC 17+)::

    pip install gemmi

Alternatively, to install the latest version directly from the repository::

    pip install git+https://github.com/project-gemmi/gemmi.git

or clone the `project <https://github.com/project-gemmi/gemmi/>`_
(or download a zip file) and from the top-level directory do::

    pip install .

On Windows Python 3.5+ should automatically find an appropriate compiler
(MSVC 2015+) . If the compiler is not installed, pip shows a message
with a download link.
For Python 2.7 pip prefers MSVC 2008, which is too old to compile gemmi.
You may still use MSVC 2015 or 2017, but before invoking pip you need to
set the compiler environment with one of these commands::

    "%VS140COMNTOOLS%\..\..\VC\vcvarsall.bat" x64
    "%VS140COMNTOOLS%\..\..\VC\vcvarsall.bat"

If you'd like to use PyPy instead of CPython
we support PyPy2.7 >= 5.7 (although we test it only occasionally).

Fortran and C bindings
----------------------

The Fortran bindings are in early stage and are not documented yet.
They use the ISO_C_BINDING module introduced in Fortran 2003.
You may see the ``fortran/`` directory to know what to expect.
The bindings and usage examples can be compiled with CMake::

    cmake -D USE_FORTRAN=1 .
    make

The C bindings are used only for making Fortan bindings,
but they should be usable on their own.
If you use cmake to build the project
you get a static library ``libcgemmi.a`` that can be used from C,
together with the :file:`fortran/*.h` headers.

Gemmi program
-------------

The library comes with a command-line program. To build it from source::

    git clone https://github.com/project-gemmi/gemmi.git
    cd gemmi
    cmake .
    make
