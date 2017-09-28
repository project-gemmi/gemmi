
Project GEMMI
=============

Gemmi is a library used primarily in **macromolecular crystallography** (MX)
programs.

Parts of the library may also be useful in structural bioinformatics
(for symmetry-aware analysis of protein models),
and in other molecular-structure sciences that use CIF files
(as we have the fastest open-source CIF parser).

It is an open-source (MPL_), portable (Linux, Windows, MacOS) project,
written in C++11, with Python (2 and 3) bindings,
as well as with (soon-to-be-added) C and Fortran 2003 interface and
(possibly in the future) with a subset of functionality translated to
JavaScript.

.. _MPL: https://www.mozilla.org/en-US/MPL/2.0/

Gemmi is a joint project of
`Global Phasing Ltd <https://www.globalphasing.com/>`_
and `CCP4 <http://www.ccp4.ac.uk>`_,
started in 2017, aiming to:

* enhance macromolecular refinement programs (Refmac and BUSTER),
* replace the CCP4 Coordinate Library (MMDB),
* improve support for PDBx/mmCIF files in the CCP4 and GPhL software suites.

More details will be added as the project progresses.

Source code repository: https://github.com/project-gemmi/gemmi


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

Python 2.7/3.x module
---------------------

To install the gemmi module you need pip, git and not too old
C++ compiler (GCC 4.8+, Clang 3.4+, MSVC 2015+, ICC 16+)::

    pip install git+https://github.com/project-gemmi/gemmi.git

(when the project is more mature and has regular releases, it will be simply
``pip install gemmi``).

Fortran 2003+
-------------

TODO

Utilities
---------

The library comes with growing number of small command-line programs.
When the project is more mature we will provide binaries for Windows, Mac
and Linux. At this moment the utilities are tested only Linux and Mac
and need to be compiled from source::

    git clone https://github.com/project-gemmi/gemmi.git
    cd gemmi/src
    make
