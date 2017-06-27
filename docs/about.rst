
About this project
==================

Gemmi is an open-source (MPL_) library aiming to:

.. _MPL: https://www.mozilla.org/en-US/MPL/2.0/

* read and write CIF, mmCIF and PDB files,
* provide other functionality needed by macromolecular refinement programs,
* replace the CCP4 Coordinate Library (MMDB).

It is a joint project of
`Global Phasing Ltd <https://www.globalphasing.com/>`_
and
`CCP4 <http://www.ccp4.ac.uk>`_,
started in 2017 and with funding secured until 2020.

More details will be added as the project progresses.

Source code repository: https://github.com/project-gemmi/gemmi


Installation
============

.. highlight:: none

For C++
-------

At this moment it is a header-only library, so you need to ensure that
the ``include`` and ``third_party`` directories are in your include path
when compiling your program. For example::

    git clone https://github.com/project-gemmi/gemmi.git
    c++ -std=c++11 -I. -Igemmi/third_party -O2 my_program.cpp

If you'd like Gemmi to uncompress gzipped files on the fly,
i.e. when you ``#include <gemmi/cifgz.hpp>`` or
``#include <gemmi/pdbgz.hpp>``,
you will also need to link your program with the zlib library.

For Python 2.7 or 3.x
---------------------

To install the gemmi module you need pip, git and not too old
C++ compiler (GCC 4.8+, Clang 3.4+, MSVC 2015+, ICC 16+)::

    pip install git+https://github.com/project-gemmi/gemmi.git

(when the project is more mature and has regular releases, it will be simply
``pip install gemmi``).

