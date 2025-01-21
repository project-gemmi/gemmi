.. meta::
   :google-site-verification: LsEfb1rjo2RL8WOSZGigV11Kgyhtk9v1Vb-6GZFnHKo

Overview
########

What is it for?
===============

Gemmi is a library, accompanied by a :ref:`set of programs <program>`,
developed primarily for use in **structural biology**,
and in particular in **macromolecular crystallography** (MX).
For working with:

* macromolecular models (content of PDB, PDBx/mmCIF and mmJSON files),
* refinement restraints (CIF files) and small molecule models,
* reflection data (MTZ and mmCIF formats),
* crystallographic symmetry,
* data on a 3D grid with crystallographic symmetry
  (electron density maps, masks, MRC/CCP4 format)

Parts of this library can be useful in structural bioinformatics
(for symmetry-aware analysis of protein models),
in chemical crystallography and in other molecular-structure sciences
that use CIF files
(we have the `fastest <https://github.com/project-gemmi/mmcif-benchmark>`_
open-source CIF parser).

Gemmi is open-source (MPL_) and portable -- it runs on Linux, Windows,
MacOS and even inside a web browser if compiled to WebAssembly
(`here <https://project-gemmi.github.io/wasm/>`__ and
`here <https://www.npmjs.com/package/mtz>`__).
It is written in C++17, with Python (3.8+) bindings,
and with partial C and Fortran 2003 interface.

.. _MPL: https://www.mozilla.org/en-US/MPL/2.0/

Occasionally, the project gets sidetracked into
`visualization of the PDB data <https://project-gemmi.github.io/pdb-stats/>`_.

Gemmi is a joint project of
`Global Phasing Ltd <https://www.globalphasing.com/>`_
and `CCP4 <http://www.ccp4.ac.uk>`_.
It is named after
`Gemmi Pass <https://goo.gl/maps/akBLbGfrGE9j1oWC7>`_.
The name can also be expanded as *GEneral MacroMolecular I/o*.

Source code repository: https://github.com/project-gemmi/gemmi

.. note::

    You can ask questions in
    `Discussions <https://github.com/project-gemmi/gemmi/discussions>`_
    or `Issues <https://github.com/project-gemmi/gemmi/issues>`_ on GitHub.
    Alternatively, send me_ an email.

.. _me: wojdyr+gemmi@gmail.com

Contents
========

.. toctree::
   :maxdepth: 1

   Overview <self>
   install
   program

.. toctree::
   :caption: Prerequisites
   :maxdepth: 2

   cif
   symmetry
   cell
   misc

.. toctree::
   :caption: Working with Molecules
   :maxdepth: 2

   chemistry
   mol
   analysis

.. toctree::
   :caption: Working with Data
   :maxdepth: 2

   grid
   hkl
   scattering

.. toctree::
   :caption: Other Docs

   ChangeLog <https://github.com/project-gemmi/gemmi/releases>
   Python API reference <https://project-gemmi.github.io/python-api/>
   C++ API reference <https://project-gemmi.github.io/cxx-api/>

Credits
=======

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

