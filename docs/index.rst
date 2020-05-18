.. meta::
   :google-site-verification: LsEfb1rjo2RL8WOSZGigV11Kgyhtk9v1Vb-6GZFnHKo

GEMMI - library for structural biology
======================================

Gemmi is a library, accompanied by a set of programs,
developed primarily for use in **macromolecular crystallography** (MX).
For working with:

* macromolecular models (content of PDB, PDBx/mmCIF and mmJSON files),
* refinement restraints (CIF files),
* reflection data (MTZ and mmCIF formats),
* data on a 3D grid (electron density maps, masks, MRC/CCP4 format)
* crystallographic symmetry.

Parts of this library can be useful in structural bioinformatics
(for symmetry-aware analysis of protein models),
and in other molecular-structure sciences that use CIF files
(we have the `fastest <https://github.com/project-gemmi/mmcif-benchmark>`_
open-source CIF parser).

Gemmi is open-source (MPL_) and portable -- it runs on Linux, Windows,
MacOS and even inside a web browser if compiled to WebAssembly.
It is written in C++11, with Python (2 and 3) bindings,
and with partial C and Fortran 2003 interface.

.. _MPL: https://www.mozilla.org/en-US/MPL/2.0/

Occasionally, the project gets sidetracked into
`visualization of the PDB data <https://project-gemmi.github.io/pdb-stats/>`_.

Gemmi is a joint project of
`Global Phasing Ltd <https://www.globalphasing.com/>`_
and `CCP4 <http://www.ccp4.ac.uk>`_.

The project is named after
`Gemmi Pass <https://goo.gl/maps/akBLbGfrGE9j1oWC7>`_.
The name can also be expanded as *GEneral MacroMolecular I/o*.

Source code repository: https://github.com/project-gemmi/gemmi

.. important::

    As of 2020 the library is under intensive development and not
    everything is documented yet. Just ask questions.

Contents
--------

.. toctree::
   :maxdepth: 2

   install
   cif
   symmetry
   mol
   analysis
   grid
   hkl
   utils
   Python API <https://project-gemmi.github.io/gemmi/api/python>