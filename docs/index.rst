.. meta::
   :google-site-verification: LsEfb1rjo2RL8WOSZGigV11Kgyhtk9v1Vb-6GZFnHKo

GEMMI - GEneral MacroMolecular I/O
==================================

Gemmi is a library developed primarily for use in
**macromolecular crystallography** (MX) programs.
For working with:

* macromolecular models (content of PDB, PDBx/mmCIF and mmJSON files),
* refinement restraints (CIF files),
* reflection data (MTZ and mmCIF formats).
* data on a 3D grid (electron density maps, masks, MRC/CCP4 format)
* crystallographic symmetry.

Parts of this library can be useful in structural bioinformatics
(for symmetry-aware analysis of protein models),
and in other molecular-structure sciences that use CIF files
(we have the `fastest <https://github.com/project-gemmi/mmcif-benchmark>`_
open-source CIF parser).

Gemmi is open-source (MPL_) and portable (Linux, Windows, MacOS).
It is written in C++11, with Python (2 and 3) bindings,
and with partial C and Fortran 2003 interface.
A WebAssembly (and perhaps also JavaScript) port is under consideration.

.. _MPL: https://www.mozilla.org/en-US/MPL/2.0/

Gemmi is a joint project of
`Global Phasing Ltd <https://www.globalphasing.com/>`_
and `CCP4 <http://www.ccp4.ac.uk>`_.

Source code repository: https://github.com/project-gemmi/gemmi

.. important::

    As of 2019 the library is under intensive development and not
    everything is documented yet. Just ask questions.

Contents
--------

.. toctree::
   :maxdepth: 2

   install
   cif-parser
   sym
   mol
   analysis
   grid
   hkl
   utils
