.. meta::
   :google-site-verification: LsEfb1rjo2RL8WOSZGigV11Kgyhtk9v1Vb-6GZFnHKo

GEMMI - GEneral MacroMolecular I/O
==================================

Gemmi is a library developed primarily for use in
**macromolecular crystallography** (MX) programs.

Parts of the library may also be useful in structural bioinformatics
(for symmetry-aware analysis of protein models),
and in other molecular-structure sciences that use CIF files
(as we have the fastest open-source CIF parser).

Gemmi is open-source (MPL_) and portable (Linux, Windows, MacOS).
It is written in C++11, with Python (2 and 3) bindings,
as well as with (soon-to-be-added) C and Fortran 2003 interface and
(possibly in the future) with a subset of functionality translated to
JavaScript for web applications.

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

.. important::

    As of 2018 the library is under intensive development and not
    everything is documented yet. Just ask questions.

Contents
--------

.. toctree::
   :maxdepth: 2

   install
   cif-parser
   sym
   mol
   grid
   hkl
   utils
