.. _api:

C++ API Reference
#################

This page documents the Gemmi C++ library API, generated from Doxygen comments
in ``include/gemmi/*.hpp``. It is updated with each pull request in the
API documentation series.

For the Python API, see the `Python API reference <https://project-gemmi.github.io/python-api/>`_.

.. note::

   Documentation coverage is being added incrementally. Headers not yet
   listed here will appear in subsequent pull requests.

Core Data Structures
--------------------

*(Full documentation added in PR 2.)*

.. doxygenfile:: model.hpp
   :project: gemmi

Map and Grid Data
-----------------

*(Stub — full documentation added in PR 6.)*

.. doxygenfile:: grid.hpp
   :project: gemmi

Calculations and Analysis
-------------------------

Geometric calculations, sequence alignment, structure superposition,
neighbour search, contact detection, biological assembly, atom selection,
structure modification, polymer heuristics, and secondary structure assignment.

.. doxygenfile:: calculate.hpp
   :project: gemmi

.. doxygenfile:: align.hpp
   :project: gemmi

.. doxygenfile:: neighbor.hpp
   :project: gemmi

.. doxygenfile:: contact.hpp
   :project: gemmi

.. doxygenfile:: assembly.hpp
   :project: gemmi

.. doxygenfile:: select.hpp
   :project: gemmi

.. doxygenfile:: modify.hpp
   :project: gemmi

.. doxygenfile:: polyheur.hpp
   :project: gemmi

.. doxygenfile:: dssp.hpp
   :project: gemmi

Structure Factor Calculations
-----------------------------

Direct structure factor summation, amplitude normalisation (F→E),
and anisotropic scaling with optional bulk-solvent correction.

.. doxygenfile:: sfcalc.hpp
   :project: gemmi

.. doxygenfile:: ecalc.hpp
   :project: gemmi

.. doxygenfile:: scaling.hpp
   :project: gemmi
