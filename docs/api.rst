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

Scattering, Math, and Geometry
-------------------------------

Form factor tables, anomalous scattering, Bessel function helpers, core
linear-algebra primitives, unit-cell reduction, and numerical tools used
throughout structure-factor and density calculations.

*(Full documentation added in PR 9.)*

.. doxygenfile:: fprime.hpp
   :project: gemmi

.. doxygenfile:: formfact.hpp
   :project: gemmi

.. doxygenfile:: it92.hpp
   :project: gemmi

.. doxygenfile:: c4322.hpp
   :project: gemmi

.. doxygenfile:: neutron92.hpp
   :project: gemmi

.. doxygenfile:: bessel.hpp
   :project: gemmi

.. doxygenfile:: math.hpp
   :project: gemmi

.. doxygenfile:: cellred.hpp
   :project: gemmi

Sequence Alignment and Twinning
--------------------------------

Sequence utilities, pairwise alignment, twinning-law discovery, and
interoperability helpers.

*(Full documentation added in PR 9.)*

.. doxygenfile:: seqtools.hpp
   :project: gemmi

.. doxygenfile:: seqalign.hpp
   :project: gemmi

.. doxygenfile:: twin.hpp
   :project: gemmi

.. doxygenfile:: serialize.hpp
   :project: gemmi

.. doxygenfile:: interop.hpp
   :project: gemmi

.. doxygenfile:: flat.hpp
   :project: gemmi

.. doxygenfile:: smarts.hpp
   :project: gemmi

Density Analysis and Numerical Methods
---------------------------------------

Electron density blob finding, isosurface extraction, Levenberg-Marquardt
least-squares minimization, and quaternion-based superposition (QCP).

*(Full documentation added in PR 9.)*

.. doxygenfile:: blob.hpp
   :project: gemmi

.. doxygenfile:: isosurface.hpp
   :project: gemmi

.. doxygenfile:: levmar.hpp
   :project: gemmi

.. doxygenfile:: qcp.hpp
Reflection Data
---------------

*(Full documentation added in PR 5.)*

.. doxygenfile:: mtz.hpp
   :project: gemmi

.. doxygenfile:: refln.hpp
   :project: gemmi

.. doxygenfile:: cif2mtz.hpp
   :project: gemmi

.. doxygenfile:: mtz2cif.hpp
   :project: gemmi

.. doxygenfile:: xds_ascii.hpp
   :project: gemmi

.. doxygenfile:: xds2mtz.hpp
   :project: gemmi

.. doxygenfile:: intensit.hpp
   :project: gemmi

.. doxygenfile:: binner.hpp
   :project: gemmi

.. doxygenfile:: asudata.hpp
   :project: gemmi

Map and Grid Data
-----------------

*(Stub — full documentation added in PR 6.)*

.. doxygenfile:: grid.hpp
   :project: gemmi
