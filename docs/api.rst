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

Structure I/O
-------------

*(Full documentation added in PR 4.)*

.. doxygenfile:: mmcif.hpp
   :project: gemmi

.. doxygenfile:: mmread.hpp
   :project: gemmi

.. doxygenfile:: mmread_gz.hpp
   :project: gemmi

.. doxygenfile:: pdb.hpp
   :project: gemmi

.. doxygenfile:: to_mmcif.hpp
   :project: gemmi

.. doxygenfile:: to_pdb.hpp
   :project: gemmi

.. doxygenfile:: crd.hpp
   :project: gemmi

.. doxygenfile:: smcif.hpp
   :project: gemmi

.. doxygenfile:: mmdb.hpp
   :project: gemmi

.. doxygenfile:: pirfasta.hpp
Map and Grid Data
-----------------

*(Stub — full documentation added in PR 6.)*

.. doxygenfile:: grid.hpp
   :project: gemmi

I/O and Filesystem Utilities
------------------------------

File and directory traversal, gzip support, stream abstractions, PDB path
utilities, and general-purpose string and container helpers.

*(Full documentation added in PR 10.)*

.. doxygenfile:: dirwalk.hpp
   :project: gemmi

.. doxygenfile:: fileutil.hpp
   :project: gemmi

.. doxygenfile:: fstream.hpp
   :project: gemmi

.. doxygenfile:: gz.hpp
   :project: gemmi

.. doxygenfile:: input.hpp
   :project: gemmi

.. doxygenfile:: glob.hpp
   :project: gemmi

.. doxygenfile:: logger.hpp
   :project: gemmi

.. doxygenfile:: pdb_id.hpp
   :project: gemmi

.. doxygenfile:: util.hpp
   :project: gemmi

Low-level Primitives
--------------------

Span and range views, custom iterators, error utilities, fast numeric parsing,
and version information.

*(Full documentation added in PR 10.)*

.. doxygenfile:: span.hpp
   :project: gemmi

.. doxygenfile:: iterator.hpp
   :project: gemmi

.. doxygenfile:: fail.hpp
   :project: gemmi

.. doxygenfile:: atof.hpp
   :project: gemmi

.. doxygenfile:: atox.hpp
   :project: gemmi

.. doxygenfile:: version.hpp
   :project: gemmi

Miscellaneous
-------------

Anomalous scattering addends, bond index, DSN6/BRIX map format, enum/string
conversions, string formatting, statistics, and PyMOL selection language.

*(Full documentation added in PR 10.)*

.. doxygenfile:: addends.hpp
   :project: gemmi

.. doxygenfile:: bond_idx.hpp
   :project: gemmi

.. doxygenfile:: dsn6.hpp
   :project: gemmi

.. doxygenfile:: enumstr.hpp
   :project: gemmi

.. doxygenfile:: sprintf.hpp
   :project: gemmi

.. doxygenfile:: stats.hpp
   :project: gemmi

.. doxygenfile:: pymol_select.hpp
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
