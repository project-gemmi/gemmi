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
   :project: gemmi
