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
   :project: gemmi

Map and Grid Data
-----------------

*(Stub — full documentation added in PR 6.)*

.. doxygenfile:: grid.hpp
   :project: gemmi
