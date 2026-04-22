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

CIF Data Reading and Writing
-----------------------------

*(Full documentation added in PR 3.)*

.. doxygenfile:: cifdoc.hpp
   :project: gemmi

.. doxygenfile:: cif.hpp
   :project: gemmi

.. doxygenfile:: read_cif.hpp
   :project: gemmi

.. doxygenfile:: to_cif.hpp
   :project: gemmi

.. doxygenfile:: to_json.hpp
   :project: gemmi

.. doxygenfile:: json.hpp
   :project: gemmi

.. doxygenfile:: numb.hpp
   :project: gemmi

.. doxygenfile:: ddl.hpp
   :project: gemmi

Map and Grid Data
-----------------

*(Stub — full documentation added in PR 6.)*

.. doxygenfile:: grid.hpp
   :project: gemmi
