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
