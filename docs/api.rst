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

Chemistry and Restraints
------------------------

Chemical component definitions, monomer library, topology of restraints
applied to a model, hydrogen placement, link hunting, and related I/O helpers.

*(Full documentation added in PR 8.)*

.. doxygenfile:: chemcomp.hpp
   :project: gemmi

.. doxygenfile:: chemcomp_xyz.hpp
   :project: gemmi

.. doxygenfile:: ener_lib.hpp
   :project: gemmi

.. doxygenfile:: monlib.hpp
   :project: gemmi

.. doxygenfile:: topo.hpp
   :project: gemmi

.. doxygenfile:: riding_h.hpp
   :project: gemmi

.. doxygenfile:: linkhunt.hpp
   :project: gemmi

.. doxygenfile:: to_chemcomp.hpp
   :project: gemmi

.. doxygenfile:: mmcif_impl.hpp
   :project: gemmi
