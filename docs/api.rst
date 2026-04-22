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

.. doxygenfile:: model.hpp
   :project: gemmi

.. doxygenfile:: unitcell.hpp
   :project: gemmi

.. doxygenfile:: symmetry.hpp
   :project: gemmi

.. doxygenfile:: metadata.hpp
   :project: gemmi

.. doxygenfile:: elem.hpp
   :project: gemmi

.. doxygenfile:: seqid.hpp
   :project: gemmi

.. doxygenfile:: resinfo.hpp
   :project: gemmi

.. doxygenfile:: small.hpp
   :project: gemmi

.. note::

   The following sections will be populated by subsequent PRs (3–10) in this series.
   See `PR #413 <https://github.com/project-gemmi/gemmi/pull/413>`_ for the full roadmap.
