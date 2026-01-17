GEMMI API Documentation
========================

This is the API documentation for GEMMI - a library for macromolecular crystallography and structural biology.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

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

.. doxygenfile:: resinfo.hpp
   :project: gemmi

.. doxygenfile:: seqid.hpp
   :project: gemmi

CIF Reading and Writing
------------------------

.. doxygenfile:: cif.hpp
   :project: gemmi

.. doxygenfile:: cifdoc.hpp
   :project: gemmi

.. doxygenfile:: read_cif.hpp
   :project: gemmi

.. doxygenfile:: ddl.hpp
   :project: gemmi

.. doxygenfile:: to_cif.hpp
   :project: gemmi

.. doxygenfile:: to_json.hpp
   :project: gemmi

.. doxygenfile:: numb.hpp
   :project: gemmi

.. doxygenfile:: json.hpp
   :project: gemmi

Structure File I/O
-------------------

.. doxygenfile:: mmcif.hpp
   :project: gemmi

.. doxygenfile:: mmread.hpp
   :project: gemmi

.. doxygenfile:: mmread_gz.hpp
   :project: gemmi

.. doxygenfile:: mmcif_impl.hpp
   :project: gemmi

.. doxygenfile:: pdb.hpp
   :project: gemmi

.. doxygenfile:: to_mmcif.hpp
   :project: gemmi

.. doxygenfile:: to_pdb.hpp
   :project: gemmi

.. doxygenfile:: crd.hpp
   :project: gemmi

.. doxygenfile:: mmdb.hpp
   :project: gemmi

Reflection Data
---------------

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

.. doxygenfile:: merge.hpp
   :project: gemmi

.. doxygenfile:: binner.hpp
   :project: gemmi

.. doxygenfile:: asudata.hpp
   :project: gemmi

.. doxygenfile:: reciproc.hpp
   :project: gemmi

Map and Grid Data
-----------------

.. doxygenfile:: grid.hpp
   :project: gemmi

.. doxygenfile:: ccp4.hpp
   :project: gemmi

.. doxygenfile:: recgrid.hpp
   :project: gemmi

.. doxygenfile:: dencalc.hpp
   :project: gemmi

.. doxygenfile:: asumask.hpp
   :project: gemmi

.. doxygenfile:: solmask.hpp
   :project: gemmi

.. doxygenfile:: floodfill.hpp
   :project: gemmi

Calculations
------------

.. doxygenfile:: calculate.hpp
   :project: gemmi

.. doxygenfile:: neighbor.hpp
   :project: gemmi

.. doxygenfile:: contact.hpp
   :project: gemmi

.. doxygenfile:: assembly.hpp
   :project: gemmi

.. doxygenfile:: fourier.hpp
   :project: gemmi

.. doxygenfile:: sfcalc.hpp
   :project: gemmi

.. doxygenfile:: ecalc.hpp
   :project: gemmi

Structure Factor and Scattering
--------------------------------

.. doxygenfile:: formfact.hpp
   :project: gemmi

.. doxygenfile:: it92.hpp
   :project: gemmi

.. doxygenfile:: c4322.hpp
   :project: gemmi

.. doxygenfile:: fprime.hpp
   :project: gemmi

.. doxygenfile:: neutron92.hpp
   :project: gemmi

Chemistry and Restraints
-------------------------

.. doxygenfile:: chemcomp.hpp
   :project: gemmi

.. doxygenfile:: monlib.hpp
   :project: gemmi

.. doxygenfile:: topo.hpp
   :project: gemmi

.. doxygenfile:: bond_idx.hpp
   :project: gemmi

.. doxygenfile:: riding_h.hpp
   :project: gemmi

.. doxygenfile:: polyheur.hpp
   :project: gemmi

.. doxygenfile:: linkhunt.hpp
   :project: gemmi

.. doxygenfile:: to_chemcomp.hpp
   :project: gemmi

Small Molecules
---------------

.. doxygenfile:: small.hpp
   :project: gemmi

.. doxygenfile:: smcif.hpp
   :project: gemmi

Sequences and Alignment
------------------------

.. doxygenfile:: pirfasta.hpp
   :project: gemmi

.. doxygenfile:: seqalign.hpp
   :project: gemmi

.. doxygenfile:: seqtools.hpp
   :project: gemmi

.. doxygenfile:: align.hpp
   :project: gemmi

Selection and Analysis
----------------------

.. doxygenfile:: select.hpp
   :project: gemmi

.. doxygenfile:: pymol_select.hpp
   :project: gemmi

.. doxygenfile:: dssp.hpp
   :project: gemmi

Modification
------------

.. doxygenfile:: modify.hpp
   :project: gemmi

Refinement
----------

.. doxygenfile:: scaling.hpp
   :project: gemmi

.. doxygenfile:: levmar.hpp
   :project: gemmi

.. doxygenfile:: addends.hpp
   :project: gemmi

Crystallography
---------------

.. doxygenfile:: cellred.hpp
   :project: gemmi

.. doxygenfile:: twin.hpp
   :project: gemmi

Utilities
---------

.. doxygenfile:: util.hpp
   :project: gemmi

.. doxygenfile:: fileutil.hpp
   :project: gemmi

.. doxygenfile:: dirwalk.hpp
   :project: gemmi

.. doxygenfile:: input.hpp
   :project: gemmi

.. doxygenfile:: gz.hpp
   :project: gemmi

.. doxygenfile:: fstream.hpp
   :project: gemmi

.. doxygenfile:: atof.hpp
   :project: gemmi

.. doxygenfile:: atox.hpp
   :project: gemmi

.. doxygenfile:: sprintf.hpp
   :project: gemmi

.. doxygenfile:: span.hpp
   :project: gemmi

.. doxygenfile:: iterator.hpp
   :project: gemmi

.. doxygenfile:: fail.hpp
   :project: gemmi

.. doxygenfile:: version.hpp
   :project: gemmi

.. doxygenfile:: utf.hpp
   :project: gemmi

.. doxygenfile:: logger.hpp
   :project: gemmi

.. doxygenfile:: stats.hpp
   :project: gemmi

.. doxygenfile:: blob.hpp
   :project: gemmi

.. doxygenfile:: pdb_id.hpp
   :project: gemmi

.. doxygenfile:: glob.hpp
   :project: gemmi

.. doxygenfile:: flat.hpp
   :project: gemmi

.. doxygenfile:: enumstr.hpp
   :project: gemmi

Mathematical Functions
----------------------

.. doxygenfile:: math.hpp
   :project: gemmi

.. doxygenfile:: eig3.hpp
   :project: gemmi

.. doxygenfile:: bessel.hpp
   :project: gemmi

.. doxygenfile:: qcp.hpp
   :project: gemmi

Serialization
-------------

.. doxygenfile:: serialize.hpp
   :project: gemmi

Interoperability
----------------

.. doxygenfile:: interop.hpp
   :project: gemmi

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
