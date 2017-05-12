
Macromolecular Models
#####################

General consideration
=====================

Naming
------

The most useful representation for working with macromolecular models
is a hierarchy of objects.
To a first approximation all macromolecular libraries present the same
hierarchy: model - chain - residue - atom.

While *chain* and *residue* are not good names when referring to
ligands and waters, we use this nomenclature as it is the most popular one.
Some libraries (clipper) call it polymer - monomer - atom.
PDBx/mmCIF uses more general (but not so obvious) terms:
*entity* and *struct_asym* (structural component in asymetric unit)
instead of chain,
and *chem_comp* (chemical component) for residue/monomer.


Disorder (altloc)
-----------------

Apart from the naming, the biggest difference between libraries is
how the disorder is presented. The main options are:

* group together atoms from the same conformer (e.g. cctbx.iotbx)

* group together alternative locations of the same atom (e.g. BioPython)

* leave it to the user (e.g. mmdb and clipper).

Handling alternative conformations may add a lot complexity.
The `iotbx.pdb <https://cci.lbl.gov/cctbx_docs/iotbx/iotbx.pdb.html>`_
documentation says that
"about 90% of the development time invested into iotbx.pdb was in some form
related to alternative conformations".

to be continued...


PDBx/mmCIF format
=================

Reading mmCIF as a generic CIF file gives access
to all the data as it is stored in in the file --
mostly tables with named columns.
While this part of Gemmi provide higher-level API to work with
molecular models, in some situation it is still necessary to understand
the design of the underlying file format.

The main characteristics of the syntax are described in the
:ref:`CIF introduction <cif_intro>`.
Here we focus on things specific to mmCIF/DDL2:

* PDBx/mmCIF dictionary is actually a relational database schema put into DDL.
  Categories correspond to tables. Data items correspond to columns.
  Key data items correspond to primary (or composite) keys in RDBMS.

  While a single block in a single file always describes a single PDB entry,
  the dictionary is designed for any number of entries in one block.
  For example, you always have a single ``_entry.id`` and a single
  ``_struct.title``,
  but if you would have many of them you would need to somehow match them
  (``_entry`` and ``_struct`` are separate tables).
  So you got also ``_struct.entry_id``.
  Should you check if ``_entry.id`` and ``_struct.entry_id`` match each other
  before reading ``_struct.title``? I do not know.

* Any category (RDBMS table) can be written as a CIF loop (table).
  If such a table would have a single row it can be (and always is in wwPDB)
  written as key-value pairs.
  So when accessing a value it is safer to use abstraction that hides this
  syntactic detail (``cif::TableView`` in Gemmi).

* The cost of the flexibility of the CIF format is that it is harder to parse
  than the old PDB format. At least in our opinion -- wwPDB FAQ states
  that the simple context-free grammar of mmCIF is much simpler to parse
  than PDB. Nevertheless, the same FAQ shows how to extract atoms from mmCIF
  with ``grep`` and ``awk`` and assures the reader that all mmCIF files in PDB
  have the same regular layout. So we believe that many in-house (and not
  only) programs are written to parse only mmCIF files in a specific layout.

* Four columns in the atoms (``_atom_site``) table have "author defined
  alternatives" (``.auth_*`` instead of ``.label_*``).
  Two of them, atom name (atom_id) and residue name (comp_id)
  :ref:`almost never <auth_label_example>` differ.
  The other two, chain name (asym_id) and sequence number (seq_id)
  may differ in a confusing way (A,B,C <-> C,A,B).
  Which one is presented to the user depends on a program (the author's
  version seems to be more common). This may lead to funny situations.

* At last, there is a formal distinction between mmCIF and PDBx/mmCIF
  dictionary. The latter is built upon the former. So we have
  the ``pdbx_`` prefix in otherwise random places.

As described above, the mmCIF format has two sets of names/numbers:
*label* and *auth* (for "author").
``atom_id`` and ``comp_id`` almost never differ, so
we ignore the author-defined alternatives.

We cannot ignore the author-defined chain names (``asym_id``) --
they are often used. We split the model into chains based on
the primary mmCIF chain name, but we keep both sets of names.
In all (almost all? to be checked)
PDB entries the author-defined chains are either kept as they are,
or the ligands or waters are moved into separate "chains" (structural units).
The names of polypeptide chains may also be changed.
Typically it happens when originally they were not A, B, C ...

It is debatable which chain name and sequence number should be presented
to the user. While the author-defined ones are "alternative",
they are used in the PDB format so may be preferred for consistency.

We also keep author-defined sequence numbers.
Unlike the primary sequence numbers, they must be used together with
the so-called PDB insertion code.


Model - Chain - Residue - Atom
==============================

TODO

Sequence
========

TODO

Chemical Component
==================

TODO

