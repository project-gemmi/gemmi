
Macromolecular Models
#####################

Reading mmCIF as a generic CIF file gives access
to all the data as it is stored in in the file --
mostly tables with named columns.

When working with macromolecular models a different representation is useful:
a hierarchy of objects.
To a first approximation all macromolecular libraries present the same
hierarchy: model - chain - residue - atom.

While *chain* and *residue* are not good names when referring to
ligands and waters, we use this nomenclature as it is the most popular one.
Some libraries (clipper) call it polymer - monomer - atom.
PDBx/mmCIF uses more general (but not so obvious) terms:
*entity* and *struct_asym* (structural component in asymetric unit)
instead of chain,
and *chem_comp* (chemical component) for residue/monomer.

The mmCIF format introduced two sets of names/numbers:
*label* and *auth* (for "author").
Let us first consider ``_atom_site.label_asym_id`` and
``_atom_site.auth_asym_id``, which are two sets of chain IDs
in the atom table.
It is debatable which one should be used. While the latter is only
an alternative ID according to the mmCIF spec, it is the one that is
used in the PDB format so it may be preferred for consistency.
The usual difference is that the ligands or waters from the author's
"chain" are split into separate (label) "chains", but sometimes
the names of polypeptide chains also change (e.g. C,A,B <-> A,B,C).
Also frequently different are the sequence identifiers
``auth_`` and ``label_seq_id``, analogically to the chain ID.
The two other couples, ``atom_id`` and ``comp_id``, rarely differ,
as we checked in the :ref:`auth_label_example` example.

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


Model - Chain - Residue - Atom
==============================

TODO

Sequence
========

TODO

Chemical Component
==================

TODO

