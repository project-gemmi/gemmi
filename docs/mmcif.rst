
Macromolecular Models
=====================

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
and *chem_comp* (chemical component) as a generalization of the "residue".

Apart from the naming, the biggest differences between libraries is
in presenting the disorder. The main options are:

* group atoms from the same conformer together (e.g. cctbx.iotbx)

* group alternatives together (e.g. BioPython)

* leave it to the user (e.g. mmdb and clipper).

Handling alternative conformations may add a lot complexity.
The `iotbx.pdb <https://cci.lbl.gov/cctbx_docs/iotbx/iotbx.pdb.html>`_
documentation says that
"about 90% of the development time invested into iotbx.pdb was in some form
related to alternative conformations".

to be continued...
