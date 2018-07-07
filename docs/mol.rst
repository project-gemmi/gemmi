
Molecular models
################

This part of the Gemmi library is for handling structural models of
biomolecules.

A model from a file (PDB, mmCIF, etc.) is stored in a ``Structure`` object,
with the usual Model-Chain-Residue-Atom hierarchy.
Apart from low-level functions to access and manipulate the structure,
we keep on expanding a toolset of higher-level functionality:
calculation of popular metrics, ready-to-use complex operations, etc.

Comparing with the usual tools rooted in bioinformatics:

* Gemmi focuses more on working with incomplete models
  (on all stages before they are published and submitted to the PDB),
* and Gemmi is aware of the neighbouring molecules that are implied by
  the crystallographic and non-crystallographic symmetry.


Elements
========

When working with molecular structures it is good to have a basic data
from the periodic table at hand.

**C++**

.. literalinclude:: code/elem.cpp
   :language: cpp


**Python**

.. doctest::

    >>> import gemmi
    >>> gemmi.Element('Mg').weight
    24.305
    >>> gemmi.Element(118).name
    'Og'
    >>> gemmi.Element('Mo').atomic_number
    42

.. _chemcomp:

Chemical Components
===================

Residues (monomers) and small molecule components are called in the PDB
*chemical components*.
Gemmi can use three sources of knowledge about the chemical components:

* built-in basic data about a couple hundreds of the most popular components,
* the Chemical Component Dictionary maintained by the PDB (25,000+ components),
* so-called CIF files from the Refmac monomer library (maintained by CCP4)
  or from other libraries in the same format.

The built-in data is accessed through the function ``find_tabulated_residue``:

**C++**

.. literalinclude:: code/resinfo.cpp
   :language: cpp

**Python**

.. doctest::

    >>> gln = gemmi.find_tabulated_residue('GLN')
    >>> gln.is_amino_acid()
    True
    >>> gln.one_letter_code
    'Q'
    >>> gemmi.find_tabulated_residue('DOD').is_water()
    True
    >>> # PDB marks "non-standard" residues as HETATM.
    >>> # Pyrrolysine is now standard - some microbes have it.
    >>> gemmi.find_tabulated_residue('PYL').is_standard()
    True
    >>> gemmi.find_tabulated_residue('MSE').is_standard()
    False

To get more complete information we need to first read either the CCD
or a monomer library.

TODO

Transformation matrices
=======================

We also need a tiny bit of linear algebra to work with 3D transformations,
such as crystallographic symmetry and NCS operations,
or fractionalization and orthogonalization of coordinates.

A tranformation is represented by class ``Transform`` that has two member
variables: ``mat`` (of type ``Mat33``) and ``vec`` (of type ``Vec3``).
In C++ these types are defined in ``gemmi/math.hpp``.

Coordinates are represented by classes derived from ``Vec3``:

* ``Position`` for coordinates in Angstroms (orthogonal coordinates),
* ``Fractional`` for coordinates relative to the unit cell
  (fractional coordinates).

The only reason to have separate types is to prevent functions that
expect fractional coordinates from accepting orthogonal ones, and vice versa.

For the same reason gemmi also has a ``FTransform``, which is like
``Transform`` but can be applied only to ``Fractional`` coordinates.

**Python**

.. doctest::

    >>> # get NCS transformation from an example pdb file 
    >>> ncs_op = gemmi.read_structure('../tests/1lzh.pdb.gz').ncs[0].tr
    >>> type(ncs_op)
    <class 'gemmi.Transform'>
    >>> ncs_op.mat
    <gemmi.Mat33 [0.97571, -0.2076, 0.06998]
                 [0.2156, 0.96659, -0.13867]
                 [-0.03885, 0.15039, 0.98786]>
    >>> _.determinant()
    1.000003887799667
    >>> ncs_op.vec
    <gemmi.Vec3(-14.1959, 0.72997, -30.5229)>
    >>> # is the 3x3 matrix above orthogonal?
    >>> mat = ncs_op.mat
    >>> identity = gemmi.Mat33()
    >>> mat.multiply(mat.transpose()).approx(identity, epsilon=1e-5)
    True
    >>> ncs_op.apply(gemmi.Vec3(20, 30, 40))
    <gemmi.Vec3(1.8895, 28.4929, 12.7262)>
    >>> ncs_op.inverse().apply(_)
    <gemmi.Vec3(20, 30, 40)>

Unit Cell
=========

When working with a structural models in a crystal, we need to work
with a unit cell, and in particular we need to be able to switch between
orthogonal and fractional coordinates.

**C++**

.. literalinclude:: code/cell.cpp
   :language: cpp

**Python**

.. doctest::

    >>> cell = gemmi.UnitCell(25.14, 39.50, 45.07, 90, 90, 90)
    >>> cell
    <gemmi.UnitCell(25.14, 39.5, 45.07, 90, 90, 90)>
    >>> cell.a, cell.b, cell.c
    (25.14, 39.5, 45.07)
    >>> cell.alpha, cell.beta, cell.gamma
    (90.0, 90.0, 90.0)
    >>> cell.volume
    44755.8621
    >>> cell.fractionalize(gemmi.Position(10, 10, 10))
    <gemmi.Fractional(0.397772, 0.253165, 0.221877)>
    >>> cell.orthogonalize(gemmi.Fractional(0.5, 0.5, 0.5))
    <gemmi.Position(12.57, 19.75, 22.535)>

The UnitCell object can also store a list of symmetry transformations.
This list is populated automatically when reading a coordinate file.
It contains crystallographic symmetry operations. In rare cases
when the file defines strict NCS operarations that are not "given"
(MTRIX record in the PDB format or _struct_ncs_oper in mmCIF)
the list contains also the NCS operations.
With this list we can use:

* ``UnitCell::volume_per_image() -> double`` - returns ``UnitCell::volume``
  divided by the number of the molecule images in the unit cell,
* ``UnitCell::find_nearest_image()`` - TODO document it
* ``UnitCell::is_special_position(const Position& pos, double max_dist=0.8) -> int`` -
  returns the number of nearby symmetry mates of an atom.
  Non-zero only for atoms on special positions.
  For example, returns 3 for an atom on 4-fold symmetry axis.

Reading coordinate files
========================

Gemmi support the following coordinate file formats:

    * mmCIF (PDBx/mmCIF),
    * PDB (with popular extensions),
    * mmJSON,
    * a binary format (MMTF, binary CIF, or own format) is to be considered.

In this section we first show how to read any of the supported formats,
then in the next sections we go into details of the individual formats,
and finally we show what can be done with the structural model.

C++
---

Any of the macromolecular coordinate files supported by Gemmi can be opened
using::

    #include <gemmi/mmread.hpp>
    // ...
    gemmi::Structure st = gemmi::read_structure_file(path);
    std::cout << "This file has " << st.models.size() << " models.\n";

``gemmi::Structure`` is documented :ref:`further on <mcra>`.

The file format above is determined from the file extension.
Alternatively, the format can be specified as the second argument
of ``read_structure_file``, one of the enumeration values::

    enum class CoorFormat { Unknown, Pdb, Mmcif, Mmjson };

Gemmi also has a templated function ``read_structure`` that you can use
to customize how you provide the data (bytes) to the parsers.
This function is used to uncompress gzipped files on the fly:

.. literalinclude:: code/maybegz.cpp
   :language: cpp

If you include the :file:`gz.hpp` header (as in the example above)
the resulting program must be linked with the zlib library.

.. code-block:: console

    $ c++ -std=c++11 -Iinclude -Ithird_party example_above.cpp -lz
    $ ./a.out 2cco.cif.gz
    This file has 20 models.

The :file:`gemmi/mmread.hpp` header includes many other headers
and is relatively slow to compile. If this matters it is recommended
to include this header in only one ``cpp`` file,
and in other files use only :file:`gemmi/model.hpp`.

Python
------

Any of the macromolecular coordinate files supported by Gemmi (gzipped
or not) can be opened using:

.. doctest::
  :hide:

  >>> path =  '../tests/1orc.pdb'


.. doctest::

  >>> gemmi.read_structure(path)  #doctest: +ELLIPSIS
  <gemmi.Structure ...>

``gemmi.Structure`` is documented :ref:`further on <mcra>`.


PDB format
==========

The PDB format evolved between 1970's and 2012. Nowadays the PDB organization
uses PDBx/mmCIF as the primary format and the legacy PDB format is frozen.
Gemmi aims to support files adhering to the official format specification_,
as well as files with commonly present deviations from the standard.
In particular, we support the following extensions:

* two-character chain IDs (columns 21 and 22),
* segment ID (columns 73-76) from PDB v2,
* hybrid-36_ encoding of sequence IDs for sequences longer than 9999
  (although we are yet to find an examples for this),
* hybrid-36_ encoding of serial numbers for more than 99,999 atoms.

The PDB format is well-documented and widely used, nonetheless some of its
features can be easily overlooked. The rest of the section describes
such features. It is for people who are interested in the details
of the PDB format.
You do not need to read it if you just want to use Gemmi and work with
molecular models.

.. note::

   The specification_ aims to describe the format of files generated by the
   wwPDB. It does not aim to specify a format that can be used for data
   exchange between third-party programs.
   Following literally the specification is neither useful nor possible.
   For example: the REVDAT record is mandatory, but using it makes sense
   only for the entries released by the PDB.
   Therefore no software generates files conforming to the specification
   except from the wwPDB software (and even this one is not strictly
   conforming: it writes ``1555`` in the LINK record for the identity operator
   while the specifications requires leaving these fields blank).
   So do not read too much into the specification.

Let us start with the the list of atoms:

.. code-block:: none

   HETATM    8  CE  MSE A   1       8.081   3.884  27.398  1.00 35.65           C
   ATOM      9  N   GLU A   2       2.464   5.718  24.671  1.00 14.40           N
   ATOM     10  CA  GLU A   2       1.798   5.810  23.368  1.00 13.26           C

Standard residues of protein, DNA or RNA are marked as ATOM. Solvent,
ligands, metals, carbohydrates and everything else is marked as HETATM.
Non-standard residues of protein, DNA or RNA are HETATM according to
the wwPDB (so it is not a good idea to remove all HETATM records when one
wants to only remove ligands and solvent),
but many programs and crystallographers prefer to mark them as ATOM.
It is better to not rely on any of the two conventions.

The second field is the serial number of an atom. The wwPDB spec limits
the serial numbers to the range 1--99,999, but the popular extension
called hybrid-36_ allows to have more atoms in the file by using
also letters in this field. If you do not need to interpret the CONECT
records the serial number can be simply ignored.

Columns 13-27 describe the atom's place in the hierarchy.
In the example above they are:

.. code-block:: none

   1      2
   345678901234567

    CE  MSE A   1
    N   GLU A   2
    CA  GLU A   2

Here the CE atom is in chain A, in residue MSE with sequence ID 1.

The atom names (columns 13-16) starts with the element name,
and as a rule columns 13-14 contain only the element name.
Therefore Cα and calcium ion, both named CA, are aligned differently:

.. code-block:: none

   1      2
   345678901234567
    CA  GLU A   2
   CA    CA A 101

This rule has an exception: when the atom name has four characters
it starts in column 13 even if it has a one-letter element code:

.. code-block:: none

   HETATM 6495  CAX R58 A 502      17.143 -29.934   7.180  1.00 58.54           C
   HETATM 6496 CAX3 R58 A 502      16.438 -31.175   6.663  1.00 57.68           C

Columns 23-27 contain a sequence ID. It consists of a number (columns 23-26)
and, optionally, also an insertion code (A-Z) in column 27:

.. code-block:: none

   ATOM  11918  CZ  PHE D 100      -6.852  76.356 -23.289  1.00107.94           C
   ATOM  11919  N   ARG D 100A     -9.676  74.726 -19.958  1.00105.71           N
   ...
   ATOM  11970  CE  MET D 100H     -8.264  83.348 -19.494  1.00107.93           C
   ATOM  11971  N   ASP D 101     -11.329  81.237 -14.804  1.00107.41           N

The insertion codes are the opposite of gaps in the numbering;
both are used to make the numbering consistent with a reference sequence
(and for the same reason the sequence number can be negative).

Another fields that is blank for most of the atoms is altloc.
It is a letter marking an alternative conformation
(columns 17, just before the residue name):

.. code-block:: none

   HETATM  557  O  AHOH A 301      13.464  41.125   8.469  0.50 20.23           O
   HETATM  558  O  BHOH A 301      12.554  42.700   8.853  0.50 26.40           O

Handling alternative conformations adds a lot of complexity,
as it will be described later on in this documentation.
These were all tricky things in the atom list.

One more thing that is sometimes forgotten. In most of the PDB entries
reading the CRYST1 record and the list of atoms is all that is needed
to construct the crystal structure. But to read correctly all PDB files
we also need to read two other records:

* MTRIX -- if marked as not-given it defines operations needed to reconstruct
  the asymmetric unit,
* SCALE -- provides fractionalization matrix. The format of this entry
  is unfortunate: for large unit cells the relative precision of numbers is
  too small. So if coordinates are given in standard settings it is better
  to calculate the fractionalization matrix from the unit cell dimensions.
  But the SCALE record needs to be checked to see if the settings are
  the standard ones.

.. _specification: https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
.. _hybrid-36: http://cci.lbl.gov/hybrid_36/

C++
---

As described in the previous section, all coordinate files can be read
using the same function calls. Additionally, in C++, you may read a selected
file format to avoid linking with the code you do not use:

::

    #include <gemmi/pdb.hpp>     // to read
    #include <gemmi/gz.hpp>      // to uncompress on the fly

    gemmi::Structure st1 = gemmi::read_pdb_file(path);
    // or
    gemmi::Structure st2 = gemmi::read_pdb_file(gemmi::MaybeGzipped(path));


The content of the file can also be read from a string or from a memory::

    Structure read_pdb_from_memory(const char* data, size_t size, const std::string& name);
    Structure read_pdb_string(const std::string& str, const std::string& name);

Python
------

.. code-block:: python

    import gemmi

    # just use interface common for all file formats
    structure = gemmi.read_structure(path)

    # if you have the content of the PDB file in a string:
    structure = gemmi.read_pdb_string(string)


PDBx/mmCIF format
=================

The mmCIF format (more formally: PDBx/mmCIF) became the primary format
used by the wwPDB. The format uses the CIF 1.1 syntax and a semantics
as described by the PDBx/mmCIF DDL2 dictionary.

While this section may clarify a few things, you do not need to read it
to work with mmCIF files.

The main characteristics of the CIF syntax are described in the
:ref:`CIF introduction <cif_intro>`.
Here we focus on things specific to mmCIF:

* PDBx/mmCIF dictionary is clearly based on a relational database schema.
  Categories correspond to tables. Data items correspond to columns.
  Key data items correspond to primary (or composite) keys in RDBMS.

  While a single block in a single file always describes a single PDB entry,
  some relations between tables seem to be designed for any number of entries
  in one block.
  For example, although a file has only one ``_entry.id`` and
  ``_struct.title``, the dictionary uses extra item called ``_struct.entry_id``
  to match the title with id.
  Is it a good practice to check ``_struct.entry_id`` before reading
  ``_struct.title``? Probably not, as I have seen files with missing
  ``_struct.entry_id`` but never (yet) with multiple ``_struct.title``.

* Any category (RDBMS table) can be written as a CIF loop (table).
  If such a table would have a single row it can be (and always is in wwPDB)
  written as key-value pairs.
  So when accessing a value it is safer to use abstraction that hides the
  difference between a loop and a key-value pair
  (``cif::Table`` in Gemmi).

* Arguably, the mmCIF format is harder to parse than the old PDB format.
  Using ``grep`` and ``awk`` to extract atoms will work only with files
  written in a specific layout, usually by a particular software.
  It is unfortunate that the wwPDB FAQ encourages it, so one may expect
  portability problems when using mmCIF.

* The atoms (``_atom_site``) table has four "author defined alternatives"
  (``.auth_*``) that have similar meaning to the "primary" identifiers
  (``.label_*``).
  Two of them, atom name (``atom_id``) and residue name (``comp_id``)
  :ref:`almost never <auth_label_example>` differ.
  The other two, chain name (``asym_id``) and sequence number (``seq_id``)
  may differ in a confusing way (A,B,C <-> C,A,B).
  Which one is presented to the user depends on a program (usually
  the author's version). This may lead to funny situations.

* There is a formal distinction between mmCIF and PDBx/mmCIF
  dictionaries (they are controlled by separate committees).
  The latter is built upon the former. So we have
  the ``pdbx_`` prefix in otherwise random places, to mark tags
  that are not in the vanilla mmCIF.

Here are example lines from a PDB file (3B9F) with the fields
numbered at the bottom:

.. code-block:: none

    ATOM   1033  OE2 GLU H  77      -9.804  19.834 -55.805  1.00 25.54           O
    ATOM   1034  N  AARG H  77A     -4.657  24.646 -55.236  0.11 20.46           N
    ATOM   1035  N  BARG H  77A     -4.641  24.646 -55.195  0.82 22.07           N
     |       |   |  | |  |  | |       |       |       |      |     |             | |
     1       2   3  4 5  6  7 8       9       10      11     12    13           14 15

and the corresponding lines from PDBx/mmCIF v5 (as served by the PDB in 2018):

.. code-block:: none

    ATOM   1032 O OE2 . GLU B 2  72  ? -9.804  19.834  -55.805 1.00 25.54 ? 77   GLU H OE2 1
    ATOM   1033 N N   A ARG B 2  73  A -4.657  24.646  -55.236 0.11 20.46 ? 77   ARG H N   1
    ATOM   1034 N N   B ARG B 2  73  A -4.641  24.646  -55.195 0.82 22.07 ? 77   ARG H N   1
     |       |  | |   |  |  | |   |  |    |       |       |     |    |    |  |    |  | |   |
     1       2 14 N   4  N  N N   N  8    9       10      11    12   13   15 7    5  6 3   N
     |       |  | |   |  label_comp_id    Cartn_x |       |     |    B_iso_or_equiv  | auth_atom_id
     |       id | |   label_alt_id|  pdbx_PDB_ins_code    |     occupancy |  |    |  auth_asym_id
     group_PDB  | label_atom_id   label_seq_id    |       Cartn_z         |  |    auth_comp_id
                type_symbol | label_entity_id     Cartn_y                 |  auth_seq_id   pdbx_PDB_model_num
                            label_asym_id                                 pdbx_formal_charge

``N`` marks columns not present in the PDB file.
Numbers in column 2 differ because in the PDB file the TER record (that mark
end of a polymer) is also assigned a number.

``auth_seq_id`` used to be the full author's sequence ID,
but currently in the wwPDB entries it is only the sequence number;
the insertion code is stored in a separate column (``pdbx_PDB_ins_code``).
Confusingly, this column is placed next to ``label_seq_id`` not ``auth_seq_id``
(``label_seq_id`` is always a positive number and has nothing to do
with the insertion code).

As mentioned above, the mmCIF format has two sets of names/numbers:
*label* and *auth* (for "author").
``label_atom_id`` and ``auth_atom_id`` almost never differ, the same
about ``label_comp_id`` and ``auth_comp_id``.
Gemmi uses author-defined atom and component IDs if they are present,
otherwise it uses *label* ones.

On the other hand, chain names (``asym_id``) and sequence numbers often
differ and in the user interface it is better to use the author-defined
names, for consistency with the PDB format and with the literature.

While this is not guaranteed by the specification, in all PDB entries
each ``auth_asym_id`` "chain" is split into one or more ``label_asym_id``
"chains"; let us call them *subchains*.
Polymer (residues before the TER record in the PDB format) goes into
one subchain, all waters go into another one, and all the other non-polymer
residues are put into separate single-residue subchains.
Non-linear polymers (sugars) are treated as a set of non-polymers.

.. note::

   Having two sets of identifiers in parallel is not a good idea.
   Making them look the same so they can be confused is a terrible design.

   Additionally, the label_* identifiers are not unique: waters have null
   ``label_seq_id`` and therefore all waters in one chain have the same
   identifier. If a water atom is referenced in another table (_struct_conn
   or _struct_site_gen) the label_* identifier is ambiguous,
   so it is necessary to use the auth_* identifier anyway.

This all is quite confusing and lacks a proper documentation.
So once again, now in a color-coded version:

.. raw:: html

 <div class="highlight"><pre style="color:#444">
 ATOM   <b>1032</b> O <span style="color:#d50">OE2 <span style="background-color:#ace">.</span> GLU B</span> 2  <span style="color:#d50">72</span>  <span style="background-color:#ace">?</span> -9.804  19.834  -55.805 1.00 25.54 ? <span style="background-color:#ace">77   GLU H OE2</span> 1
 ATOM   <b>1033</b> N <span style="color:#d50">N   <span style="background-color:#ace">A</span> ARG B</span> 2  <span style="color:#d50">73</span>  <span style="background-color:#ace">A</span> -4.657  24.646  -55.236 0.11 20.46 ? <span style="background-color:#ace">77   ARG H N  </span> 1
 ATOM   <b>1034</b> N <span style="color:#d50">N   <span style="background-color:#ace">B</span> ARG B</span> 2  <span style="color:#d50">73</span>  <span style="background-color:#ace">A</span> -4.641  24.646  -55.195 0.82 22.07 ? <span style="background-color:#ace">77   ARG H N  </span> 1
 </pre></div>

and a couple lines for another file (6any):

.. raw:: html

 <div class="highlight"><pre style="color:#444">
 ATOM   <b>1   </b> N <span style="color:#d50">N   <span style="background-color:#ace">.</span> PHE A</span> 1 <span style="color:#d50">1  </span> <span style="background-color:#ace">?</span> 21.855 30.874 0.439  1.00 29.16 ? <span style="background-color:#ace">17  PHE A N  </span> 1 
 ATOM   <b>2   </b> C <span style="color:#d50">CA  <span style="background-color:#ace">.</span> PHE A</span> 1 <span style="color:#d50">1  </span> <span style="background-color:#ace">?</span> 20.634 31.728 0.668  1.00 26.60 ? <span style="background-color:#ace">17  PHE A CA </span> 1

 ATOM   <b>1630</b> C <span style="color:#d50">CD2 <span style="background-color:#ace">.</span> LEU A</span> 1 <span style="color:#d50">206</span> <span style="background-color:#ace">?</span> 23.900 18.559 1.006  1.00 16.97 ? <span style="background-color:#ace">222 LEU A CD2</span> 1 
 HETATM <b>1631</b> C <span style="color:#d50">C1  <span style="background-color:#ace">.</span> NAG B</span> 2 <span style="color:#d50">.  </span> <span style="background-color:#ace">?</span> 5.126  22.623 37.322 1.00 30.00 ? <span style="background-color:#ace">301 NAG A C1 </span> 1 
 HETATM <b>1632</b> C <span style="color:#d50">C2  <span style="background-color:#ace">.</span> NAG B</span> 2 <span style="color:#d50">.  </span> <span style="background-color:#ace">?</span> 5.434  21.608 38.417 1.00 30.00 ? <span style="background-color:#ace">301 NAG A C2 </span> 1

 HETATM <b>1709</b> O <span style="color:#d50">O   <span style="background-color:#ace">.</span> HOH I</span> 6 <span style="color:#d50">.  </span> <span style="background-color:#ace">?</span> -4.171 14.902 2.395  1.00 33.96 ? <span style="background-color:#ace">401 HOH A O  </span> 1 
 HETATM <b>1710</b> O <span style="color:#d50">O   <span style="background-color:#ace">.</span> HOH I</span> 6 <span style="color:#d50">.  </span> <span style="background-color:#ace">?</span> 9.162  43.925 8.545  1.00 21.30 ? <span style="background-color:#ace">402 HOH A O  </span> 1
 </pre></div>

Each atom site has three independent identifiers:

1. The number in bold is a short and simple one (it does not need to
   be a number according to the mmCIF spec).
2. The hierarchical identifier from the PDB format (blue background)
   is what people usually use. Unfortunately, the arbitrary ordering
   of columns makes it harder to interpret.
3. The new mmCIF identifier (orange) is confusingly similar to 2,
   but it cannot uniquely identify water atoms,
   so it cannot be used in every context.

How other tables in the mmCIF file refer to atom sites?
Some use both 2 and 3 (e.g. _struct_conn), some use only 2 (e.g. _struct_site),
and _atom_site_anisotrop uses all 1, 2 and 3.

C++
---

All coordinate files can be read using the same interface.

Internally, reading mmCIF files is done in two stages:
file → ``cif::Document`` → ``Structure``.

::

    #include <gemmi/cif.hpp>       // file -> cif::Document
    #include <gemmi/gz.hpp>        // uncompressing on the fly
    #include <gemmi/mmcif.hpp>     // cif::Document -> Structure

    namespace cif = gemmi::cif;

    cif::Document doc = cif::read(gemmi::MaybeGzipped(mmcif_file));
    gemmi::Structure structure =  gemmi::make_structure(doc);

The same with writing: first a ``cif::Document`` is created
and then it is written to disk::

    #include <gemmi/to_mmcif.hpp>  // Structure -> cif::Document
    #include <gemmi/to_cif.hpp>    // cif::Document -> file

    gemmi::write_cif_to_file(gemmi::make_mmcif_document(structure), "new.cif");

``cif::Document`` can be used to access meta-data,
such as the details of the experiment or software used for data processing.
The examples are provided in the :ref:`CIF parser <cif_examples>` section.

Python
------

.. doctest::
  :hide:

  >>> mmcif_path =  '../tests/5i55.cif'


.. doctest::

  >>> # just use interface common for all file formats
  >>> structure = gemmi.read_structure(mmcif_path)
  >>>
  >>> # but if you use the cif.read() or cif.read_string() for low level
  >>> # access you may make the structure directly from the cif block:
  >>> # read the content of the PDB file in a string:
  >>> cif_block = gemmi.cif.read(mmcif_path)[0]
  >>> structure = gemmi.make_structure_from_block(cif_block)
  >>>
  >>> # similarly, writing can also be done in two stages
  >>> structure.make_mmcif_document().write_file('new.cif')


mmJSON format
=============

The mmJSON_ format is a JSON representation of the mmCIF data.
This format can be easily parsed with any JSON parser (Gemmi uses sajson).
It is a good alternative to PDBML -- easier to parse
and twice smaller (gzipped).

.. _mmJSON: https://pdbj.org/help/mmjson?lang=en

Files in this format are available from PDBj:

.. code-block:: none

    curl -o 5MOO.json.gz 'https://pdbj.org/rest/downloadPDBfile?id=5MOO&format=mmjson-all'

Gemmi can reads mmJSON files into ``cif::Document``,
as it does with mmCIF files.

C++
---

Reading::

    #include <gemmi/json.hpp>     // JSON -> cif::Document
    #include <gemmi/mmcif.hpp>    // cif::Document -> Structure
    #include <gemmi/gz.hpp>       // to uncompress on the fly

    namespace cif = gemmi::cif;

    cif::Document doc = cif::read_mmjson_file(path);
    // or, to handle gzipped files:
    cif::Document doc = cif::read_mmjson(gemmi::MaybeGzipped(path));
    // and then:
    gemmi::Structure structure =  gemmi::make_structure(doc);

Writing::

    #include <gemmi/to_json.hpp>  // to write

    // cif::Document doc = gemmi::make_mmcif_document(structure); 
    gemmi::write_mmjson_to_stream(ostream, doc);

Python
------

.. doctest::
  :hide:

  >>> mmjson_path =  '../tests/1pfe.json'


.. doctest::

  >>> # just use interface common for all file formats
  >>> structure = gemmi.read_structure(mmjson_path)
  >>>
  >>> # but you can do it in two steps if you wish
  >>> cif_block = gemmi.cif.read_mmjson(mmjson_path)[0]
  >>> structure = gemmi.make_structure_from_block(cif_block)
  >>>
  >>> # similarly, Structure -> mmJSON can also be done in two stages
  >>> json_str = structure.make_mmcif_document().as_json(mmjson=True)


.. _mcra:

Model - Chain - Residue - Atom
==============================

Hierarchy
---------

The most useful representation for working with macromolecular models
is a hierarchy of objects.
To a first approximation all macromolecular libraries present the same
hierarchy: model - chain - residue - atom.

Naming
~~~~~~

While *chain* and *residue* are not good names when referring to
ligands and waters, we use this nomenclature as it is the most popular one.
Some libraries (clipper) call it polymer - monomer - atom.
PDBx/mmCIF uses more general (but not so obvious) terms:
*entity* and *struct_asym* (structural component in asymetric unit)
instead of chain,
and *chem_comp* (chemical component) for residue/monomer.

Alternative conformations
~~~~~~~~~~~~~~~~~~~~~~~~~

Apart from the naming, the biggest difference between libraries is
how the disorder is presented. The main options are:

* group together atoms from the same conformer

* group together alternative locations of the same atom
  (cctbx.iotbx has residue-groups and atom-groups)

* leave it to the user (e.g. mmdb and clipper).

Handling alternative conformations may add a lot complexity.
The `iotbx.pdb <https://cci.lbl.gov/cctbx_docs/iotbx/iotbx.pdb.html>`_
documentation says that
"about 90% of the development time invested into iotbx.pdb was in some form
related to alternative conformations".

Gemmi exposes the *altloc* field to the user (like mmdb).
On top of it it offers simple utilities that make working with conformers
easier:

- functions that ignore all but the main conformations (inspired by BioPython),
- and lightweight proxy objects ResidueSpan and AtomSpan that group
  alternative conformers (inspired by iotbx).

Structure
---------

The object of type Structure that we get from reading a PDB or mmCIF file
may contain multiple models. This adds one more level to the hierarchy.
At this point it may be good to show an example code.
Let us mutate all methionine residues (MET) to selenomethionine (MSE).

**C++**

.. literalinclude:: code/mutate.cpp
   :language: cpp

**Python**

.. testcode::

  import gemmi

  def met_to_mse(st: gemmi.Structure) -> None:
      for model in st:
          for chain in model:
              for residue in chain:
                  if residue.name == 'MET':
                      s_atom = residue['SD']
                      residue.name = 'MSE'
                      s_atom.name = 'SE'
                      s_atom.element = gemmi.Element('Se')

.. doctest::
  :hide:

  >>> st = gemmi.read_structure('../tests/1orc.pdb')
  >>> st[0].sole_residue('A', 12, ' ')
  <gemmi.Residue 12(MET) with 8 atoms>
  >>> met_to_mse(st)
  >>> st[0].sole_residue('A', 12, ' ')
  <gemmi.Residue 12(MSE) with 8 atoms>

Model
-----

Chain
-----

.. doctest::

  >>> st = gemmi.read_structure('../tests/1pfe.cif.gz')
  >>> st.remove_ligands_and_waters()
  >>> chain_b = st[0]['B']
  >>> # iteration goes through all residues and atom sites
  >>> [res.name for res in chain_b]
  ['DSN', 'ALA', 'N2C', 'NCY', 'MVA', 'DSN', 'ALA', 'NCY', 'N2C', 'MVA']
  >>> # The two pairs N2C/NCY above are alternative conformations.
  >>> # Sometimes we want to ignore alternative conformations:
  >>> [res.name for res in chain_b.first_conformer()]
  ['DSN', 'ALA', 'N2C', 'MVA', 'DSN', 'ALA', 'NCY', 'MVA']
  >>>
  >>> chain_a = st[0]['A']
  >>> # the first residue has sequence id '1', but since one sequence id
  >>> # may belong to more than one residue (microheterogeneity)
  >>> # this getter returns ResidueSpan:
  >>> chain_a['1']
  <gemmi.ResidueSpan [1(DG)]>
  >>> first_residue = _[0]
  >>> first_residue
  <gemmi.Residue 1(DG) with 23 atoms>
  >>> sum(atom.altloc == 'B' for atom in first_residue)
  4
  >>> # first_conformer() is for iteration only
  >>> first_residue.first_conformer()  #doctest: +ELLIPSIS
  <gemmi.FirstConformerAtoms object at 0x...>
  >>> len(list(_))  # 23 - 4
  19

TODO

In case of microheterogeneity (two different residues as alternative
conformations at the same position) the sequence ID does not identify
the residue, we also need to supply

Residue
-------

Atom
----

C++
---

Everything will be documented later on.

Here is a minimal example that shows the hierarchy and properties
of each object.

.. literalinclude:: code/structure.cpp
   :language: cpp


TODO

Python
------

.. code-block:: python

    import gemmi

TODO

Neighbor search
===============

Fixed-radius near neighbor search is usually implemented using
a `cell lists <https://en.wikipedia.org/wiki/Cell_lists>`_ method,
also known as binning, bucketing or cell technique.
It has been used in the context of macromolecular structures
since 1960's (in 1966 it was
`described <https://web.stanford.edu/class/sbio228/public/readings/Molecular_Simulation_I_Lecture4/Levinthal_SCIAM_66_Protein_folding.pdf>`_
under the name cubing).

In Gemmi it is implemented in a class named ``SubCells``.
The implementation works with both crystal and non-crystal system and:

* handles crystallographic symmetry (including non-standard settings with
  origin shift that are present in a couple hundreds of PDB entries),
* handles strict NCS (MTRIX record in the PDB format that is not "given";
  in mmCIF it is the _struct_ncs_oper category),
* can find neighbors any number of unit cells apart; surprisingly,
  molecules from different and not neighboring unit cells can be
  in contact, either because of shape (a single chain can be
  :ref:`longer then four unit cells <long_chain>`) or because of
  non-optimal choice of symmetric images in the model
  (some PDB entries have even links between chains more than
  10 unit cells away which cannot be expressed in the 1555 type of notation).

Note that while an atom can be bonded with its own symmetric image,
it sometimes happens that an atom meant to be on a special position
is slightly off, and its symmetric images represent the same atom
(so we may have four nearby images each with occupancy 0.25).
Such images will be returned by the SubCells class as neighbors
and need to be filtered out by the users.

TODO

Selections
==========

TODO

MMDB selection language.

And perhaps something else (selections used in JMol or PyMOL?)


Sequence
========

TODO

.. _pdb_dir:

Local copy of the PDB archive
=============================

In examples that work with the Protein Data Bank archive
we use a local copy of the archive. Like in BioJava,
we assume that the ``$PDB_DIR`` environment variable
points to a directory that contains ``structures/divided/mmCIF`` -- the same
arrangement as on the
`PDB's FTP <ftp://ftp.wwpdb.org/pub/pdb/data/structures/>`_ server.

.. code-block:: console

    $ cd $PDB_DIR
    $ du -sh structures/*/*  # as of Jun 2017
    34G    structures/divided/mmCIF
    25G    structures/divided/pdb
    101G   structures/divided/structure_factors
    2.6G   structures/obsolete/mmCIF

A traditional way to keep an up-to-date local archive is to rsync it
once a week:

.. code-block:: shell

    #!/bin/sh -x
    set -u  # PDB_DIR must be defined
    rsync_subdir() {
      mkdir -p "$PDB_DIR/$1"
      # Using PDBe (UK) here, can be replaced with RCSB (USA) or PDBj (Japan),
      # see https://www.wwpdb.org/download/downloads
      rsync -rlpt -v -z --delete \
	  rsync.ebi.ac.uk::pub/databases/pdb/data/$1/ "$PDB_DIR/$1/"
    }
    rsync_subdir structures/divided/mmCIF
    #rsync_subdir structures/obsolete/mmCIF
    #rsync_subdir structures/divided/pdb
    #rsync_subdir structures/divided/structure_factors

Examples
========

.. _long_chain:

Chain longer than cell
----------------------

Is it possible for a single chain to exceed the size of the unit cell
in one of the directions? How much longer can it be than the cell?

.. literalinclude:: ../examples/long_geom.py
   :language: python
   :lines: 2-

When run on the PDB database (on a local copy of either pdb or mmCIF files)
this script prints too many lines to show here.

.. code-block:: console

  $ ./examples/long_geom.py $PDB_DIR/structures/divided/pdb/
  105M   chain:A   deltaY = 1.225
  208L   chain:A   deltaZ = 1.203
  11BA   chain:A   deltaX = 1.227
  11BA   chain:B   deltaX = 1.202
  ...
  3NWH   chain:A   deltaX = 3.893
  3NWH   chain:B   deltaX = 3.955
  3NWH   chain:C   deltaX = 4.093
  3NWH   chain:D   deltaX = 3.472
  ...
  5XG2   chain:A   deltaX = 4.267
  5XG2   chain:A   deltaZ = 1.467
  ...

As we see, a single chain may be even longer than four unit cells in one
of the directions.  How such chains look like?

For example, here is 3NWH -- a homo-4-mer in P2
(4 x 2 chains per unit cell) -- colored by chain id:

.. image:: img/3nwh.png
    :align: center
    :scale: 100
    :target: https://www.rcsb.org/3d-view/3NWH/

And here is 5XG2 -- a monomer in P21 -- with two copies of the
rainbow-colored chain:

.. image:: img/5xg2.png
    :align: center
    :scale: 100
    :target: https://www.rcsb.org/3d-view/5XG2/


B-factor analysis
-----------------

TODO
