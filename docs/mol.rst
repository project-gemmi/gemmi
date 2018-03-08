
Molecular models
################

This part of the Gemmi library is for handling structural models of
biomolecules.

The first prerequisite for this is reading and writing popular file
formats. The content of a file is stored as a ``Structure`` object,
with the usual Model-Chain-Residue-Atom hierarchy.
Then we have functions to access and manipulate the structure.
Finally, with time we will grow a toolset of higher-level functionality:
calculation of various metrics, ready-to-use complex operations, etc.

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

Chemical Components
===================

Atoms form residues (monomers) and small molecule components
that the PDB calls chemical components.
The components are described in:

* the Chemical Component Dictionary maintained by the PDB,
* and in so-called CIF files (that often do not conform to the IUCr CIF
  syntax): in the Refmac monomer library (maintained by CCP4)
  and in other libraries.

Gemmi has a built-in mini-database of popular components that is often
sufficient.

**C++**

.. literalinclude:: code/resinfo.cpp
   :language: cpp

**Python**

.. doctest::

    >>> gemmi.find_tabulated_residue('ALA').is_amino()
    True
    >>> gemmi.find_tabulated_residue('DOD').is_water()
    True
    >>> # PDB marks "non-standard" residues as HETATM.
    >>> # Pyrrolysine is now standard - some microbes have it.
    >>> gemmi.find_tabulated_residue('PYL').pdb_standard
    True
    >>> gemmi.find_tabulated_residue('MSE').pdb_standard
    False

To get more complete information we need to first read either the CCD
or a monomer library.

TODO

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

When the unit cell is created by reading a coordinate file
it gets information about symmetry transformations and it can do
a few other things, such as finding symmetry images of a point
or determining if the given coordinates are close to a special position.
This will be covered later on.

Coordinate files
================

Gemmi support the following coordinate file formats:

    * mmCIF (PDBx/mmCIF),
    * PDB (with popular extensions),
    * mmJSON,
    * a binary format (MMTF, binary CIF, or own format) is to be considered.

In this section we first show how to read any of the supported formats,
then we go into details of the individual formats,
and finally we show what can be done with the structural model.

Any format
----------

C++
~~~

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

    enum class CoorFormat { Unknown, Pdb, Cif, Json };

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

TODO: generic write functions

Python
~~~~~~

TODO

``gemmi.Structure`` is documented :ref:`further on <mcra>`.


PDBx/mmCIF format
-----------------

Although the model-chain-residue-atom API abstracts details of
the underlying files, sometimes one needs to know the quirks of the file
format that is being used.

The main characteristics of the CIF syntax are described in the
:ref:`CIF introduction <cif_intro>`.
Here we focus on things specific to mmCIF/DDL2:

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
Gemmi ignores author-defined alternatives.

On the other hand, chain names (``asym_id``) and sequence numbers often
differ and usually the author-defined names should be presented to the user,
as they are the ones used in the PDB format.

In Gemmi, we split the model into chains based on the primary mmCIF
chain name, but we keep both sets of names.
Apart from chain renaming (when the original naming was not A, B, C ...),
it is common that ligands and waters are moved into separate *label* "chains"
(structural units).

Note that unlike the primary sequence numbers,
*author* sequence numbers must be used together with the so-called
PDB insertion code.


C++
~~~

Reading mmCIF files is done in two stages:
file → ``cif::Document`` → ``Structure``.
The same with writing: first a ``cif::Document`` is created or updated
using the data from ``Structure``, and then it is written to disk.

::

    #include <gemmi/cif.hpp>       // file -> cif::Document
    #include <gemmi/gz.hpp>        // uncompressing on the fly
    #include <gemmi/mmcif.hpp>     // read_atoms(cif::Document&) -> Structure
    #include <gemmi/to_mmcif.hpp>  // update_cif_block()
    #include <gemmi/to_cif.hpp>    // cif::Document -> file

    // TODO

``cif::Document`` can be used to access meta-data,
such as the details of the experiment or software used for data processing.
The examples are provided in the :ref:`CIF parser <cif_examples>` section.

Python
~~~~~~

.. code-block:: python

    import gemmi

TODO


mmJSON format
-------------

The mmJSON_ format is a JSON representation of the mmCIF data.
It is available from PDBj:

.. code-block:: none

    curl -o 5MOO.json.gz 'https://pdbj.org/rest/downloadPDBfile?id=5MOO&format=mmjson-all'

Gemmi can read and write files in this format into ``cif::Document``,
which then can be used to prepare ``Structure`` as described
in the mmCIF section.

.. _mmJSON: https://pdbj.org/help/mmjson?lang=en

C++
~~~

::

    #include <gemmi/json.hpp>     // to read
    #include <gemmi/gz.hpp>       // to uncompress on the fly
    #include <gemmi/to_json.hpp>  // to write

    namespace cif = gemmi::cif;

    cif::Document a = cif::read_mmjson_file(path_a);
    cif::Document b = cif::read_mmjson(gemmi::MaybeGzipped(path_b));
    // ... the next steps are the same as for mmCIF

    // TODO: writing

Python
~~~~~~

.. code-block:: python

    import gemmi

TODO



PDB format
----------

Gemmi parses a subset of records from the official
`format specification`__.

__ https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html


Additionally, it supports popular extension that are not included
in this spec:

* two-character chain IDs (columns 21 and 22)
* segment ID (columns 73-76)
* hybrid-36_ encoding of sequence IDs for sequences longer than 9999
  (although we are yet to find an examples for this)
* hybrid-36_ encoding of serial numbers for 99,999+ atoms.

.. _hybrid-36: http://cci.lbl.gov/hybrid_36/

C++
~~~

::

    #include <gemmi/pdb.hpp>     // to read
    #include <gemmi/gz.hpp>      // to uncompress on the fly
    #include <gemmi/to_pdb.hpp>  // to write

TODO

Python
~~~~~~

.. code-block:: python

    import gemmi

TODO


.. _mcra:

Model - Chain - Residue - Atom
==============================

Hierarchy
---------

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
  >>> st[0].residues('A', 12, ' ')
  <gemmi.ResidueGroup [ 12/MET ]>
  >>> met_to_mse(st)
  >>> st[0].residues('A', 12, ' ')
  <gemmi.ResidueGroup [ 12/MSE ]>

Structure
---------

Model
-----

Chain
-----

Disorder (altloc)
~~~~~~~~~~~~~~~~~

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

In Gemmi we expose the *altLoc* field to the user (like mmdb),
and on top of it we offer simple utilities that make working with conformers
easier (similarly to BioPython).


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

Selections
==========

TODO

MMDB selection language.

And perhaps something else (selections used in JMol or PyMOL?)


Sequence
========

TODO

.. _chemcomp:

Chemical Component
==================

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
