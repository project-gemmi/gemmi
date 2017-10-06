
Molecular models
################

This part of the Gemmi library is for handling structural models of
biomolecules.

The first prerequisite for this is reading and writing popular file
formats. Then we have functions to access and manipulate
the structure in the usual Model-Chain-Residue-Atom hierarchy.
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
from periodic table at hand.

C++
---

.. literalinclude:: doc_elem.cpp
   :language: cpp


Python
------

.. code-block:: pycon

    >>> from gemmi import mol
    >>> mol.Element('Mg').weight
    24.305
    >>> mol.Element(118).name
    'Og'
    >>> mol.Element('Mo').atomic_number
    42


Coordinate files
================

Gemmi support the following coordinate file formats:

    * mmCIF (PDBx/mmCIF),
    * PDB (with popular extensions),
    * mmJSON,
    * a binary format (MMTF, binary CIF, or own format) is to be considered.

The next sections discuss details of the individual formats
and functions specific to them.
We also have a generic ``read`` and ``write`` functions that handle
all the supported formats, with the format either specified
or inferred from the file extension.

Generic access
--------------

TODO

The data read from a file is kept in a ``Structure``,
which is documented :ref:`further on <mcra>`.

PDBx/mmCIF format
-----------------

While this part of Gemmi provide higher-level API to work with
molecular models, in some situation it is still necessary to understand
the design of the underlying file format.

The main characteristics of the syntax are described in the
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
  (``cif::TableView`` in Gemmi).

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

CIF-based interface

Structure-based interface

combination of both

Structural models from a file are stored in ``Structure``.

Other details from mmCIF, such as literature citations,
can be accessed through a generic CIF interface (``cif::Block``).

::

    #include <gemmi/mmcif.hpp>     // to read
    #include <gemmi/gz.hpp>        // to uncompress on the fly
    #include <gemmi/to_mmcif.hpp>  // to write

TODO

Python
~~~~~~

.. code-block:: python

    from gemmi import mol

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

    from gemmi import mol

TODO


mmJSON format
-------------

The mmJSON_ format is a JSON representation of the mmCIF data.
It is available from PDBj:

.. code-block:: none

    curl -o 5MOO.json.gz 'https://pdbj.org/rest/downloadPDBfile?id=5MOO&format=mmjson-all'

Gemmi can read and write files in this format in a similar way as it reads
and write mmCIF files.

.. _mmJSON: https://pdbj.org/help/mmjson?lang=en

C++
~~~

::

    #include <gemmi/json.hpp>     // to read
    #include <gemmi/gz.hpp>       // to uncompress on the fly
    #include <gemmi/to_json.hpp>  // to write
    namespace mol = gemmi::mol;

    mol::Document = mol::read_mmjson_file(path);
    mol::Document = mol::read_mmjson(MaybeGzipped(path));

TODO

Python
~~~~~~

.. code-block:: python

    from gemmi import mol

TODO

.. _mcra:

Model - Chain - Residue - Atom
==============================

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

C++
---

::

    #include <gemmi/mm.hpp>     // to read/write all the supported formats
    #include <gemmi/gz.hpp>     // to uncompress on the fly
    namespace mol = gemmi::mol;

    mol::Content = mol::read(path);

TODO

Python
------

.. code-block:: python

    from gemmi import mol

TODO

Selections
==========

TODO


Sequence
========

TODO

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

B-factor analysis
-----------------

TODO
