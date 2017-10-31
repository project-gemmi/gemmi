Utilities
#########

The library comes with a few small command-line programs;
running a program is easier than calling a library function.

gemmi-validate
==============

A CIF validator. Apart from checking the syntax it can check most of the rules
imposed by DDL1 and DDL2 dictionaries.

.. literalinclude:: validate-help.txt
   :language: console

gemmi-grep
==========

.. highlight:: console

Searches for a specified tag in CIF files and prints the associated values,
one value per line::

    $ gemmi-grep _refine.ls_R_factor_R_free 5fyi.cif.gz
    5FYI:0.2358
    $ gemmi-grep _refine.ls_R_factor_R_free mmCIF/mo/?moo.cif.gz
    1MOO:0.177
    3MOO:0.21283
    4MOO:0.22371
    5MOO:0.1596
    5MOO:0.1848
    $ gemmi-grep -b _software.name 5fyi.cif.gz
    DIMPLE
    PHENIX

Some of the command-line options are the options of GNU grep (``-l``, ``-H``,
``-n``). As with other utilities, option ``--help`` shows the usage:

.. literalinclude:: grep-help.txt
   :language: console

This is a minimalistic program designed to be used together with Unix
text-processing utilities. For example, it cannot filter values itself,
but one may use grep::

    $ gemmi-grep _pdbx_database_related.db_name /pdb/mmCIF/aa/* | grep EMDB
    4AAS:EMDB
    5AA0:EMDB

Gemmi-grep tries to be simple to use like Unix grep, but at the same time
it is aware of the CIF syntax rules. In particular, ``gemmi-grep _one``
will give the same output for both ``_one 1`` and ``loop_ _one _two 1 2``.
This is helpful in surprising corner cases. For example, when a PDB entry
has two Rfree values (see the 5MOO example above).

If the searched tag is near the beginning of the file, the option ``-O``
makes gemmi-grep much faster. This option tells the program that the file
has only a single block; when the tag is found the program does not need
to parse the rest of the file.

Searching the whole compressed mmCIF archive from the PDB
(35GB of gzipped files) should take on an average computer
between 10 and 30 minutes, depending where the searched tag is located.
This is much faster than with other CIF parsers (to my best knowledge)
and it makes the program useful for ad-hoc PDB statistics::

    $ gemmi-grep -O -b _entity_poly.type /pdb/mmCIF | sort | uniq -c
          1 cyclic-pseudo-peptide
          4 other
          2 peptide nucleic acid
       9905 polydeoxyribonucleotide
        156 polydeoxyribonucleotide/polyribonucleotide hybrid
         57 polypeptide(D)
     168923 polypeptide(L)
       4559 polyribonucleotide
         18 polysaccharide(D)

Going back to moo, we may want to know to what method the Rfree values
correspond::

    $ gemmi-grep _refine.ls_R_factor_R_free -a _refine.pdbx_refine_id mmCIF/mo/?moo.cif.gz
    1MOO:0.177;X-RAY DIFFRACTION
    3MOO:0.21283;X-RAY DIFFRACTION
    4MOO:0.22371;X-RAY DIFFRACTION
    5MOO:0.1596;X-RAY DIFFRACTION
    5MOO:0.1848;NEUTRON DIFFRACTION

Option ``-a`` (``--and``) can be specified many times.
If we would add ``-a _pdbx_database_status.recvd_initial_deposition_date``
we would get the deposition date in each line. In this case it would be
repeated for 5MOO::

    5MOO:0.1596;X-RAY DIFFRACTION;2016-12-14
    5MOO:0.1848;NEUTRON DIFFRACTION;2016-12-14

The option ``-a`` can be used (with some further processing) to generate
relatively sophisticated reports. Here is a little demo:
https://project-gemmi.github.io/pdb-stats/

The major limitation here is that gemmi-grep cannot match
corresponding values from different tables (it is not possible to do this
on the syntax level).
In the example above we have two values from the same table (``_refine``)
and a deposition date (a single value). This works well.
But we are not able to add corresponding wavelengths from ``_diffrn_source``.
If an extra tag (specified with ``-a``) is not in the same table
as the main tag, gemmi-grep uses only the first value for this tag.

Unless we just count the number of value. Counting works for any combination
of tags::

    $ gemmi-grep -c _refln.intensity_meas -a _diffrn_refln.intensity_net r5paysf.ent.gz
    r5paysf:63611;0
    r5payAsf:0;356684

(The file used in this example is structure factor (SF) mmCIF.
Strangely these files in the PDB have extension ``ent`` not ``cif``.)

The first number in the output above is the number of specified intensities.
If you would like to count in also values ``?`` and ``.`` specify
the option ``--raw``::

    $ gemmi-grep --raw -c _refln.intensity_meas r5paysf.ent.gz
    r5paysf:63954
    r5payAsf:0

Gemmi-grep can work with any CIF files but it has one feature
specific to the PDB data. When :ref:`$PDB_DIR <pdb_dir>` is set
one may use PDB codes: just ``5moo`` or ``5MOO`` instead of the path
to ``5moo.cif.gz``. And for convenience, using a PDB code implies
option ``-O``.

The file paths or PDB codes can be read from a file.
For example, if we want to analyse PDB data deposited in 2016
we may first make a file that lists all such files::

    $ gemmi-grep -H -O _pdbx_database_status.recvd_initial_deposition_date $PDB_DIR/structures/divided/mmCIF | \
            grep 2016 >year2016.txt

The 2016.txt file file has lines that start with the filename::

    /hdd/structures/divided/mmCIF/ww/5ww9.cif.gz:5WW9:2016-12-31
    /hdd/structures/divided/mmCIF/ww/5wwc.cif.gz:5WWC:2016-12-31

and a command such as::

    $ gemmi-grep -f year2016.out _diffrn.ambient_temp

will grep only the listed cif files.

Exit status of gemmi-grep has the same meaning as in GNU grep:
0 if a line is selected, 1 if no lines were selected,
and 2 if an error occurred.

gemmi-convert
=============

.. literalinclude:: convert-help.txt
   :language: console

This programs combines a few functions.

CIF -- JSON
-----------

Syntax-level conversion. The JSON representation of the CIF data
can be customized. In particular we support CIF-JSON_ standard from COMCIFS
and mmJSON_ standard from PDBj (the latter is specific to mmCIF files).

The major difference between the two is that CIF-JSON is dictionary-agnostic:
it cannot recognize categories (mmJSON groups by categories),
and it cannot recognize numbers (so it quotes the numbers).
CIF-JSON adds also two extra objects: "CIF-JSON" and "Metadata".
The minor differences are:

 =========== =========== ===========
    CIF        CIF-JSON    mmJSON
 =========== =========== ===========
  data_a      a           data_a
  _tag        _tag        tag
  _CasE       _case       CasE
  .           false       null
  ?           null        null
 =========== =========== ===========


.. _CIF-JSON: http://comcifs.github.io/cif-json
.. _mmJSON: https://pdbj.org/help/mmjson?lang=en

mmCIF -- PDB -- mmJSON
----------------------

Conversion between macromolecular coordinate formats.

We made an effort to format the PDB files we write in the same way as the
software used internally by the PDB, apart from writing fewer records.
Thanks to this, in some scenarios a ``diff`` tool can be used to compare
a PDB file written by Gemmi with an official PDB file from PDB.

The library and the converter also have an option (``--iotbx-compat``)
that formats PDB files similarly to iotbx from cctbx (for example,
in this mode ``TER`` records have no numbers).
We do not aim to be fully compatible with CCTBX, but in many cases
the difference will be only in whitespace.

NCS expansion
-------------
The option ``--expand-ncs`` expands strict NCS,
defined in the ``MTRIX`` record (PDB)
or in the ``_struct_ncs_oper`` table (mmCIF).
By default, new chains have different names than the original ones.
But when used together with ``--iotbx-compat``,
the program mimicks ``iotbx.pdb.expand_ncs`` and leaves the same chain names
while adding distinct segment IDs.


gemmi-map
=========

(work in progress)

.. literalinclude:: map-help.txt
   :language: console


gemmi-mask
==========

(work in progress)

.. literalinclude:: mask-help.txt
   :language: console

gemmi-sg
=========

Prints information about given space group.
