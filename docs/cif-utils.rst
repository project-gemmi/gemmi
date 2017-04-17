Utilities and Examples
######################

The CIF-handling part of the Gemmi library comes with a couple of utilities,
and several examples that show how you can use the library in your scripts.

Currently, all the examples perform PDB-wide analysis on a local copy
of the mmCIF archive.

Utilities
=========

(the names are tentative and will be changed to gemmi-something,
or perhaps we'll have a single executable ``gemmi`` with subcommands).

validate
--------

A validator that checks the syntax and, optionally, also ontology
using a corresponding DDL1/DDL2 dictionary.
(checking with DDL1 is mostly finished, DDL2 is only started).

It has a few options, run ``validate -h`` for details.

to_json
-------

Converts CIF to JSON. It does not try to preserve all the information
(like the converter included in cod-tools_),
but rather aims for simple output (similar to the converter `in Jmol`_).
It is useful for testing the parser.  Usage::

    to_json my.cif > my.json

.. _cod-tools: https://github.com/sauliusg/cod-tools
.. _in Jmol: https://sourceforge.net/p/jmol/mailman/message/35622017/


.. _cif_examples:

Examples
========

.. highlight:: python

The examples here use Python, as it is the most popular language
for this kind of tasks.
Full working code code can be found in the examples__ directory.

__ https://github.com/project-gemmi/gemmi/tree/master/examples

Amino acid frequency
--------------------

Let say we see in a `paper <https://doi.org/10.1093/nar/gkw978>`_
amino acid frequency averaged over 5000+ proteomes
(Ala 8.76%, Cys 1.38%, Asp 5.49%, etc),
and we want to compare it with the frequency in the PDB database.
So we write a little script

.. literalinclude:: ../examples/aafreq.py

We can run this script on a
`local copy <https://www.wwpdb.org/download/downloads>`_ of the PDB database
in the mmCIF format (30GB+ gzipped, don't uncompress),
and get such an output:

.. code-block:: none

    200L ALA:17 LEU:15 LYS:13 ARG:13 THR:12 ASN:12 GLY:11 ILE:10 ASP:10 VAL:9 GLU:8 SER:6 TYR:6 GLN:5 PHE:5 MET:5 PRO:3 TRP:3 HIS:1
    ...
    4ZZZ LEU:40 LYS:33 SER:30 ASP:25 GLY:25 ILE:25 VAL:23 ALA:19 PRO:18 GLU:18 ASN:17 THR:16 TYR:16 GLN:14 ARG:11 PHE:10 HIS:9 MET:9 CYS:2 TRP:2
    TOTAL LEU:8.90% ALA:7.80% GLY:7.34% VAL:6.95% GLU:6.55% SER:6.29% LYS:6.18% ASP:5.55% THR:5.54% ILE:5.49% ARG:5.35% PRO:4.66% ASN:4.16% PHE:3.88% GLN:3.77% TYR:3.45% HIS:2.63% MET:2.14% CYS:1.38% TRP:1.35%

On my laptop it takes about an hour, using a single core.
Most of this hour is spent on tokenizing the CIF files and copying
the content into a DOM structure, what could be largely avoided given
that we use only sequences not atoms.
But it is not worth to optimize one-off script.
The same goes for using multiple core.

Search PDB by elements
----------------------

Let say we want to be able to search PDB by specifying a set of elements
present in the model. First we write down elements present in each
PDB entry::

    block = cif.read_any(path).sole_block()
    elems = set(block.find_loop("_atom_site.type_symbol"))
    print(name + ' ' + ' '.join(elems))

This example ended up overdone a bit. The code resides in a
`separate repository <https://github.com/project-gemmi/periodic-table>`_.

Demo: `<https://project-gemmi.github.io/periodic-table/>`_


Solvent content vs resolution
-----------------------------

.. |Vm| replace:: *V*\ :sub:`M`
.. |Vs| replace:: *V*\ :sub:`S`
.. |dmin| replace:: *d*\ :sub:`min`

Let say that we come across one of the articles
(`Acta Cryst D <https://www.ncbi.nlm.nih.gov/pubmed/24914969>`_
or `CNN <http://www.phenix-online.org/newsletter/CCN_2015_01.pdf#page=14>`_)
by C. X. Weichenberger and B. Rupp
about solvent content probabilities as a function of |dmin|.
We like the solvent-content-vs-resolution plot there
and would like to make a more detailed version of it.

In macromolecular crystallography solvent content is conventionally
estimated as:

    |Vs| = 1 -- 1.230 / |Vm|

where |Vm| is Matthews coefficient defined as |Vm|\ =\ *V*\ /\ *m*
(volume of the asymmetric unit over the molecular weight of all
molecules in this volume).

Weichenberger and Rupp calculated |Vs| themselves, but we will simply
use the values present in mmCIF files.
So first we extract |Vm|, |Vs| and |dmin|,
as well as deposition date and group ID (which will be explained later).

.. literalinclude:: ../examples/matthews.py
   :pyobject: gather_data

After a few checks (see :file:`examples/matthews.py`) we may decide to use
only those PDB entries in which |Vm| and |Vs| are consistent.
To make things more interesting we plot separately depositions from before
and after 1/1/2015.

.. figure:: img/solvent-content-old.png
   :alt: A plot of solvent content vs resolution. PDB entries up to 2014.
   :align: center
   :scale: 100

   About 90,000 PDB entries (up to 2014)
   plotted as intentionally undersmoothed kernel density estimate.
   The code used to produce this plot is in :file:`examples/matthews.py`.

Ripples in the right subplot show that in many entries |Vs|
is reported as integer, so we should calculate it ourselves
(like Weichenberger & Rupp) for better precision.
Ripples in the top plot show that we should use a less arbitrary metric
of resolution than |dmin| (but it's not so easy).
Or we could just smooth them out by changing parameters of this plot.

.. figure:: img/solvent-content-recent.png
   :alt: A plot of solvent content vs resolution. PDB entries since 2015.
   :align: center
   :scale: 100

   About 17,000 PDB entries (2015 - April 2017)
   plotted as intentionally undersmoothed kernel density estimate.
   The code used to produce this plot is in :file:`examples/matthews.py`.


On the left side of the yellow egg you can see dark stripes
caused by *group depositions*, which were introduced by PDB in 2016.
As of Apr 2017 only a few groups have been deposited (it's just the beginning).
They came from two European high-throughput beamlines and
serve as an illustration of how automated software can analyze hundreds
of similar samples (fragment screening) and submit them quickly to the PDB.

We can easily filter out group depositions -- that is why we stored
group IDs when reading mmCIF files.

A systematic research of the same protein, with dozens PDB submissions,
can also make a spot on our plot.
For example, at |Vs|\ ≈66.5%, |dmin| 2.5-3A we can see the 20S proteasome
studied by Huber *et al*.
But this should not skew the final statistics significantly.

Entity weight
-------------

TODO
