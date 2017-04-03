What are STAR, CIF, DDL and mmCIF?
==================================

(in case someone comes here when looking for a serialization format)

STAR is a human-readable data serialization format
(think XML or YAML, but also ASN.1)
widely used in molecular-structure sciences.

CIF (Crystallographic Information File) -- a standard file format
in crystallography -- is basically a restricted variant of STAR.
It is restricted in features (to make implementation easier),
but also imposes arbitrary limits -- for example on the line length.

DDL is a schema language for STAR/CIF.

All of them (STAR, CIF and DDL) have multiple versions.
We will be more specific in the following sections.

The STAR/CIF syntax is relatively simple, but may be confusing at first.
(Note that the initial version of STAR was published by Sydney Hall in 1991 --
before XML and much before JSON and YAML, not to mention TOML).

.. highlight:: default

Key-value pairs have the form::

    _year 2017  # means year = 2017

Blocks (sections) begin with the ``data_`` keyword
with a block name glued to it::

    data_tomato   # [tomato]
    _color red    # color = "red"

Importantly, unlike XML/JSON/YAML/TOML, STAR is designed
to concisely represent tabular data. This example from TOML docs::

    # TOML
    points = [ { x = 1, y = 2, z = 3 },
               { x = 7, y = 8, z = 9 },
               { x = 2, y = 4, z = 8 } ]

could be expressed as a so-called *loop* in STAR (keyword ``loop_``
followed by column names followed by values)::

    # STAR/CIF
    loop_
    _points.x _points.y _points.z
    1 2 3
    7 8 9
    2 4 8

Typically, long tables (loops) make most of the CIF content::

    1    N N   . LEU A 11  ? 0.5281 0.5618 0.5305 -0.0327 -0.0621 0.0104
    2    C CA  . LEU A 11  ? 0.5446 0.5722 0.5396 -0.0317 -0.0632 0.0080
    # many, many lines...
    5331 S SD  . MET C 238 ? 2.2952 2.3511 2.3275 -0.0895 0.0372  -0.0230
    5332 C CE  . MET C 238 ? 1.5699 1.6247 1.6108 -0.0907 0.0388  -0.0244

Note: "STAR File" is trademarked by IUCr, and it used to be patented_.

.. _patented: https://patents.google.com/patent/WO1991016682A1

The mmCIF (a.k.a. PDBx/mmCIF) format, which replaced the PDB format as
the primary format of the Protein Data Bank,
is the CIF syntax + a huge dictionary (ontology/schema) in DDL2.
The dictionary defines relations between columns in different tables,
which makes it resemble a relational database.
But it would be definitely not a normalized database - there is a lot
of redundancy in mmCIF.

**Where are the specs?**

International Tables for Crystallography
`Vol. G (2006) <http://it.iucr.org/Ga/contents/>`_
describes all of the STAR, CIF 1.1, DDL1 and DDL2.
If you don't have access to it -- IUCr website has specs of
`CIF1.1 <http://www.iucr.org/resources/cif/spec/version1.1>`_
and `DDLs <http://www.iucr.org/resources/cif/ddl>`_.
As far as I can tell all versions of the STAR spec are behind paywalls.

Later versions of the formats (hardly used as of 2017)
are described in these articles:
`STAR <https://dx.doi.org/10.1021/ci300074v>`_ (2012),
`DDLm <http://pubs.acs.org/doi/abs/10.1021/ci300075z>`_ and
`dREL <http://pubs.acs.org/doi/abs/10.1021/ci300076w>`_ (2012),
and `CIF 2.0 <http://journals.iucr.org/j/issues/2016/01/00/aj5269/>`_ (2016).
Only the last one is freely available.

mmCIF is documented at `mmcif.pdb.org <http://mmcif.pdb.org/>`_.

