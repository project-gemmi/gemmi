.. _cif_intro:

What are STAR, CIF, DDL, mmCIF?
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

The dot and question mark in the example above are two null types:
``?`` = *unknown* and ``.`` = *not applicable*.

The CIF syntax has a serious flaw resulting from historical trade-offs:
a string that can be interpreted as a number does not need to be quoted.
Therefore, the type of ``5332`` above is not certain:
the JSON equivalent can be either ``5332`` or ``"5332"``.

Note: "STAR File" is trademarked by IUCr, and it used to be patented_.

.. _patented: https://patents.google.com/patent/WO1991016682A1

The mmCIF (a.k.a. PDBx/mmCIF) format
is the CIF syntax + a huge dictionary (ontology/schema) in DDL2.
The dictionary defines relations between columns in different tables,
which makes it resemble a relational database.

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

PDBx/mmCIF is documented at `mmcif.pdb.org <http://mmcif.pdb.org/>`_.

What is parsed?
===============

The parser supports CIF 1.1 spec and some extras.

Currently, it is available as:

* C++11 header-only library, and
* Python (2 and 3) extension module.

We use it to read:

* mmCIF files
* Chemical Component Dictionary from PDB
* DDL1 and DDL2 dictionaries from IUCr and PDB
* monomer library a.k.a. Refmac dictionary
* a subset of COD

The parser handles:

* all constructs of CIF 1.1 (including *save frames*),
* the ``global_`` and ``stop_`` keywords from STAR -- needed for Refmac
  monomer library and ``mmcif_nmr-star.dic``, respectively.

It could be extended to handle also the new features of CIF 2.0
or even the full STAR format, but we don't have a good reason to do this.
The same goes for DDLm/dREL.

Additionally:

* names and lines can have any length like in STAR
  (the CIF spec imposes the limit of 2048 characters, but some mmCIF files
  from PDB exceed it, e.g. 3j3q.cif),
* unquoted strings cannot start with keywords (STAR spec is ambiguous
  about this -- see
  `StarTools doc <http://www.globalphasing.com/startools/>`_ for details)

Not supported (yet):

* the line wrapping convention (``eol;\eol``)

* Greek letters (``\m`` -> µ), accented letters (``\'o`` -> ó),
  special alphabetic characters (``\%A`` -> Å) and some other codes
  (``\\infty`` -> ∞) as per CIF 1.1 markup conventions

CIF parsers in the small-molecules field often face incorrect syntax.
The papers about `iotbx.cif <https://doi.org/10.1107/S0021889811041161>`_
and `COD::CIF::Parser <http://dx.doi.org/10.1107/S1600576715022396>`_
enumerate 12 and 10 common error categories, respectively.
As the MX community embraced the CIF format later, the problem may be less
severe for us. So we start with relatively strict parser and will be
pragmatically relaxing it when needed.


C++ Library
===========

The CIF parser is implemented in header files (``#include <gemmi/cif.hpp>``),
so you do not need to compile Gemmi.
It has a single dependency: PEGTL (also header-only),
which is included in the ``third_party`` directory.
All you need is to make sure that Gemmi and PEGTL headers are in your
project's include path, and compile your program as C++11 or later.
For example:

.. code-block:: none

    git clone https://github.com/project-gemmi/gemmi.git
    c++ -std=c++11 -Igemmi/include -Igemmi/third_party -O2 my_program.cpp

If you'd like Gemmi to uncompress gzipped (.cif.gz) files on the fly,
add ``#include <gemmi/gz.hpp>`` and link your program with the zlib library.

The CIF parsing functionality can be used in two ways that roughly
correspond to DOM and SAX parsing:

1. Parse a file (or C++ stream or string) into a document
   that can be easily accessed and manipulated.

2. Define own `PEGTL Actions <https://github.com/taocpp/PEGTL/blob/master/doc/Actions-and-States.md>`_
   corresponding to the grammar rules from ``cif.hpp``.
   These actions will be triggered while reading a CIF file.

This documention focuses on the former (DOM parsing).
The hierarchy in the DOM reflects the structure of CIF 1.1:

* Document contains blocks.
* Block can contain name-value pairs, loops and frames.
* Frame can contain name-value pairs and loops.
* Loop (*m*\ ×\ *n* table) contains *n* column names and *m*\ ×\ *n* values.

Names are often called *tags*. The leading ``_`` is usually treated
as part of the tag, not just a syntactic feature, so we store it in DOM
as ``std::string`` with value ``_my_tag``.

Values have types, but while the type can be inferred correctly in almost all
cases, it is the corresponding dictionary that specifies the type.
Additionally, DDL2 dictionaries can specify subtypes of the standard CIF types,
for example ``int`` and ``float`` are mmCIF subtypes of a generic CIF ``numb``.

Let us start with a simple example.
This code reads mmCIF file and shows weights of the chemical components::

    #include <iostream>
    #include <gemmi/cif.hpp>

    int main() {
      gemmi::cif::Document doc("1mru.cif");
      for (const gemmi::cif::Block& block : doc.blocks)
        for (const auto& cc : block.find("_chem_comp.", {"id", "formula_weight"}))
          std::cout << cc[0] << " weights " << cc[1] << std::endl;
    }

Reading a file
--------------

The ``gemmi::cif`` namespace has a few functions that return Document::

    Document read_file(const std::string& filename)
    Document read_memory(const char* data, const size_t size, const char* name)
    Document read_cstream(std::FILE *f, size_t maximum, const char* name)
    Document read_istream(std::istream &is, size_t maximum, const char* name)

Parameter ``name`` is used only when reporting errors.
Parameter ``maximum`` determines the buffer size and only affects performance.
Regardless of the buffer size, the last two options are slower
than ``read_file()`` -- they were not optimized for.

Additional header ``<gemmi/gz.hpp>`` is needed to transparently open
a gzipped file (by uncompressing it first into a memory buffer)
if the filename ends with ``.gz``::

    gemmi::cif::Document doc = gemmi::cif::read(gemmi::MaybeGzipped(path));

And if the ``path`` above is ``-``, the standard input is read.

Document
--------

``Document`` has the following variables that can be accessed directly::

  std::string source;  // filename or the explicitly provided name
  std::vector<Block> blocks;

As the mmCIF files are expected to have only a single block,
we have a function::

  const Block& sole_block() const;

to express the intention of accessing the only block in the file
(it throws an exception if the number of blocks is not one).

At last, function ``clear()`` works in the same way as in C++ containers.

Block
-----

.. warning::
    The API still evolves and for now this documentation lists only
    the most used functions.

Value corresponding to a particular tag can is read using::

    const std::string* find_value(const std::string& tag) const;

which returns ``nullptr`` if there is no such tag in the block.
The result is a raw string (possibly with quotes) that can be fed into
``as_string()`` or ``as_number()`` or ``as_int()``.
For example::

    const std::string *rf = block.find_value("_refine.ls_R_factor_R_free");
    // here you may check for rf == nullptr, possibly also for "?" and "."
    double rfree = gemmi::cif::as_number(*rf); // NaN if '?' or '.'

If you do not need to distinguish between missing tag, unknown (?)
and n/a (.), we have convenience functions::

    // returns empty string if not found or unknown or n/a.
    const std::string& find_string(const std::string& tag) const;
    // returns NaN if not found or unknown or n/a, throws if not numeric
    double find_number(const std::string& tag) const;
    // returns default_ if not found or unknown or n/a, throws if not int
    int find_int(const std::string& tag, int default_) const;

To read values from a single column for a loop (table) use::

    LoopColumn find_loop(const std::string& tag) const;

The values can be iterated over using a C++11 range-based ``for``::

    for (const std::string &s : block.find_loop("_atom_site.type_symbol"))
      std::cout << gemmi::cif::as_string(s) << std::endl;

TODO: document LoopColumn

Most often, we want to access multiple (but not necessarily all) columns
from a table.
Additionally, some values can be given either in a loop or, if the loop
would have only a single row, as tag-value pairs.
So we want our access function to handle transparently both cases.
These requirements led to a functions ``find``::

    TableView find(const std::string& tag) const;
    TableView find(const std::vector<std::string>& tags) const;

which returns a lightweight, iterable (by C++11 range-based ``for``) view
of the data.

As a rule, columns from the same loop have a common prefix,
so we added a third overload::

    TableView find(const std::string& prefix,
                   const std::vector<std::string>& tags) const;

so we can write::

    block.find("_entity_poly_seq.", {"entity_id", "num", "mon_id"})

instead of::

    block.find({"_entity_poly_seq.entity_id", "_entity_poly_seq.num", "_entity_poly_seq.mon_id"})

TODO: document TableView methods (``ok()``, ``width()``, ``length()``,
``operator[](size_t)``, ``at(size_t)``, ``find_row(const std::string&)``.

The first example in this section shows how this function can be used.


Python Module
=============

.. highlight:: python

Both Python 2.7 and 3.x are supported (thanks to the excellent
`pybind11 <https://github.com/pybind/pybind11>`_ project).
To install the gemmi module you need pip, git and not too old
C++ compiler (GCC 4.8+, Clang 3.4+, MSVC 2015+, ICC 16+)::

    pip install git+https://github.com/project-gemmi/gemmi.git

(when the project is more mature and has regular releases, it will be simply
``pip install gemmi``).

TODO: documentation

(in the meantime see ``pydoc gemmi.cif`` and examples)

Example (says hello to each element found in mmCIF):

.. literalinclude:: ../examples/hello.py
   :lines: 2-

More complex examples are shown in the :ref:`cif_examples` section.

Performance
===========

We owe the good performance to the excellent
`PEGTL <https://github.com/taocpp/PEGTL/>`_ project.

In my testing (with GCC 5 and Clang 3.8) Gemmi CIF parser is
3x faster than cif_api (``validate -f`` vs ``cif2_syncheck -f``),
which in turn `is reported <https://doi.org/10.1107/S1600576715021883>`_
to be several times faster than iotbx.cif (ucif).

On the other hand the ChimeraX readcif library
(which is not publicly available?) is likely even faster.
`This benchmark <http://www.cgl.ucsf.edu/chimerax/docs/devel/core/atomic/readcif_cpp/docs/compare.html>`_
reports that readcif run in the "tokenized" mode reads
3j3q.cif (250MB) in 1.8 sec (20x faster than cifparse-obj and ucif).
On my computer (with similar spec) Gemmi parses the same file in <2s
in the validation-only mode,
but in >4s when copying all the strings into a DOM structure.
Doing the same what that benchmark does should be somewhere between 2 and 4s.

While big (10x) differences between programs may be surprising,
it is the same with
`JSON parsers <https://github.com/miloyip/nativejson-benchmark>`_,
and in many cases it does not matter.


Design rationale
================

Parser
------

Parsing formal languages is a well-researched topic in computer science.
The first versions of lex and yacc - popular tools that generate lexical
analyzers and parsers - were written in 1970's. Today many tools exist
to translate grammar rules into C/C++ code.
On the other hand many compilers and high-profile tools use hand-coded
parsers - as it is more flexible.

Looking at other STAR and CIF parsers -- some use parser generators
(COD::CIF::Parser and starlib2 use yacc/bison, iotbx.cif uses Antlr3),
others are hand-coded.

I had experience with flex/bison and Boost.Spirit
(and I wanted to try also Lemon and re2c)
but I decided to use PEGTL for this task. I was convinced by the
`TAOC++ JSON <https://github.com/taocpp/json>`_
parser that is based on PEGTL and has a good balance of simplicity
and performance.

PEGTL is a C++ library (not a generator) for creating PEG parsers.
PEG stands for Parsing Expression Grammar -- a simpler approach than
tradional Context Free Grammar.

As a result, our parser depends on a third-party (header-only) library,
but the parser itself is pretty simple.

Data structures
---------------

The next thing is how we store the data read from file.
We decided rely on the C++ standard library where we can.

Generally, storage of such data involves (in C++ terms) some containers
and a union/variant type for storing values of mixed types.

We use primarily ``std::vector`` as a container,
and ``std::unordered_map`` when quicker access is needed.

Custom structures with (unrestricted) unions are used where variants
are needed.

Strings are stored in ``std::string`` and it is fast enough.
Mainstream C++ standard libraries have short string optimization (SSO)
for up to 15 or 22 characters, which covers most of the values in mmCIF files.

