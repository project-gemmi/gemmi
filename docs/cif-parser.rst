Gemmi CIF Parser
################

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


How to use it from C++
======================

The CIF parser is implemented in header files (``#include <gemmi/cif.hh>``),
so you do not need to compile Gemmi.
It has a single dependency: PEGTL (also header-only).
All you need is to make sure that Gemmi and PEGTL headers are in your
project's include path, and compile your program as C++11 or later.

If you'd like Gemmi to uncompress gzipped (.cif.gz) files on the fly,
add ``#include <gemmi/cifgz.hh>`` and link your program with the zlib library.

The CIF parsing functionality can be used in two ways that roughly
correspond to DOM and SAX parsing:

1. Parse a file (or C++ stream or string) into a document
   that can be easily accessed and manipulated.

2. Define own `PEGTL Actions <https://github.com/taocpp/PEGTL/blob/master/doc/Actions-and-States.md>`_
   corresponding to the grammar rules from ``cif.hh``.
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

We also keep track of comments, so we can write them back when saving file,
but they are in a separate list.

Let us start with a simple example.
This code reads mmCIF file and shows weights of the chemical components::

    #include <iostream>
    #include <gemmi/cif.hh>

    int main() {
      gemmi::cif::Document doc("1mru.cif");
      for (const gemmi::cif::Block& block : doc.blocks)
        for (const auto& cc : block.find("_chem_comp.", {"id", "formula_weight"}))
          std::cout << cc.as_str(0) << " weights " << cc.as_num(1) << std::endl;
    }

Reading a file
--------------

``struct gemmi::cif::Document`` has a few functions that populate it::

    void read_file(const std::string& filename)
    void read_memory(const char* data, const size_t size, const char* name)
    void read_cstream(std::FILE *f, size_t maximum, const char* name)
    void read_istream(std::istream &is, size_t maximum, const char* name)

Parameter ``name`` is used only when reporting errors.
Parameter ``maximum`` determines the buffer size and only affects performance.
Regardless of the buffer size, the last two options are slower
than ``read_file()`` -- they were not optimized for.

The constructor can take the filename as the argument::

    gemmi::cif::Document doc("1mru.cif")

which is equivalent to::

    gemmi::cif::Document doc;
    doc.read_file("1mru.cif");

Additional header ``cifgz.hh`` has a function::

    inline Document read_any(const std::string& path)

that can transparently open a gzipped file
(by uncompressing it first into a memory buffer).
Convenient when working with a local copy of the PDB archive::
    
    gemmi::cif::Document doc = read_any("mmCIF/pe/5pep.cif.gz");


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


How to use it from Python
=========================

You may get the project from github and compile the extension yourself,
or wait a few weeks and then it will be simply: ``pip install gemmi``.

Both Python 2 and 3 are supported.

Python bindings use `pybind11 <https://github.com/pybind/pybind11>`_.

.. highlight:: python

TODO: documentation

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

Many CIF readers do not store comments, as they are not part of the data.
In our case we may want to read, manipulate and save mmCIF file preserving
as much of the original file as possible. So we record the order of items,
comments and line breaks (but for now not any other whitespaces).

