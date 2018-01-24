CIF Parser
##########

.. _cif_intro:

What are STAR, CIF, DDL, mmCIF?
==================================

(in case someone comes here when looking for a serialization format)

STAR is a human-readable data serialization format (think XML or JSON)
that happens to be known and used only in molecular-structure sciences.

CIF (Crystallographic Information File) -- a file format used
in crystallography -- is a restricted derivative of STAR.
It is restricted in features (to make implementation easier),
but also imposes arbitrary limits -- for example on the line length.

DDL is a schema language for STAR/CIF.

All of them (STAR, CIF and DDL) have multiple versions.
We will be more specific in the following sections.

The STAR/CIF syntax is relatively simple, but may be confusing at first.
(Note that the initial version of STAR was published by Sydney Hall in 1991 --
before XML and long before JSON and YAML, not to mention TOML).

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

The dot and question mark in the example above are two null types.
In the CIF spec: ``?`` = *unknown* and ``.`` = *not applicable*.
In mmCIF files ``.`` is used for mandatory items, ``?`` for not mandatory.

The CIF syntax has a serious flaw resulting from historical trade-offs:
a string that can be interpreted as a number does not need to be quoted.
Therefore, the type of ``5332`` above is not certain:
the JSON equivalent can be either ``5332`` or ``"5332"``.

Note: "STAR File" is trademarked by IUCr, and it used to be patented_.

.. _patented: https://patents.google.com/patent/WO1991016682A1

The mmCIF format (by mmCIF we mean what is more formally called PDBx/mmCIF)
is the CIF syntax + a huge dictionary (ontology/schema) in DDL2.
The dictionary defines relations between columns in different tables,
which makes it resemble a relational database (it was designed at the
height of popularity of RDBMSs).

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

* mmCIF files (both coordinates and structure factors)
* CIF files from Crystallography Open Database (COD)
* Chemical Component Dictionary from PDB
* DDL1 and DDL2 dictionaries from IUCr and PDB
* monomer library a.k.a. Refmac dictionary

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
* quoted strings may contain non-ascii characters (if nothing has changed
  one entry in the PDB has byte A0 corresponding to non-breaking space)
* unquoted strings cannot start with keywords (STAR spec is ambiguous
  about this -- see
  `StarTools doc <http://www.globalphasing.com/startools/>`_ for details)

The parser does *not* interpret automatically CIF 1.1 conventions
(that work around syntactic limitations):

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


DOM and SAX
===========

The CIF parsing functionality can be used in two ways that roughly
correspond to DOM and SAX parsing:

* The usual way is to parse a file or a string into a document (DOM)
  that can be easily accessed and manipulated.

* Alternatively, from C++ only, one can define own
  `PEGTL Actions <https://github.com/taocpp/PEGTL/blob/master/doc/Actions-and-States.md>`_
  corresponding to the grammar rules from ``cif.hpp``.
  These actions will be triggered while reading a CIF file.

This documention focuses on the former (DOM parsing).
The hierarchy in the DOM reflects the structure of CIF 1.1:

* Document contains blocks.
* Block can contain name-value pairs, loops and frames.
* Frame can contain name-value pairs and loops.
* Loop (*m*\ ×\ *n* table) contains *n* column names and *m*\ ×\ *n* values.

Names are often called *tags*. The leading ``_`` is usually treated
as part of the tag, not just a syntactic feature. So we store tag string
with the underscore (``_my_tag``), and function that take tags as arguments
expect strings starting with ``_``.

Values have types. CIF 1.1 defines four base types:

* char (string)
* uchar (ughh.. case-insensitive string)
* numb (number that cannot be recognized as number on the syntax level)
* null (one of two possible nulls: ``?`` and ``.``)

Additionally, DDL2 dictionaries specify subtypes.
For example, ``int`` and ``float`` are mmCIF subtypes of ``numb``.

Since in general it is not possible to infer numeric type without
a corresponding dictionary, the numbers are stored in the DOM as strings
and they are converted to numeric values later, if needed.

The case of all strings is preserved in the DOM.

C++ Library
===========

.. highlight:: cpp

CIF parser is implemented in header files,
so you do not need to compile Gemmi.
It has a single dependency: PEGTL (also header-only),
which is included in the ``third_party`` directory.
All you need is to make sure that Gemmi and PEGTL headers are in your
project's include path, and compile your program as C++11 or later.

Let us start with a simple example.
This little program reads mmCIF file and shows weights of the chemical
components:

.. literalinclude:: doc_cif_cc.cpp
   :language: cpp
   :lines: 1-11

To compile it on Unix system you need to fetch Gemmi source code
and run a compiler:

.. code-block:: none

    git clone https://github.com/project-gemmi/gemmi.git
    c++ -std=c++11 -Igemmi/include -Igemmi/third_party -O2 my_program.cpp

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

    // in this and all the next examples: namespace cif = gemmi::cif;
    cif::Document doc = cif::read(gemmi::MaybeGzipped(path));

If you use it, you must also link the program with zlib. On Unix systems
it usually means adding ``-lz`` to the compiler invocation.

And if the ``path`` above is ``-``, the standard input is read.

Document
--------

``Document`` has the two member variables that can be accessed directly::

  std::string source;  // filename or the name passed to read_memory()
  std::vector<Block> blocks;

Each ``Block`` corresponds to a data block in the CIF file and must have
a unique name. To access a block with known name use::

  Block* find_block(const std::string& name);

As the coordinate mmCIF files are expected to have only a single block,
we have a function::

  const Block& sole_block() const;

to express the intention of accessing the only block in the file
(it throws an exception if the number of blocks is not one).

.. note::

    We first cover only accessing the data. Making changes is described
    later in the "Editing" subsection.

Block
-----

Each block has a name and a list of items::

    std::string name;
    std::vector<Item> items;

Item contains an unrestricted (C++11) union that holds one of:

* name-value pair (Pair),
* table/loop (Loop)
* or save frame (also ``struct Block``).

Pair
~~~~

Name-value pairs from the CIF file are stored as Pair::

    using Pair = std::array<std::string, 2>;

A pair with a particular tag can be located using::

    const Pair* find_pair(const std::string& tag) const;

As a short-cut, the value of a pair can be read with function::

    const std::string* find_value(const std::string& tag) const;

Both functions return ``nullptr`` if the tag is not found (and they
do not search in CIF loops).

The value string is a raw string (possibly with quotes) that can be fed into
``as_string()`` or ``as_number()`` or ``as_int()``.

For example::

    const std::string *rf = block.find_value("_refine_ls_R_factor_all");
    // here we should check for rf == nullptr, possibly also for "?" and "."
    double rfree = cif::as_number(*rf); // NaN if '?' or '.'

Loop
~~~~

To read values from a single column for a loop (table) use::

    Column find_loop(const std::string& tag);

``Column`` is a lightweight class with a few functions::

    // Number of rows in the loop. 0 means that the tag was not found.
    int length() const;

    // Returns pointer to the column name in the DOM.
    // The name is the same as argument to find_loop() or find_values().
    std::string* get_tag();

    // Returns pointer to the DOM structure containing the whole table.
    Loop* get_loop() const;

    // Get raw value (no bounds checking).
    std::string& operator[](int n);

    // Get raw value (after bounds checking).
    std::string& at(int n);

    // Short-cut for cif::as_string(column.at(n)).
    std::string str(int n) const;

``Column`` also provides support for C++11 range-based ``for``::

    for (const std::string &s : block.find_loop("_atom_site_type_symbol"))
      std::cout << cif::as_string(s) << std::endl;

Pair or Loop
~~~~~~~~~~~~

Since some values (and in mmCIF files -- all values) can be given
either as a single item or in a loop, it is convenient to handle
transparently both cases. To to this, just replace ``find_loop`` with
a similar function::

    Column find_values(const std::string& tag);

For example::

    // mmCIF _chem_comp_atom is usually a table, but not always
    for (const std::string &s : block.find_values("_chem_comp_atom.type_symbol"))
      std::cout << s << std::endl;

The ``find_values`` function also returns ``struct Column``,
but unlike ``find_loop`` it finds name-value pairs as well.
In such case ``Column::get_length()`` returns 1
and ``Column::get_loop()`` return ``nullptr``.


Often, we want to access multiple columns at once,
so the library has another abstraction (``Table``)
that can be used with multiple tags
(this function was used in the first example in this section)::

    Table find(const std::vector<std::string>& tags);

Since columns from the same loop tend to have common prefix (category name),
the library provides a second form::

    Table find(const std::string& prefix, const std::vector<std::string>& tags);

These two calls are equivalent::

    block.find({"_entity_poly_seq.entity_id", "_entity_poly_seq.num", "_entity_poly_seq.mon_id"})

    block.find("_entity_poly_seq.", {"entity_id", "num", "mon_id"})

Note that ``find`` is not aware of dictionaries and categories,
therefore the category name should end with a separator
(dot for mmCIF files, as shown above).

If one or more of the tags is not found in the block, the returned
``Table`` is empty.

TODO: document optional tags: ``{"_required_tag", "?_optional_tag"}``

The ``Table`` above is, like ``Column``,  a lightweight, iterable view
of the data. It has functions to check the shape of the table::

    bool ok() const;  // true if the table is not empty
    size_t length() const;  // number of rows
    size_t width() const;  // number of values in each row
    bool has_column(int n) const; // for use with optional columns

Most importantly, it provides access to rows (``Table::Row``)
that in turn provide access to value strings (``std::string``)::

    Row operator[](int n);  // access Row
    Row at(int n);  // the same but with bounds checking
    // and also begin() and end() to work in range-for.

as well as to the tags::

    Row tags();  // pseudo-row that contains tags

Additionally, it has two convenience functions to access rows::

    // Returns the first row that has the specified string in the first column.
    Row find_row(const std::string& s);

    // Make sure that the table has only one row and return it.
    Row one();

``Table::Row`` has as few basic functions::

    size_t size() const;  // the width of the table

    // Get raw value.
    std::string& operator[](int n);  // no bounds checking
    std::string& at(int n);          // with bounds checking

    std::string str(int n) const; // short-cut for cif::as_string(row.at(n))
    bool has(int n) const; // the same as Table::has_column(n)

and also supports iterators.

It is also possible to work with the table column-wise::

    Column column(int n);
    // alternatively, specify the full tag name or its unique ending 
    Column find_column(const std::string& suffix);

Both ``Column`` and ``Table::Row`` have functions ``begin()`` and ``end()``
in const and non-const variants, returning ``iterator`` and
``const_iterator`` types, respectively. These types satisfy requirements
of the BidirectionalIterator concept.
Conversely, the iterator over the rows of ``Table`` is a minimalistic
structure -- just enough get the range-for work.

As an example, let us convert mmCIF
to the `XYZ format <https://en.wikipedia.org/wiki/XYZ_file_format>`_:

.. literalinclude:: doc_cif_cc.cpp
   :language: cpp
   :lines: 15-


mmCIF categories
~~~~~~~~~~~~~~~~

mmCIF files group data into categories. All mmCIF tags have a dot
(e.g. ``_entry.id``) and the category name is the part before the dot.

We have two functions to work with categories::

    std::vector<std::string> get_mmcif_category_names() const;
    Table find_mmcif_category(std::string cat);

The first one returns a list of all categories in the block.
The second one returns a Table described above, with all tags belonging
the specified category.

Frame
~~~~~

The named save frames (keyword ``save_``) from the STAR specification
are used in CIF files only as sub-sections of a block.
The only place where they are enountered are mmCIF dictionaries.
Because of this, we expect that few users will ever need to access
frames and we do not have special function on this occasion.
Still, one can easily work with save frames::

    for (cif::Item& item : block.items)
      if (item.type == cif::ItemType::Frame)
        // doing something with item.frame which is a (nested) Block
        cif::Block& frame = item.frame;

A full example is in the header ``gemmi/ddl.hpp``.

Editing
-------

TODO: describe how to create a new cif::Document and modify existing one,
(changing values using iterators,
Document::clear(), Column::item().erase(), Table::erase(),
Block::set_pair(), Block::init_loop(), ...)

Writing
-------

The functions writing ``cif::Document`` to C++ stream or to a file
are in a separate header file ``gemmi/to_cif.hpp``::

    void write_to_file(const Document& doc, const std::string& filename,
                       Style style=Style::Simple);

    void write_out_document(std::ostream& os, const Document& doc, Style style);

    // the same as write_out_document() with Style::Simple
    std::ostream& operator<<(std::ostream& os, const gemmi::cif::Document& doc);

Three "styles" are supported:

* ``Style::Simple`` simply writes out the DOM structure,
* ``Style::PreferPairs`` writes single-row loops as pairs,
* ``Style::Pdbx`` additionally puts ``#`` (empty comments) between categories,
  mimicking the official PDBx/mmCIF files from wwPDB.

JSON
----

Header ``gemmi/to_json.hpp`` provides code for serializing
``cif::Document`` as JSON. It has a number of options for customizing
the translation. In particular, both mmCIF and CIF-JSON flavours are supported.
More details about the flavours are given in the description of
:ref:`gemmi-convert <json>`.

Such JSON files can be read back into the ``cif::Document`` structure
using function from ``gemmi/json.hpp``.

Python Module
=============

.. highlight:: python

Python module for Python 2.7 and 3.x
can be installed with pip, as described in the
:ref:`Installation <install_py>` section.
After installation ``pydoc gemmi.cif`` should list all classes and methods.

To start with a simple example, here is a program that says hello to each
element found in mmCIF:

.. literalinclude:: ../examples/hello.py
   :lines: 2-

More complex examples are shown in the :ref:`cif_examples` section.

Internally, Python bindings use
`pybind11 <https://github.com/pybind/pybind11>`_.

Reading a file
--------------

The content of a CIF file is stored in class ``Document``.

.. testcode::

  from gemmi import cif
  # read and parse a CIF file
  doc = cif.read_file("components.cif")
  # the same, but if the filename ends with .gz it is uncompressed on the fly
  doc = cif.read("../tests/1pfe.cif.gz")
  # read content of a CIF file from string
  doc = cif.read_string("data_this _is valid _cif content")

Document
--------

``Document`` contains blocks with data.
It can be iterated, accessed by block index and by block name:

.. doctest::

  >>> doc = cif.read_file("components.cif")
  >>> len(doc)  #doctest: +SKIP
  25219
  >>> doc[0]
  <gemmi.cif.Block 000>
  >>> doc[-1]
  <gemmi.cif.Block ZZZ>
  >>> doc['MSE']
  <gemmi.cif.Block MSE>

  >>> # The function block.find_block(name) is like block[name] ...
  >>> doc.find_block('MSE')
  <gemmi.cif.Block MSE>
  >>> # ... except when the block is not found:
  >>> doc.find_block('no such thing')  # -> None
  >>> # doc['no such thing'] # -> KeyError

  # As many cif files (including mmCIF) have only a single block,
  # we have a dedicated function that check if it is indeed a file with
  # a single block (throws exception if not) and returns the block
  >>> cif.read("../tests/1pfe.cif.gz").sole_block()
  <gemmi.cif.Block 1PFE>

Functions that modify the content of the block are described in the section
"Editing" below.

.. note::

    We first cover only accessing the data. Making changes and writing
    CIF files is described later in the "Editing" subsection.

Block
-----

Each block has a name:

.. doctest::

  >>> doc = cif.read("../tests/1pfe.cif.gz")
  >>> block = doc.sole_block()
  >>> block.name
  '1PFE'

and a list of items. Each item corresponds to one of:

* name-value pair (``cif.Pair``)
* loop/table (``cif.Table``)
* or save frame (also an instance of ``cif.Block``).

We have functions to access name-value pairs:

.. doctest::

  >>> block.find_pair('_cell.length_a')
  ['_cell.length_a', '39.374']
  >>> block.find_value('_cell.length_b')
  '39.374'
  >>> print(block.find_value('_cell.length_y'))
  None

and values in loop:

.. doctest::

  >>> block.find_loop('_atom_type.symbol')
  <gemmi.cif.Column _atom_type.symbol length 6>
  >>> list(_)
  ['C', 'CL', 'N', 'O', 'P', 'S']

but it is often better to use functions that work with both pairs and loops:

.. doctest::

  >>> block.find_values('_cell.length_a')
  <gemmi.cif.Column _cell.length_a length 1>
  >>> block.find_values('_atom_type.symbol')
  <gemmi.cif.Column _atom_type.symbol length 6>

Both ``find_loop()`` and ``find_values()`` return ``Column``.
Column's special method ``__bool__`` tells if the tag was found.
``__len__`` returns the number of corresponding values.
``__iter__``,  ``__getitem__`` and ``__setitem__``
get or set a raw string (i.e. string with quotes, if applicable).
To get the actual string content one may use the method ``str``:

.. doctest::

  >>> column = block.find_values('_chem_comp.formula')
  >>> column[7]
  "'H2 O'"
  >>> cif.as_string(column[7])
  'H2 O'
  >>> column.str(7)  # the same with less writing
  'H2 O'

If the tag was found in a loop, method ``get_loop`` returns a reference
to this ``Loop`` in the DOM. Otherwise it returns ``None``.

Table
~~~~~

Often, we want to access multiple columns at once,
so the library has another abstraction (``Table``)
that can be used with multiple tags:

.. doctest::

  >>> block.find(['_entity_poly_seq.entity_id', '_entity_poly_seq.num', '_entity_poly_seq.mon_id'])
  <gemmi.cif.Table 18 x 3>

Since tags in one loop tend to have a common prefix (category name),
the library provides a second form that takes the common prefix as the first
argument:

.. doctest::

  >>> block.find('_entity_poly_seq.', ['entity_id', 'num', 'mon_id'])
  <gemmi.cif.Table 18 x 3>

**TODO**: document Table and Row

mmCIF categories
~~~~~~~~~~~~~~~~

mmCIF files group data into categories. All mmCIF tags have a dot
(e.g. ``_entry.id``) and the category name is the part before the dot.

.. doctest::

  >>> block.get_mmcif_category_names()[:3]
  ['_entry.', '_audit_conform.', '_database_2.']
  >>> block.find_mmcif_category('_entry.')
  <gemmi.cif.Table 1 x 1>
  >>> _.tags[0], _[0][0]
  ('_entry.id', '1PFE')
  >>> cat = block.find_mmcif_category('_database_2.')
  >>> cat
  <gemmi.cif.Table 4 x 2>
  >>> len(cat.tags)
  2
  >>> list(cat.tags)
  ['_database_2.database_id', '_database_2.database_code']
  >>> for row in cat:
  ...   print('%s: %s' % tuple(row))
  PDB: 1PFE
  NDB: DD0057
  RCSB: RCSB019291
  WWPDB: D_1000019291
  >>> cat[3][1]
  'D_1000019291'

Editing
-------

Changes in-place
~~~~~~~~~~~~~~~~

Both Column and Row implement ``__setitem__`` that can be used to change
raw string values.
As an example, we can shift all the atoms by 1A in the x direction:

.. doctest::

  >>> col_x = block.find_values('_atom_site.Cartn_x')
  >>> for n, x in enumerate(col_x):
  ...   col_x[n] = str(float(x) + 1)

or swap two tag names (these two tend to have identical values, so no one
will notice):

.. doctest::

  >>> tags = block.find('_atom_site.', ['label_atom_id', 'auth_atom_id']).tags
  >>> tags[0], tags[1] = tags[1], tags[0]

Pairs
~~~~~

To add a name-value pair, replacing current item if it exists,
use function ``set_pair``:

.. doctest::

  >>> block.set_pair('_tree', 'willow')

Loops
~~~~~

To add a row to an existing table (loop) use ``add_row``:

.. doctest::

  >>> loop = block.find_loop('_atom_type.symbol').get_loop()
  >>> loop.add_row(['Au'], pos=0)
  >>> loop.add_row(['Zr'])  # appends
  >>> list(block.find_loop('_atom_type.symbol'))
  ['Au', 'C', 'CL', 'N', 'O', 'P', 'S', 'Zr']

``set_all_values`` sets all the data in a table. It takes as an argument
a list of lists of string. The lists of strings correspond to columns.

.. doctest::

  >>> loop = block.find_loop('_citation_author.citation_id').get_loop()
  >>> loop.tags
  ['_citation_author.citation_id', '_citation_author.name', '_citation_author.ordinal']
  >>> loop.set_all_values([['primary']*2, ['Alice', 'Bob'], ['1', '2']])
  >>> for row in block.find_mmcif_category('_citation_author.'):
  ...   print(list(row))
  ['primary', 'Alice', '1']
  ['primary', 'Bob', '2']
      

To add a new table (replacing old one if it exists) use ``init_loop``:

.. doctest::

  >>> loop = block.init_loop('_ocean_', ['id', 'name'])
  >>> # empty table is invalid in CIF, we need to add something
  >>> loop.add_row(['1', 'Atlantic'])

In the above example, if the block already has tags ``_ocean_id``
and/or ``_ocean_name`` and

* if they are in a table: the table will be cleared and re-used,
* if they are in name-value pairs: the pairs will be removed
  and a table will be created at the position of the first pair.

mmCIF
~~~~~

We also have a variant of ``init_loop`` for working with mmCIF categories:

.. doctest::

  >>> loop = block.init_mmcif_loop('_ocean.', ['id', 'name'])
  >>> loop.add_row(['1', 'Atlantic'])

The subtle difference is that if the block has other name-value pairs
in the same category (say, ``_ocean.depth 8.5``)
``init_loop`` leaves them untouched, but ``init_mmcif_loop`` removes them.

Additionally, like in other ``_mmcif_`` functions, the trailing dot
in the category name may be omitted (but the leading underscore is required).

Python programs often keep mmCIF data in dictionaries.
Thus we have a function that takes a category name and a dictionary
mapping tags to lists, and creates (or replaces) a corresponding category:

.. doctest::

  >>> ocean_data = {'id': [1, 2, 3], 'name': ['Atlantic', 'Pacific', 'Indian']}
  >>> block.set_mmcif_category('_ocean', ocean_data)

This function has a third parameter ``raw`` with default value False.
Specify ``raw=True`` if the values are already quoted, otherwise
they will be quoted twice.

.. doctest::

  >>> block.set_mmcif_category('_raw', {'a': ['?', '.', "'a b'"]}, raw=True)
  >>> list(block.find_values('_raw.a'))
  ['?', '.', "'a b'"]
  >>> block.set_mmcif_category('_nop', {'a': ['?', '.', "'a b'"]}, raw=False)
  >>> list(block.find_values('_nop.a'))
  ["'?'", "'.'", '"\'a b\'"']
  >>> block.set_mmcif_category('_ok', {'a': [None, False, 'a b']}, raw=False)
  >>> list(block.find_values('_ok.a'))
  ['?', '.', "'a b'"]

Creating and deleting
~~~~~~~~~~~~~~~~~~~~~

If needed, one can remove block from a document (``Document.__delitem__``),
or delete pair or table from a block (``Table.erase()``).
An example for this is provided in the :ref:`CCD section <ccd_example>` below.

Lastly, the CIF file can be created from scratch.

.. doctest::

  >>> d = cif.Document()
  >>> d.add_new_block('oak')
  <gemmi.cif.Block oak>
  >>> _.set_pair('_nut', 'acorn')
  >>> print(d.as_string().strip())
  data_oak
  _nut acorn


Writing
-------

The modifications can be written using function ``Document.write_file(path)``.

.. doctest::

  >>> doc.write_file('1pfe-modified.cif')

The formatting of the output can be customized to some degree.
Currently, we have an option to write single-row loops as pairs:

.. doctest::

  >>> doc.write_file('1pfe-modified.cif', cif.Style.PreferPairs)


and to, additionally, put ``#`` (empty comments) between mmCIF categories:

.. doctest::

  >>> doc.write_file('1pfe-styled.cif', cif.Style.Pdbx)

This mimicks the peculiar formatting of PDBx/mmCIF files in the official
wwPDB archive, and it is helpful when diff-ing original and modified
files (you also need diff option ``--ignore-space-change``).

Other functions:

* ``Document.as_string()`` -- returns the text that
  would be written by ``write_file()``,
* ``Document.as_json()`` -- returns the document serialized to JSON string.

Performance
===========

Gemmi CIF parser is based on `PEGTL <https://github.com/taocpp/PEGTL/>`_ --
an excellent "parser combinator library".

In my testing (with GCC 5 and Clang 3.8) our parser is
3x faster than cif_api (``validate -f`` vs ``cif2_syncheck -f``),
which in turn `is reported <https://doi.org/10.1107/S1600576715021883>`_
to be several times faster than iotbx.cif (ucif).

While further improvement would be possible (the fastest JSON parsers are
`much faster <https://github.com/project-gemmi/benchmarking-json>`_),
it is not a priority.

Directory walking
=================

Many of the utilities and examples developed for this project
work with archives of CIF files such as wwPDB or COD.
To make it easier to iterate over all CIF files in a directory tree
we provide a C++11 ``CifWalk`` class with a simple API:

.. code-block:: cpp

  #include <gemmi/dirwalk.hpp>

  // ...
  // throws std::runtime_error if top_dir doesn't exist
  for (const char* cif_file : gemmi::CifWalk(top_dir)) {
    cif::Document doc = cif::read(gemmi::MaybeGzipped(cif_file));
    // ...
  }

This header file contains also a more general ``DirWalk`` class.
Both are based on the `tinydir <https://github.com/cxong/tinydir>`_ library.

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
We decided to rely on the C++ standard library where we can.

Generally, storage of such data involves (in C++ terms) some containers
and a union/variant type for storing values of mixed types.

We use primarily ``std::vector`` as a container,
and ``std::unordered_map`` when quicker access is needed.

A custom structure ``Item`` with (unrestricted) union is used as
a variant-like class.

Strings are stored in ``std::string`` and it is fast enough.
Mainstream C++ standard libraries have short string optimization (SSO)
for up to 15 or 22 characters, which covers most of the values in mmCIF files.


.. _cif_examples:

Examples
========

The examples here use C++11 or Python 2.7/3.x.
Full working code code can be found in the examples__ directory.

If you have the Python ``gemmi`` module installed you should also have
the Python examples -- try ``python -m gemmi-examples`` if not sure
where they are.

The examples below can be run on one or more PDBx/mmCIF files.
We normally run them on a
`local copy <https://www.wwpdb.org/download/downloads>`_ of the mmCIF
archive (30GB+ gzipped, don't uncompress!) to perform PDB-wide analyses,

__ https://github.com/project-gemmi/gemmi/tree/master/examples

.. _auth_label_example:

auth vs label
-------------

When you look at the list of atoms (``_atom_site.*``) in mmCIF files
some columns seem to be completely redundant. Are they?

.. raw:: html

    <table style="font-size:11px; border-spacing:5px 0; margin-top:60px; width:100%; border-collapse:separate;">
    <colgroup>
    <col span="3">
    <col style="background-color:#fc8">
    <col>
    <col style="background-color:#ffa">
    <col style="background-color:#aff">
    <col>
    <col style="background-color:#bd8">
    <col span="12">
    <col style="background-color:#bd8">
    <col style="background-color:#ffa">
    <col style="background-color:#aff">
    <col style="background-color:#fc8">
    <col>
    </colgroup>
    <tr style="font-size:10px;">
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">group_PDB</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">id</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">type_symbol</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">label_atom_id</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">label_alt_id</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">label_comp_id</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">label_asym_id</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">label_entity_id</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">label_seq_id</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">pdbx_PDB_ins_code</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">Cartn_x</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">Cartn_y</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">Cartn_z</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">occupancy</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">B_iso_or_equiv</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">Cartn_x_esd</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">Cartn_y_esd</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">Cartn_z_esd</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">occupancy_esd</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">B_iso_or_equiv_esd</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">pdbx_formal_charge</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">auth_seq_id</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">auth_comp_id</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">auth_asym_id</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">auth_atom_id</th>
    <th style="max-width:12px; transform-origin:left; transform:rotate(315deg);">pdbx_PDB_model_num</th>
    <tr>
    <td>ATOM</td><td>1</td><td>N</td><td>N</td><td>.</td><td>GLY</td><td>A</td>
    <td>1</td><td>1</td><td>?</td><td>-9.009</td><td>4.612</td><td>6.102</td>
    <td>1.00</td><td>16.77</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>1</td><td>GLY</td><td>A</td><td>N</td><td>1</td>
    </tr><tr>
    <td>ATOM</td><td>2</td><td>C</td><td>CA</td><td>.</td><td>GLY</td><td>A</td>
    <td>1</td><td>1</td><td>?</td><td>-9.052</td><td>4.207</td><td>4.651</td>
    <td>1.00</td><td>16.57</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>1</td><td>GLY</td><td>A</td><td>CA</td><td>1</td>
    </tr><tr>
    <td>ATOM</td><td>3</td><td>C</td><td>C</td><td>.</td><td>GLY</td><td>A</td>
    <td>1</td><td>1</td><td>?</td><td>-8.015</td><td>3.140</td><td>4.419</td>
    <td>1.00</td><td>16.16</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>1</td><td>GLY</td><td>A</td><td>C</td><td>1</td>
    </tr><tr>
    <td>ATOM</td><td>4</td><td>O</td><td>O</td><td>.</td><td>GLY</td><td>A</td>
    <td>1</td><td>1</td><td>?</td><td>-7.523</td><td>2.521</td><td>5.381</td>
    <td>1.00</td><td>16.78</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>1</td><td>GLY</td><td>A</td><td>O</td><td>1</td>
    </tr><tr>
    <td>ATOM</td><td>5</td><td>N</td><td>N</td><td>.</td><td>ASN</td><td>A</td>
    <td>1</td><td>2</td><td>?</td><td>-7.656</td><td>2.923</td><td>3.155</td>
    <td>1.00</td><td>15.02</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>2</td><td>ASN</td><td>A</td><td>N</td><td>1</td>
    </tr><tr>
    <td>ATOM</td><td>6</td><td>C</td><td>CA</td><td>.</td><td>ASN</td><td>A</td>
    <td>1</td><td>2</td><td>?</td><td>-6.522</td><td>2.038</td><td>2.831</td>
    <td>1.00</td><td>14.10</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>2</td><td>ASN</td><td>A</td><td>CA</td><td>1</td>
    </tr><tr>
    <td>ATOM</td><td>7</td><td>C</td><td>C</td><td>.</td><td>ASN</td><td>A</td>
    <td>1</td><td>2</td><td>?</td><td>-5.241</td><td>2.537</td><td>3.427</td>
    <td>1.00</td><td>13.13</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>2</td><td>ASN</td><td>A</td><td>C</td><td>1</td>
    </tr><tr>
    <td>ATOM</td><td>8</td><td>O</td><td>O</td><td>.</td><td>ASN</td><td>A</td>
    <td>1</td><td>2</td><td>?</td><td>-4.978</td><td>3.742</td><td>3.426</td>
    <td>1.00</td><td>11.91</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>2</td><td>ASN</td><td>A</td><td>O</td><td>1</td>
    </tr><tr>
    <td>ATOM</td><td>9</td><td>C</td><td>CB</td><td>.</td><td>ASN</td><td>A</td>
    <td>1</td><td>2</td><td>?</td><td>-6.346</td><td>1.881</td><td>1.341</td>
    <td>1.00</td><td>15.38</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>2</td><td>ASN</td><td>A</td><td>CB</td><td>1</td>
    </tr><tr>
    <td>ATOM</td><td>10</td><td>C</td><td>CG</td><td>.</td><td>ASN</td><td>A</td>
    <td>1</td><td>2</td><td>?</td><td>-7.584</td><td>1.342</td><td>0.692</td>
    <td>1.00</td><td>14.08</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>2</td><td>ASN</td><td>A</td><td>CG</td><td>1</td>
    </tr><tr>
    <td>ATOM</td><td>11</td><td>O</td><td>OD1</td><td>.</td><td>ASN</td><td>A</td>
    <td>1</td><td>2</td><td>?</td><td>-8.025</td><td>0.227</td><td>1.016</td>
    <td>1.00</td><td>17.46</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>2</td><td>ASN</td><td>A</td><td>OD1</td><td>1</td>
    </tr><tr>
    <td>ATOM</td><td>12</td><td>N</td><td>ND2</td><td>.</td><td>ASN</td><td>A</td>
    <td>1</td><td>2</td><td>?</td><td>-8.204</td><td>2.155</td><td>-0.169</td>
    <td>1.00</td><td>11.72</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>2</td><td>ASN</td><td>A</td><td>ND2</td><td>1</td>
    </tr><tr>
    <td colspan="25">...</td>
    </tr><tr>
    <td>ATOM</td><td>58</td><td>O</td><td>OH</td><td>.</td><td>TYR</td><td>A</td>
    <td>1</td><td>7</td><td>?</td><td>3.766</td><td>0.589</td><td>10.291</td>
    <td>1.00</td><td>14.39</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>7</td><td>TYR</td><td>A</td><td>OH</td><td>1</td>
    </tr><tr>
    <td>ATOM</td><td>59</td><td>O</td><td>OXT</td><td>.</td><td>TYR</td><td>A</td>
    <td>1</td><td>7</td><td>?</td><td>11.358</td><td>2.999</td><td>7.612</td>
    <td>1.00</td><td>17.49</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>7</td><td>TYR</td><td>A</td><td>OXT</td><td>1</td>
    </tr><tr>
    <td>HETATM</td><td>60</td><td>O</td><td>O</td><td>.</td><td>HOH</td><td>B</td>
    <td>2</td><td>.</td><td>?</td><td>-6.471</td><td>5.227</td><td>7.124</td>
    <td>1.00</td><td>22.62</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>8</td><td>HOH</td><td>A</td><td>O</td><td>1</td>
    </tr><tr>
    <td>HETATM</td><td>61</td><td>O</td><td>O</td><td>.</td><td>HOH</td><td>B</td>
    <td>2</td><td>.</td><td>?</td><td>10.431</td><td>1.858</td><td>3.216</td>
    <td>1.00</td><td>19.71</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>9</td><td>HOH</td><td>A</td><td>O</td><td>1</td>
    </tr><tr>
    <td>HETATM</td><td>62</td><td>O</td><td>O</td><td>.</td><td>HOH</td><td>B</td>
    <td>2</td><td>.</td><td>?</td><td>-11.286</td><td>1.756</td><td>-1.468</td>
    <td>1.00</td><td>17.08</td><td>?</td><td>?</td><td>?</td><td>?</td><td>?</td>
    <td>?</td><td>10</td><td>HOH</td><td>A</td><td>O</td><td>1</td>
    </tr><tr>
    </tr>
    </table>

It is hard to manually find an example
where ``_atom_site.auth_atom_id`` differs from ``_atom_site.label_atom_id``, or
where ``_atom_site.auth_comp_id`` differs from ``_atom_site.label_comp_id``.
So here is a small C++ program:

.. literalinclude:: ../examples/auth_label.cc

We compile it, run it, and come back after an hour:

.. code-block:: none

    $ g++-6 -O2 -I.. -Ithird_party examples/auth_label.cc -lstdc++fs -lz
    $ ./a.out pdb_copy/mmCIF
    3D3W: atom_id  O1 -> OD
    1TNI: atom_id  HN2 -> HN3
    1S01: comp_id  CA -> UNL
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    2KI7: atom_id  H -> H1
    4E1U: atom_id  H101 -> H103
    4WH8: atom_id  H3 -> H391
    4WH8: atom_id  H4 -> H401
    3V2I: atom_id  UNK -> UNL
    1AGG: atom_id  H3 -> H
    1AGG: atom_id  H3 -> H
    1AGG: atom_id  H3 -> H

So, as of April 2017, only a single author's residue name was changed,
and atom names were changed in 7 PDB entries.

Amino acid frequency
--------------------

.. highlight:: python

Let say we see in a `paper <https://doi.org/10.1093/nar/gkw978>`_
amino acid frequency averaged over 5000+ proteomes
(Ala 8.76%, Cys 1.38%, Asp 5.49%, etc),
and we want to compare it with the frequency in the PDB database.
So we write a little script

.. literalinclude:: ../examples/aafreq.py

We can run this script on our copy of the PDB database
(mmCIF format) and get such an output:

.. code-block:: none

    200L ALA:17 LEU:15 LYS:13 ARG:13 THR:12 ASN:12 GLY:11 ILE:10 ASP:10 VAL:9 GLU:8 SER:6 TYR:6 GLN:5 PHE:5 MET:5 PRO:3 TRP:3 HIS:1
    ...
    4ZZZ LEU:40 LYS:33 SER:30 ASP:25 GLY:25 ILE:25 VAL:23 ALA:19 PRO:18 GLU:18 ASN:17 THR:16 TYR:16 GLN:14 ARG:11 PHE:10 HIS:9 MET:9 CYS:2 TRP:2
    TOTAL LEU:8.90% ALA:7.80% GLY:7.34% VAL:6.95% GLU:6.55% SER:6.29% LYS:6.18% ASP:5.55% THR:5.54% ILE:5.49% ARG:5.35% PRO:4.66% ASN:4.16% PHE:3.88% GLN:3.77% TYR:3.45% HIS:2.63% MET:2.14% CYS:1.38% TRP:1.35%

On my laptop it takes about an hour, using a single core.
Most of this hour is spent on tokenizing the CIF files and copying
the content into a DOM structure, what could be largely avoided given
that we use only sequences not atoms.
But it is not worth to optimize one-off scripts.
The same goes for using multiple processor cores.

Custom PDB search
-----------------

We may need to go through a local copy of the PDB archive to find entries
according to criteria that cannot be queried in RCSB/PDBe/PDBj web interfaces.
For example, if for some reason we need PDB entries with large number
of anisotropic B-factors, we may write a quick script:

.. literalinclude:: ../examples/simple_search.py

and then run it for an hour or so.

.. code-block:: none

    $ ./examples/simple_search.py pdb_copy/mmCIF
    5AEW 60155
    5AFI 98325
    3AIB 52592
    ...

Search PDB by elements
----------------------

.. image:: img/periodic-table-thumb.png
   :align: right
   :scale: 100
   :target: https://project-gemmi.github.io/periodic-table/

Let say we want to be able to search the PDB by specifying a set of elements
present in the model. First we write down elements present in each
PDB entry::

    block = cif.read(path).sole_block()
    elems = set(block.find_loop("_atom_site.type_symbol"))
    print(name + ' ' + ' '.join(elems))

This example ended up overdone a bit and it was put into a
`separate repository <https://github.com/project-gemmi/periodic-table>`_.

Here is a demo: `<https://project-gemmi.github.io/periodic-table/>`_


Solvent content vs resolution
-----------------------------

.. |Vm| replace:: *V*\ :sub:`M`
.. |Vs| replace:: *V*\ :sub:`S`
.. |dmin| replace:: *d*\ :sub:`min`

Let say that we would like to generate a plot of solvent content
as a function of |dmin|, similar to the plots
by C. X. Weichenberger and B. Rupp
in `Acta Cryst D <https://www.ncbi.nlm.nih.gov/pubmed/24914969>`_
and `CNN <http://www.phenix-online.org/newsletter/CCN_2015_01.pdf#page=14>`_),
but less smoothed and with duplicate entries included.

In macromolecular crystallography solvent content is conventionally
estimated as:

    |Vs| = 1 -- 1.230 / |Vm|

where |Vm| is Matthews coefficient defined as |Vm|\ =\ *V*\ /\ *m*
(volume of the asymmetric unit over the molecular weight of all
molecules in this volume).

Weichenberger and Rupp calculated |Vs| themselves, but we will simply
use the values present in mmCIF files.
So first we extract |Vm|, |Vs| and |dmin|,
as well as number of DNA/RNA chains (protein-only entries will have 0 here),
deposition date and group ID (which will be explained later).

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
They came from two European high-throughput beamlines and
serve as an illustration of how automated software can analyze hundreds
of similar samples (fragment screening) and submit them quickly to the PDB.

We can easily filter out group depositions -- either using the group IDs
extracted from mmCIF files (we do this for pdb-stats_)
or, like W&B, by excluding redundant entries based on the unit cell and |Vm|.

.. _pdb-stats: https://project-gemmi.github.io/pdb-stats/

Not all the dark spots are group depositions.
For example, the one at |Vs|\ ≈66.5%, |dmin| 2.5-3A is proteasome 20S
studied over years by Huber *et al*, with dozens PDB submissions.

Weights
-------

Let say we would like to verify consistency of molecular weights.
First, let us look at chem_comp tables:

.. code-block:: none

    loop_
    _chem_comp.id
    _chem_comp.type
    _chem_comp.mon_nstd_flag
    _chem_comp.name
    _chem_comp.pdbx_synonyms
    _chem_comp.formula
    _chem_comp.formula_weight
    ALA 'L-peptide linking' y ALANINE         ?                 'C3 H7 N O2'     89.093
    ARG 'L-peptide linking' y ARGININE        ?                 'C6 H15 N4 O2 1' 175.209
    ASN 'L-peptide linking' y ASPARAGINE      ?                 'C4 H8 N2 O3'    132.118
    ASP 'L-peptide linking' y 'ASPARTIC ACID' ?                 'C4 H7 N O4'     133.103
    CYS 'L-peptide linking' y CYSTEINE        ?                 'C3 H7 N O2 S'   121.158
    EDO non-polymer         . 1,2-ETHANEDIOL  'ETHYLENE GLYCOL' 'C2 H6 O2'       62.068
    ...

We expect that by using molecular weights of elements and a simple arithmetic
we can recalculate ``_chem_comp.formula_weight`` from ``_chem_comp.formula``.
The full code is in :file:`examples/weights.py``.
It includes a function that converts ``'C2 H6 O2'`` to ``{C:2, H:6, O:2}``.
Here we only show the few lines of code that sum the element weights
and compare the result:

.. literalinclude:: ../examples/weight.py
   :pyobject: check_chem_comp_formula_weight

This script prints differences above 0.1 u:

.. code-block:: none

    3A1R DOD  D2 O                20.028 -    18.015 = +2.013
    5AI2 D8U  D 1                  2.014 -         ? = +nan
    5AI2 DOD  D2 O                20.028 -    18.015 = +2.013
    4AR3 D3O  O 1                 15.999 -    22.042 = -6.043
    4AR4 D3O  O 1                 15.999 -    22.042 = -6.043
    4AR5 DOD  D2 O                20.028 -    18.015 = +2.013
    4AR6 DOD  D2 O                20.028 -    18.015 = +2.013
    4B1A K3G  MO12 O40 P 1      1822.350 -  1822.230 = +0.120
    4BVO E43  H2 O40 W12 6      2848.072 -  2848.192 = -0.120
    ...
    5FHW HFW  Hf O61 P2 W17     4341.681 -  4343.697 = -2.016
    ...

The differences are few and minor.
We see a few PDB entries with the weight of D\ :sub:`2`\ O
set to the weight H\ :sub:`2`\ O. The second line shows missing weight.
The differences +0.12 and -0.12 next to Mo12 and W12 probably come from
the 0.01u difference in the input masses of the elements.
In two entries D3O is missing D in the formula, and -2.016 in HFW
suggests two missing hydrogens.

Now let us try to re-calculate ``_entity.formula_weight`` from the chem_comp
weights and the sequence.
The PDB software calculates it as a sum of components in the chain,
minus the weight of N-1 waters.
In case of nucleic acids also PO\ :sub:`2` is subtracted
(why not PO\ :sub:`3`\ ? -- to be checked).
And in case of microheterogeneity only the main conformer is taken into
account. As the PDB software uses single precision for these computations,
we ignore differences below 0.003%, which we checked to be enough to
account for numerical errors.

.. literalinclude:: ../examples/weight.py
   :pyobject: check_entity_formula_weight

Running this script on a local copy of the PDB database prints 26 lines,
and the difference is always (except for 4PMN) the mass of PO\ :sub:`2`
showing that we have not fully reproduced the rule when to subtract this group.

.. code-block:: none

    ...
    4D67 entity_id: 44 1625855.77 - 1625917.62 =  -61.855
    1HZS entity_id:  1    1733.88 -    1796.85 =  -62.973
    2K4G entity_id:  1    2158.21 -    2221.18 =  -62.974
    3MBS entity_id:  1    2261.35 -    2324.33 =  -62.974
    1NR8 entity_id:  2    2883.07 -    2946.05 =  -62.974
    3OK2 entity_id:  1    3417.14 -    3354.17 =  +62.968
    ...

Disulfide bonds
---------------

If we were curious what residues take part in disulfide bonds we could
write a little script that inspects annotation in the _struct_conn category.
But to show something else, here we will use ``gemmi-grep``, a little
utility that is documented in a :ref:`separate section <grep>`.

First we try how to extract interesting data from a single entry:

.. code-block:: console

    $ ./gemmi-grep --delimiter=' ' _struct_conn.conn_type_id \
    > -a _struct_conn.ptnr1_label_comp_id -a _struct_conn.ptnr1_label_atom_id \
    > -a _struct_conn.ptnr2_label_comp_id -a _struct_conn.ptnr2_label_atom_id \
    > 5CBL
    5CBL disulf CYS SG BME S2
    5CBL covale ILE CD1 BME C2

Next we will use traditional Unix shell utilities to process the text.
``| grep disulf`` limits the output to the disulfide bonds.
Then we add ``| awk '{ print $3, $4 "\n" $5, $6 }'`` to change the first
line of the output above into two lines:

.. code-block:: none

    CYS SG
    BME S2

Then we run it on the whole PDB archive, sort, count and print the
statistics. The complete command is:

.. code-block:: console

    $ ./gemmi-grep --delimiter=' ' _struct_conn.conn_type_id \
    > -a _struct_conn.ptnr1_label_comp_id -a _struct_conn.ptnr1_label_atom_id \
    > -a _struct_conn.ptnr2_label_comp_id -a _struct_conn.ptnr2_label_atom_id \
    > /hdd/mmCIF \
    > | grep disulf | awk '{ print $3, $4 "\n" $5, $6 }' \
    > | sort | uniq -c | sort -nr
       367980 CYS SG
        274 DCY SG
         28 BME S2
         23 SO4 S
         19 CY3 SG
         14 MET SD
         13 MTN S1
         12 V1A S3
         10 MPT SG
          8 NCY SG
          7 LE1 SG
          7 CSX SG
          6 CSO SG
          5 GSH SG2
          5 DTT S1
          4 SC2 SG
          4 MRG S24
          3 H2S S
          3 DTT S4
          3 1WD S7
          2 P8S S1
          2 MTE S2'
          2 MTE S1'
          2 EPE S
          2 DHL SG
          2 DCD S1
          2 CSD SG
          2 COM S1
          2 6LN S2
          2 6LN S1
          2 4K3 S4
          2 3C7 S01
          2 2ON S2
          1 SCN S
          1 RXR SD
          1 RXR S10
          1 MEE S
          1 G47 SG
          1 COA S1P
          1 6ML S1
          1 5O8 SBH

And what other bond types are annotated in ``_struct_conn``?

.. code-block:: console

    $ ./gemmi-grep -O -b _struct_conn.conn_type_id /hdd/mmCIF/ | sort | uniq -c
     501851 covale
     184231 disulf
    4035115 hydrog
    1394732 metalc


mmJSON-like data
----------------

Gemmi has a built-in support for mmJSON and comes with
a :ref:`converter <json>` utility, but just as an excercise
let us convert mmJSON to mmCIF in Python:

.. literalinclude:: ../examples/from_json.py
   :lines: 6-

.. _ccd_example:

Chemical Component Dictionary
-----------------------------

For something a bit different, let us look at the data from
`CCD <https://www.wwpdb.org/data/ccd>`_.
The :file:`components.cif` file describes all the monomers (residues,
ligands, solvent molecules) from the PDB entries.
It contains information about bonds that is absent in the PDB entries.

.. note::

    The absence of this information in mmCIF files is a
    `well-known problem <https://www.cgl.ucsf.edu/chimera/data/mmcif-oct2013/mmcif.html>`_,
    mitigated somewhat by PDBe which in parallel to the wwPDB archive has also
    `mmCIF files with connectivity <https://doi.org/10.1093/nar/gkv1047>`_
    and bond-order information;
    and by RCSB which has this information in the
    `MMTF format <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5473584/#__tag_618683893>`_.

As an exercise, let us check heavy atoms in selenomethionine:

.. doctest::

    >>> from gemmi import cif
    >>> cif.read('components.cif')  #doctest: +ELLIPSIS
    <gemmi.cif.Document object at 0x...>
    >>> _.find_block('MSE')
    <gemmi.cif.Block MSE>
    >>> _.find('_chem_comp_atom.', ['atom_id', 'type_symbol'])
    <gemmi.cif.Table 20 x 2>
    >>> [':'.join(a) for a in _ if a[1] != 'H']
    ['N:N', 'CA:C', 'C:C', 'O:O', 'OXT:O', 'CB:C', 'CG:C', 'SE:SE', 'CE:C']

One may wonder what is the heaviest CCD component:

.. doctest::

    >>> ccd = cif.read('components.cif')
    >>> max(ccd, key=lambda b: float(b.find_value('_chem_comp.formula_weight')))
    <gemmi.cif.Block WO2>
    >>> _.find_value('_chem_comp.formula')
    '"O62 P2 W18"'

The :file:`components.cif` file is big, so we may want to split it into
multiple file. Here we write all components on letter A to a new file.
Additionally, we remove the descriptor category:

.. doctest::

    >>> ccd = cif.read('components.cif')
    >>> len(ccd)  #doctest: +SKIP
    25219
    >>> to_be_deleted = [n for n, block in enumerate(ccd) if block.name[0] != 'A']
    >>> for index in reversed(to_be_deleted):
    ...   del ccd[index]
    ...
    >>> len(ccd)  #doctest: +SKIP
    991
    >>> for block in ccd:
    ...     block.find_mmcif_category('_pdbx_chem_comp_descriptor.').erase()
    ...
    >>> ccd.write_file('A.cif')

The examples here present low-level, generic handling of CCD as a CIF file.
If we would need to look into coordinates, it would be better to use
higher level representation of the chemical component, which is provided
by :ref:`gemmi::ChemComp <chemcomp>`.

Monomer libraries
-----------------

Macromolecular refinement programs need to know more about monomers
than CCD can tell: they need to know how to restrain the structure.
Therefore, they have own dictionaries, such as Refmac dictionary
(a.k.a. monomer library or cif files).
:file:`examples/monomers.py` has a simple consistency checks
between the CCD and Refmac dictionary. This example will be expanded
after the coming revolution in the Refmac monomer library.
