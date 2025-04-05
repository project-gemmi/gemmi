CIF Parser
##########

This section covers working with CIF files on a syntactic level.

Higher-level functions that understand semantics of different types
of CIF files are documented in other sections. These types include:

* :ref:`small molecule <small_molecules>` or inorganic CIF files,
* macromolecular :ref:`PDBx/mmCIF with coordinates <mmcif_format>`,
* PDBx/mmCIF :ref:`with reflection data <sf_mmcif>`,
* :ref:`monomer/ligand files <CCD_etc>` used to supplement
  macromolecular coordinates.

.. _cif_intro:

What are STAR, CIF, DDL, mmCIF?
==================================

(in case someone comes here looking for a serialization format)

STAR is a human-readable data serialization format (think XML or JSON)
that happens to be known and used only in molecular-structure sciences.

CIF (Crystallographic Information File) -- a file format used
in crystallography -- is a restricted derivative of STAR.
It is restricted in features (to make implementation easier)
but also imposes arbitrary limits -- for example, on line length.

DDL is a schema language for STAR/CIF.

All of them (STAR, CIF and DDL) have multiple versions.
We will be more specific in the following sections.

The STAR/CIF syntax is relatively simple but may be confusing at first.
(Note that the initial version of STAR was published by Sydney Hall in 1991 --
before XML and long before JSON and YAML, not to mention TOML.)

.. highlight:: default

Key-value pairs have the form::

    _year 2017  # means year = 2017

Blocks (sections) begin with the `data_` keyword
with a block name glued to it::

    data_tomato   # [tomato]
    _color red    # color = "red"

Importantly, unlike XML/JSON/YAML/TOML, STAR is designed
to concisely represent tabular data. This example from TOML docs::

    # TOML
    points = [ { x = 1, y = 2, z = 3 },
               { x = 7, y = 8, z = 9 },
               { x = 2, y = 4, z = 8 } ]

could be expressed as a so-called *loop* in STAR (keyword `loop_`
followed by column names followed by values)::

    # STAR/CIF
    loop_
    _points.x _points.y _points.z
    1 2 3
    7 8 9
    2 4 8

Typically, long tables (loops) make up most of the CIF content::

    1    N N   . LEU A 11  ? 0.5281 0.5618 0.5305 -0.0327 -0.0621 0.0104
    2    C CA  . LEU A 11  ? 0.5446 0.5722 0.5396 -0.0317 -0.0632 0.0080
    # many, many lines...
    5331 S SD  . MET C 238 ? 2.2952 2.3511 2.3275 -0.0895 0.0372  -0.0230
    5332 C CE  . MET C 238 ? 1.5699 1.6247 1.6108 -0.0907 0.0388  -0.0244

The dot and question mark in the example above are two null types.
In the CIF spec: `?` = *unknown* and `.` = *not applicable*.
In mmCIF files `.` is used for mandatory items, `?` for non-mandatory.

The CIF syntax has a flaw resulting from historical trade-offs:
strings interpretable as numbers do not need to be quoted.
Therefore, while the CIF specification defines a formal grammar for
numbers, we can't recognize them. The type of `5332` is not certain:
the JSON equivalent could be either `5332` or `"5332"`.

Note: "STAR File"
`used to be <https://www.iucr.org/__data/iucr/lists/comcifs-l/msg00710.html>`_
trademarked by IUCr
and `patented <https://patents.google.com/patent/WO1991016682A1>`_,
but not anymore.

The mmCIF format (by mmCIF we mean what is more formally called PDBx/mmCIF)
is the CIF syntax + a huge dictionary (ontology/schema) in DDL2.
The dictionary defines relations between columns in different tables,
which makes it resemble a relational database (it was designed at the
height of popularity of RDBMSs).

**Where are the specs?**

International Tables for Crystallography
`Vol. G (2006) <http://it.iucr.org/Ga/contents/>`_
describes STAR, CIF 1.1, DDL1 and DDL2.
If you don't have access to it -- the IUCr website has specs of
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
(More about it :ref:`later <ddl2>`.)

.. _what_is_parsed:

What is parsed?
===============

The parser supports the CIF 1.1 spec and some extras.

Currently, it is available as:

* C++14 header-only library, and
* Python (3.8+) extension module.

We use it to read:

* mmCIF files (both coordinates and structure factors)
* CIF files from the Crystallography Open Database (COD)
* Chemical Component Dictionary from PDB
* DDL1 and DDL2 dictionaries from IUCr and PDB
* monomer library a.k.a. Refmac dictionary

The parser handles:

* all constructs of CIF 1.1 (including *save frames*),
* the `global_` and `stop_` keywords from STAR -- needed for the Refmac
  monomer library and `mmcif_nmr-star.dic`, respectively.

It could be extended to handle the new features of CIF 2.0
or even the full STAR format, but we don't have a good reason to do this.
The same goes for DDLm/dREL.

The parser does not handle CIF 1.1 conventions, as they are not part
of the syntax:
line wrapping (``eol;\eol``),
Greek letters (``\m`` → µ),
accented letters (``\'o`` → ó),
special alphabetic characters (``\%A`` → Å)
and other codes (``\\infty`` → ∞).

CIF parsers in the small-molecules field need to deal with incorrect syntax.
The papers about `iotbx.cif <https://doi.org/10.1107/S0021889811041161>`_
and `COD::CIF::Parser <http://dx.doi.org/10.1107/S1600576715022396>`_
enumerate 12 and 10 common error categories, respectively.
Nowadays the problem is less severe, especially in the MX community
that embraced the CIF format later. So we've decided to add only
the following rules to relax the syntax:

* names and lines can have any length like in STAR
  (the CIF spec imposes a limit of 2048 characters, but some mmCIF files
  from PDB exceed it, e.g. 3j3q.cif),
* quoted strings may contain non-ASCII characters (if nothing has changed
  one entry in the PDB has byte A0 corresponding to non-breaking space),
* a table (loop) can have no values if followed by a keyword or EOF
  (such files were written by old versions of Refmac and SHELXL,
  and were also present in the CCP4 monomer library),
* block name (*blockcode*) can be empty, i.e. the block can start
  with bare `data_` keyword (RELION and Buccaneer write such files),
* unquoted strings cannot start with keywords (STAR spec is ambiguous
  about this -- see
  `StarTools doc <http://www.globalphasing.com/startools/>`_ for details;
  this rule is actually about simplifying, not relaxing),
* missing value in a key-value pair is optionally allowed
  if whitespace after the tag ends with a new-line character.
  More specifically, the parsing step allows for missing values in such cases,
  but the next validation step (which can be skipped -- more about it
  :ref:`later <cif_low_level>`) throws an error.


DOM and SAX
===========

The terms DOM and SAX originate from XML parsing, but they also became
names of general parsing styles.
Gemmi can parse CIF files in two ways that correspond to DOM and SAX:

* The usual way is to parse a file or a string into a document (DOM)
  that can be easily accessed and manipulated.

* Alternatively, from C++ only, one can define own
  `PEGTL Actions <https://github.com/taocpp/PEGTL/blob/master/doc/Actions-and-States.md>`_
  for the grammar rules from `cif.hpp`.
  These actions will be triggered while reading a CIF file.

This documentation covers the DOM parsing only.
The hierarchy in the DOM reflects the structure of CIF 1.1:

* Document contains blocks.
* Block can contain name-value pairs, loops and frames.
* Frame can contain name-value pairs and loops.
* Loop (*m*\ ×\ *n* table) contains *n* column names and *m*\ ×\ *n* values.

Names are often called *tags*. The leading `_` is usually treated
as part of the tag, not just a syntactic feature. So we store tag strings
with the underscore (`_my_tag`), and functions that take tags as arguments
expect strings starting with `_`.
The case of tags is preserved.

At this point it might help to show an example of DOM parsing.
Here is a simple program that reads an mmCIF file and prints
weights of components that are listed as `_chem_comp.formula_weight`.
The functions used here are documented later on.

.. tab:: C++

 .. literalinclude:: code/cif_cc.cpp
    :language: cpp
    :lines: 1-16

.. tab:: Python

 .. literalinclude:: ../examples/read_cif.py
    :language: python
    :lines: 2-


Values have types. CIF 1.1 defines four base types:

* char (string)
* uchar (case-insensitive string)
* numb (number that cannot be recognized as number on the syntax level)
* null (one of two possible nulls: `?` and `.`)

Additionally, DDL2 dictionaries specify subtypes.
For example, `int` and `float` are mmCIF subtypes of `numb`.

Since in general it is not possible to infer type without
the corresponding dictionary, the DOM stores raw strings (including quotes).
They can later be converted to the required type using the following helper
functions:

* `as_string()` -- gets unquoted string
* `as_number()` -- gets floating-point number
* `as_int()` -- gets integer
* `as_char()` -- gets single character
* `is_null()` -- check if the value is null (i.e. `?` or `.`)

and we also have:

* `quote()` -- the opposite of `as_string()` -- adds quotes appropriate
  for the content of the string (usually, quotes are not necessary and
  not added).


.. highlight:: cpp

.. tab:: C++

 ::

  #include <gemmi/cifdoc.hpp> // for as_string, as_int, quote
  #include <gemmi/numb.hpp>   // for as_number

  double rfree = cif::as_number(raw_rfree_string); // NaN if it's '?' or '.'

  cif::as_int("123");      // 123
  cif::quote("two words"); // "'two words'"
  cif::is_null("?");       // true

.. tab:: Python

 .. doctest::

  >>> from gemmi import cif
  >>> cif.as_number('123')
  123.0
  >>> cif.as_int('123')
  123
  >>> cif.quote('two words')
  "'two words'"
  >>> cif.as_string(_)
  'two words'
  >>> cif.as_number(_)
  nan


Reading a file
==============

We have a few reading functions that read a file (or a string, or a stream)
and return a document (DOM) -- an instance of class `Document`.

.. tab:: C++

 .. tab:: libgemmi_cpp

  ::

   #include <gemmi/read_cif.hpp>

   // functions in namespace gemmi, from -lgemmi_cpp, usually linked with zlib or zlib-ng

   // similar to cif::read_file, but uncompresses *.gz files on the fly
   cif::Document read_cif_gz(const std::string& path);

   // reads the content of CIF file from a memory buffer (name is used when reporting errors)
   cif::Document read_cif_from_memory(const char* data, size_t size, const char* name);

   // reads from a string; it's in namespace gemmi::cif (for backward compatibility)
   Document read_string(const std::string& data);

 .. tab:: header-only

  ::

   #include <gemmi/cif.hpp>

   // Header-only functions in namespace gemmi::cif.
   // No linking, but slower compilation.
   // Parameter name is used only when reporting errors.

   Document read_file(const std::string& filename);
   Document read_memory(const char* data, const size_t size, const char* name);
   // Parameter bufsize determines the buffer size and only affects performance.
   // These functions are slower than the ones above.
   Document read_cstream(std::FILE *f, size_t bufsize, const char* name);
   Document read_istream(std::istream &is, size_t bufsize, const char* name);

.. tab:: Python

 .. testcode::

  from gemmi import cif

  # read and parse a CIF file; if the filename ends with .gz it is uncompressed on the fly
  doc = cif.read_file('components.cif')

  # the same (except that it can also read mmJSON)
  doc = cif.read('components.cif')

  # read content of a CIF file from string or bytes
  doc = cif.read_string('data_this _is valid _cif content')

The functions that can read gzipped files also understand `path` specified
as `-` to mean the standard input.

.. _cif_low_level:

Low-level functions
-------------------

The `read_file()` call is equivalent to the following sequence:

.. testcode::

  path = 'components.cif'
  doc = cif.Document()
  doc.source = path
  doc.parse_file(path)
  doc.check_for_missing_values()
  doc.check_for_duplicates()

The last two functions check for, respectively, missing values in tag-value
pairs and duplicated names.
It is possible to read erroneous CIF files by skipping these checks.

Analogically, function `read_string()` can be replaced by a similar
sequence with `Document.parse_string()` doing the main job:

.. doctest::

  >>> doc = cif.Document()
  >>> doc.parse_string('data_this _tag_has_no_value\n')
  >>> doc[0].find_value('_tag_has_no_value')
  ''
  >>> try: doc.check_for_missing_values()
  ... except RuntimeError: print('missing value')
  ...
  missing value


Writing a file
==============

Document (or a single Block) can be written to a file or to a string.
Reading from a file and then writing to a file does not preserve whitespace.
The formatting of the output file is controlled by `cif::WriteOptions`,
which contains the following fields (by default, they are all false or 0):

* `prefer_pairs` (bool) -- if set to true, write single-row loops as pairs,
* `compact` (bool) -- if set to true, do not add blank lines between categories,
* `misuse_hash` (bool) -- if set to true, put `#` (empty comments) between
  categories -- the peculiar formatting used in the wwPDB archive;
  enables diff-ing (`diff --ignore-space-change`) with other such files,
* `align_pairs` (int) -- pad tags in tag-value pairs to this width
  (if set to 33, the values will be aligned at column 35,
  except where tags have more than 33 characters),
* `align_loops` (int) -- if non-zero, columns in loops are aligned
  to the maximum string width in each column, but not more than
  this value; if one string in the column is very wide, that row will be
  misaligned, which is usually preferable to excessive padding.

C++
---

The functions writing `cif::Document` and `cif::Block` to C++ stream
is in a separate header `gemmi/to_cif.hpp`::

  void write_cif_to_stream(std::ostream& os, const Document& doc, WriteOptions options)
  void write_cif_block_to_stream(std::ostream& os, const Block& block, WriteOptions options)

Python
------

In Python, the function that writes the document to a file is a method
of the `Document` class:

.. doctest::

  >>> doc.write_file('1pfe-modified.cif')

The output can be styled by passing `WriteOptions` as a second argument:

.. doctest::

  >>> options = cif.WriteOptions()
  >>> options.align_pairs = 33
  >>> options.align_loops = 30
  >>> doc.write_file('1pfe-aligned.cif', options)

The `Document` class also has a method `as_string()` which returns
the text that would be written by `write_file()`.

The `Block` class also has methods `write_file()` and `as_string()`
with the same arguments.

Document
========

`Document` is made of blocks with data. The blocks can be iterated over,
accessed by index or by name (each CIF block must have a unique name).

As it is common for cif files to contain only a single block,
gemmi has a method `sole_block()` that returns the first block
if the document has only one block; otherwise it throws an exception.

At last, is also has a member variable `source` that contains
the path of the file from which the document was read (if it was read
from a file).

C++
---

`Document` has the two member variables::

  std::string source;  // filename or the name passed to read_memory()
  std::vector<Block> blocks;

Each `Block` corresponds to a data block.
To access a block with known name use::

  Block* Document::find_block(const std::string& name)

To access the only block in the file you may use::

  Block& Document::sole_block()

A new `Document` instance can be created with default constructor.
To modify a document you need to access directly its member variables.
With one exception: when adding new blocks you can use a function that
additionally checks if the new name is unique::

  Block& Document::add_new_block(const std::string& name, int pos=-1)


Python
------

`Document` can be iterated, accessed by block index and by block name:

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

It has two other ways of accessing a block:

.. doctest::

  >>> # The function block.find_block(name) is like block[name] ...
  >>> doc.find_block('MSE')
  <gemmi.cif.Block MSE>
  >>> # ... except when the block is not found:
  >>> doc.find_block('no such thing')  # -> None
  >>> # doc['no such thing'] # -> KeyError

  >>> # Get the only block; throws exception if the document has more blocks.
  >>> cif.read("../tests/1pfe.cif.gz").sole_block()
  <gemmi.cif.Block 1PFE>

Blocks can be inserted (by default -- appended after existing blocks)
using one of the two functions:

* `Document.add_new_block(name, pos=-1)`
* `Document.add_copied_block(block, pos=-1)`

As an example, here is how to start a new document:

.. doctest::

  >>> d = cif.Document()
  >>> block_one = d.add_new_block('block-one')
  >>> # populate block_one

To delete a block use `Document.__delitem__` (for example: `del doc[1]`).

Document has also one property

.. doctest::

  >>> doc.source
  'components.cif'


.. warning::

   Adding and removing blocks may invalidate references to other blocks
   in the same Document. This is expected when working with a C++ vector,
   but when using Gemmi from Python it is a flaw.
   The same applies to functions that add/remove items in a block.
   More precisely:

   * functions that add items (e.g. `add_new_block`) may cause memory
     re-allocation invalidating references to all other items (blocks),
   * functions that remove items (`__delitem__`) invalidate references to
     all items after the removed one.

   This means that you need to update a reference before using it:

    .. code-block:: python

       block = doc[0]
       st.add_new_block(...)     # block gets invalidated
       block = st[0]             # block is valid again

.. _cif_block:

Block
=====

Each block has a name and a list of items.
Each item is one of:

* name-value pair (Pair),
* table, a.k.a loop (Loop)
* or save frame (Block -- the same data structure as for block).

A block headed by the word `global_`, part of the STAR syntax, although
`not allowed <https://onlinelibrary.wiley.com/iucr/itc/Ga/ch2o2v0001/sec2o2o7o1o9o5/>`_
in CIF, is :ref:`parsed <what_is_parsed>` into a Block with empty name.

A block headed by bare `data_`, although not allowed neither in CIF nor
in STAR, is parsed into a Block with name set to " " (the space character).

C++
---

Each block contains::

    std::string name;
    std::vector<Item> items;

where `Item` is implemented as an unrestricted union
that holds one of Pair, Loop or Block.

Python
------

Each block has a name:

.. doctest::

  >>> doc = cif.read("../tests/1pfe.cif.gz")
  >>> block = doc.sole_block()
  >>> block.name
  '1PFE'

and a list of items (class Item):

.. doctest::

  >>> for item in block:
  ...    if item.line_number > 1670:
  ...        if item.pair is not None:
  ...            print('pair', item.pair)
  ...        elif item.loop is not None:
  ...            print('loop', item.loop)
  ...        elif item.frame is not None:
  ...            print('frame', item.frame)
  ...
  pair ('_ndb_struct_conf_na.entry_id', '1PFE')
  pair ('_ndb_struct_conf_na.feature', "'double helix'")
  loop <gemmi.cif.Loop 8 x 25>
  loop <gemmi.cif.Loop 5 x 43>
  loop <gemmi.cif.Loop 3 x 3>
  loop <gemmi.cif.Loop 83 x 10>


Frame
=====

(Very few people need it, skip this section.)

The *named save frames* (keyword `save_`) from the STAR specification
are used in CIF files only as sub-sections of a block.
The only place where they are encountered are mmCIF dictionaries.

The save-frame is stored as `Block` and can be accessed with::

  Block* Block::find_frame(std::string name)

.. doctest::

  >>> frame = block.find_frame('my_frame')

Additionally, in C++ one may iterate over all Block's items,
check each item type and handle all the save frames::

  for (cif::Item& item : block.items)
    if (item.type == cif::ItemType::Frame)
      // doing something with item.frame which is a (nested) Block
      cif::Block& frame = item.frame;


Pairs and Loops
===============

The functions in this section can be considered low-level, because they are
specific to either name-value pairs or to loops.


.. warning::

    Assuming that a certain value is contained in a pair, or that it is
    in a loop, is a common source of bugs when working with mmCIF files.
    To get the data, instead of using the functions in this section,
    it is recommended to use the abstractions (introduced in the next sections)
    that work regardless of whether the values are in pairs or loops.

    For instance, when working with proteins one could assume that anisotropic
    ADP values are in a loop, but wwPDB also has entries with anisotropic ADP
    for only one atom -- as name-value pairs.
    Conversely, one could think that R-free is always given as
    a name-value pair, but in entries from the joint X-ray and neutron
    refinement it is in a loop.
    Moreover, some CIF representations, such as mmJSON, do not preserve
    the distinction between pairs and loops at all.

The next sections introduce function that work with both pairs and loops.

C++
---

Pair is simply defined as::

    using Pair = std::array<std::string, 2>;

A pair with a particular tag can be located using::

  const Pair* Block::find_pair(const std::string& tag) const

or, if you want just the value::

  const std::string* Block::find_value(const std::string& tag) const

Both functions return `nullptr` if the tag is not found
(but the latter also searches inside CIF loops and if there
is a matching tag with only a single value, that value is returned).

To add a pair to the block, or modify an existing one, use::

  void Block::set_pair(const std::string& tag, std::string value)

If the value needs quoting, the passed argument needs to be already quoted
(you may pass `cif::quote(value)`).

----

Loop is defined as::

  struct Loop {
    std::vector<std::string> tags;
    std::vector<std::string> values;
    // and a number of functions
  };

To get values corresponding to a tag in a loop (table) you may use::

  Column Block::find_loop(const std::string& tag)

`struct Column`, which is documented further on, has method `get_loop()`
which gives access to `struct Loop`.

A new loop can be added using function::

  Loop& Block::init_loop(const std::string& prefix, std::vector<std::string> tags)

Then it can be populated by either setting directly `tags` and `values`,
or by using Loop's methods such as `add_row()` or `set_all_values()`.

Loop has a few other methods for editing its content, such as
`move_row()`, `add_columns()` and `remove_column()`.

Python
------

Accessing name-value pairs:

.. doctest::

  >>> # (1) tag and value
  >>> block.find_pair('_cell.length_a')
  ('_cell.length_a', '39.374')
  >>> block.find_pair('_no_such_tag')  # return None

  >>> # (2) only value
  >>> block.find_value('_cell.length_b')
  '39.374'
  >>> block.find_value('_cell.no_such_tag')  # returns None

  >>> # (3) Item
  >>> item = block.find_pair_item('_cell.length_c')
  >>> item.pair
  ('_cell.length_c', '79.734')
  >>> item.line_number
  72
  >>> block.find_pair_item('_nothing')  # return None

To add a name-value pair, replacing current item if it exists,
use function `set_pair`:

.. doctest::

  >>> block.set_pair('_year', '2030')

If a new item is added, it is placed at the end of the block.
Then you can move it to a more appropriate position with `move_item()`.
Alternatively, you can use function `set_pairs` that takes a prefix and,
when adding a new item, places it after the last item with the given prefix:

.. doctest::

  >>> block.set_pairs('_cell.', {'length_a_esd': '?',
  ...                            'length_b_esd': '?',
  ...                            'length_c_esd': '?'}, raw=True)

(In C++ the equivalent is `cif::ItemSpan(block, "_cell.").set_pair(…)`).
This is recommended when editing mmCIF files, because all name-value pairs
in the same category must be consecutive (an unwritten rule of the PDB).
The argument `raw` has the same meaning as in the function
:ref:`set_mmcif_category <set_mmcif_category>`.

If the value needs quoting, it must be passed quoted:

.. doctest::

  >>> block.set_pair('_title', cif.quote('Goldilocks and the Three Bears'))

Now we can create a CIF file can from scratch:

.. doctest::

  >>> d = cif.Document()
  >>> d.add_new_block('oak')
  <gemmi.cif.Block oak>
  >>> _.set_pair('_nut', 'acorn')
  >>> print(d.as_string().strip())
  data_oak
  _nut acorn

----

To access values in loop:

.. doctest::

  >>> # (1) get a Column in Loop
  >>> block.find_loop('_atom_type.symbol')
  <gemmi.cif.Column _atom_type.symbol length 6>
  >>> list(_)
  ['C', 'CL', 'N', 'O', 'P', 'S']

  >>> # (2) get Item containing the Loop
  >>> block.find_loop_item('_atom_type.symbol')  # doctest: +ELLIPSIS
  <gemmi.cif.Item object at 0x...>
  >>> _.loop
  <gemmi.cif.Loop 6 x 1>

To add a row to an existing table (loop) use `add_row`:

.. doctest::

  >>> loop = block.find_loop('_atom_type.symbol').get_loop()
  >>> loop.add_row(['Au'], pos=0)
  >>> loop.add_row(['Zr'])  # appends
  >>> list(block.find_loop('_atom_type.symbol'))
  ['Au', 'C', 'CL', 'N', 'O', 'P', 'S', 'Zr']

`add_row` takes as an argument a list of strings, which should be quoted
if necessary. If you have a list Python values use `quote_list` first:

.. doctest::

  >>> cif.quote_list([None, False, 3, -2.5, 'word', 'two words'])
  ['?', '.', '3', '-2.5', 'word', "'two words'"]

Columns can be added and removed:

.. doctest::

  >>> loop.add_columns(['_atom_type.description', '_atom_type.oxidation_number'], value='?')
  >>> loop
  <gemmi.cif.Loop 8 x 3>
  >>> loop.remove_column('_atom_type.description')
  >>> loop
  <gemmi.cif.Loop 8 x 2>

`set_all_values` sets all the data in a table. It takes as an argument
a list of lists of string. The lists of strings correspond to columns.

.. doctest::

  >>> loop = block.find_loop_item('_citation_author.citation_id').loop
  >>> loop.tags
  ['_citation_author.citation_id', '_citation_author.name', '_citation_author.ordinal']
  >>> loop.set_all_values([['primary']*2, [cif.quote('Alice A.'), cif.quote('Bob B.')], ['1', '2']])
  >>> for row in block.find_mmcif_category('_citation_author.'):
  ...   print(list(row))
  ['primary', "'Alice A.'", '1']
  ['primary', "'Bob B.'", '2']

Individual values can be accessed with (row,column) tuples:

.. doctest::

  >>> loop[0,1]
  "'Alice A.'"
  >>> loop[0,1] = "'Carol C.'"
  >>> loop.values  # Loop.values is a read-only list (a copy of all values)
  ['primary', "'Carol C.'", '1', 'primary', "'Bob B.'", '2']

To add a new loop (replacing old one if it exists) use `init_loop`:

.. doctest::

  >>> loop = block.init_loop('_ocean_', ['id', 'name'])
  >>> # empty table is invalid in CIF, we need to add something
  >>> loop.add_row(['1', cif.quote('Atlantic Ocean')])

In the above example, if the block already has tags `_ocean_id`
and/or `_ocean_name` and

* if they are in a table: the table will be cleared and re-used,
* if they are in name-value pairs: the pairs will be removed
  and a table will be created at the position of the first pair.

----

To reorder items, use `block.move_item(old_pos, new_pos)`.
The current position of a tag can be obtained using `block.get_index(tag)`.

.. doctest::

  >>> block.get_index('_entry.id')
  0
  >>> block.get_index('_ocean_id')
  385
  >>> block.move_item(0, -1)  # move first item to the end
  >>> block.get_index('_entry.id')
  385
  >>> block.get_index('_ocean_id')
  384
  >>> block.move_item(385, 0)  # let's move it back (385 == -1 here)

To copy an item:

.. doctest::

  >>> item = block.find_loop_item('_ocean_id')
  >>> new_block = cif.Block('new_block')
  >>> new_block.add_item(item)  # doctest: +ELLIPSIS
  <gemmi.cif.Item object at 0x...>


Column
======

`Column` is a lightweight proxy class for working with both loop columns
and name-value pairs.

It was returned from `find_loop` above, but it is also returned from
a more general function `find_values()`, which searches for a given
tag in both loops and pairs.

C++
---
The C++ signature of `find_values` is::

  Column Block::find_values(const std::string& tag)

`Column` has a few member functions::

  // Number of rows in the loop. 0 means that the tag was not found.
  int length() const;

  // Returns pointer to the column name in the DOM.
  // The name is the same as argument to find_loop() or find_values().
  std::string* get_tag();

  // Returns underlying Item (which contains either Pair or Loop).
  Item* item();

  // If it's in Loop, returns pointer to Loop; otherwise, nullptr.
  Loop* get_loop() const;

  // Get raw value (no bounds checking).
  std::string& operator[](int n);

  // Get raw value (after bounds checking).
  std::string& at(int n);

  // Short-cut for cif::as_string(column.at(n)).
  std::string str(int n) const;

  // Erases item for name-value pair; removes column for Loop
  void erase();

`Column` also provides support for range-based `for`::

  // mmCIF _chem_comp_atom is usually a table, but not always
  for (const std::string &s : block.find_values("_chem_comp_atom.type_symbol"))
    std::cout << s << std::endl;

If the column represents a name-value pair,
`Column::get_loop()` return `nullptr`
(and `Column::get_length()` returns 1).

Python
------

.. doctest::

  >>> block.find_values('_cell.length_a')  # name-value pair
  <gemmi.cif.Column _cell.length_a length 1>
  >>> block.find_values('_atom_site.Cartn_x')  # column in a loop
  <gemmi.cif.Column _atom_site.Cartn_x length 342>

Column's special method `__bool__` tells if the tag was found.
`__len__` returns the number of corresponding values.
`__iter__`,  `__getitem__` and `__setitem__`
get or set a raw string (i.e. string with quotes, if applicable).

As an example, let us shift all the atoms by 1A in the x direction:

.. doctest::

  >>> col_x = block.find_values('_atom_site.Cartn_x')
  >>> for n, x in enumerate(col_x):
  ...   col_x[n] = str(float(x) + 1)

To get the actual string content (without quotes) one may use
the method `str`:

.. doctest::

  >>> column = block.find_values('_chem_comp.formula')
  >>> column[7]
  "'H2 O'"
  >>> column.str(7)  # short-cut for cif.as_string(column[7])
  'H2 O'

The tag is accessible (read/write) through the `tag` property:

.. doctest::

  >>> column.tag
  '_chem_comp.formula'

If the tag was found in a loop, method `get_loop` returns a reference
to this `Loop` in the DOM. Otherwise it returns `None`.

.. doctest::

  >>> column.get_loop()
  <gemmi.cif.Loop 12 x 7>

`erase()` removes the column from a loop, if this column is in a loop;
if it's a tag-value pair it erases the containing item.

.. doctest::

  >>> loop = column.get_loop()
  >>> column.erase()
  >>> loop  # one column less
  <gemmi.cif.Loop 12 x 6>

Table
=====

Usually we want to access multiple columns at once,
so the library has another abstraction (`Table`)
that can be used with multiple tags.

`Table` is returned by `Block.find()`.
Like column, it is a lightweight, iterable view of the data,
but it is for querying multiple related tags at the same time.

The first form of `find()` takes a list of tags::

  Table Block::find(const std::vector<std::string>& tags)

.. doctest::

  >>> block.find(['_entity_poly_seq.entity_id', '_entity_poly_seq.num', '_entity_poly_seq.mon_id'])
  <gemmi.cif.Table 18 x 3>

Since tags in one loop tend to have a common prefix (category name),
the library provides also a second form that takes the common prefix
as the first argument::

  Table Block::find(const std::string& prefix, const std::vector<std::string>& tags)

  // These two calls are equivalent:

  block.find({"_entity_poly_seq.entity_id", "_entity_poly_seq.num", "_entity_poly_seq.mon_id"});
  block.find("_entity_poly_seq.", {"entity_id", "num", "mon_id"});

.. doctest::

  >>> block.find('_entity_poly_seq.', ['entity_id', 'num', 'mon_id'])
  <gemmi.cif.Table 18 x 3>


Note that `find` is not aware of dictionaries and categories,
therefore the category name should end with a separator
(dot for mmCIF files, as shown above).


In the example above, all the tags are required. If one of them is absent,
the returned Table is empty.

.. doctest::

  >>> block.find('_entity_poly_seq.', ['entity_id', 'num', 'non_existent'])
  <gemmi.cif.Table nil>

Tags (all except the first one) can be marked
as *optional* by adding prefix `?`::

  Table table = block.find({"_required_tag", "?_optional_tag"})

.. doctest::

  >>> table = block.find(['_required_tag', '?_optional_tag'])

In such case the returned table may contain either one or two columns.
Before accessing column corresponding to an optional tag one must check
if the column exists with `Table::has_column()` (or, alternatively,
with equivalent function `Table::Row::has()` which will be introduced
later)::

  bool Table::has_column(int n) const

.. doctest::

  >>> block.find('_entity_poly_seq.', ['entity_id', '?num', '?bleh'])
  <gemmi.cif.Table 18 x 3>
  >>> _.has_column(0), _.has_column(1), _.has_column(2)
  (True, True, False)

The `Table` has functions to check its shape::

  bool ok() const;  // true if the table is not empty
  size_t width() const;  // number of columns
  size_t length() const;  // number of rows

.. doctest::

  >>> table = block.find('_entity_poly_seq.', ['entity_id', 'num', 'mon_id'])
  >>> # instead of ok() in Python we use __bool__()
  >>> assert table, "table.__bool__() is expected to return True"
  >>> table.width()  # number of columns
  3
  >>> len(table)  # number of rows
  18

If Table's data is in Loop, the Loop class can be accessed using::

  Loop* get_loop();  // nullptr for tag-value pairs

.. doctest::

  >>> table.loop     # None for tag-value pairs
  <gemmi.cif.Loop 18 x 4>

Tag-value pairs can be converted to a loop using:

.. doctest::

  >>> table.ensure_loop()

but make sure to call it only when the table represents the whole category
(see example :ref:`below <append_row>`).

If a prefix was specified when calling find, the prefix length is stored
and the prefix can be retrieved::

  size_t prefix_length;
  std::string get_prefix() const;

.. doctest::

  >>> table.prefix_length
  17
  >>> table.get_prefix()
  '_entity_poly_seq.'

Table also has function `erase()` that deletes all tags and data
associated with the table.
See an example in the :ref:`CCD section <ccd_example>` below.


Row-wise access
---------------

Most importantly, the Table provides access to data in rows and columns.
We can get a row (`Table::Row`) that in turn provides access to value strings::

  Row Table::operator[](int n)  // access Row
  Row Table::at(int n)          // the same but with bounds checking
  // and also begin(), end() and iterator that enable range-for.

  // Returns the first row that has the specified string in the first column.
  Row Table::find_row(const std::string& s)

  // Makes sure that the table has only one row and returns it.
  Row Table::one()

.. doctest::

  >>> table[0]
  <gemmi.cif.Table.Row: 1 1 DG>

  >>> # and of course it's iterable
  >>> for row in table: pass

  >>> # Returns the first row that has the specified string in the first column.
  >>> table.find_row('2')
  <gemmi.cif.Table.Row: 2 1 DSN>

as well as to the tags::

  Row tags();  // pseudo-row that contains tags

.. doctest::

  >>> table.tags
  <gemmi.cif.Table.Row: _entity_poly_seq.entity_id _entity_poly_seq.num _entity_poly_seq.mon_id>

Tags in the table (in the C++ layer) are references to tags in the CIF block.
If the requested tags are not found, the returned *nil* table has no tags.
Such table has 0 rows, but it can be iterated like an empty list:

.. doctest::

  >>> block.find(['_absent_tag']).tags
  <gemmi.cif.Table.Row:>
  >>> for row in block.find(['_absent_tag']): print(row)
  ... # nothing gets printed, but there is no error


`Table::Row` has functions for accessing the values::

  // Get raw value.
  std::string& Table::Row::operator[](int n)  // no bounds checking
  std::string& Table::Row::at(int n)          // with bounds checking
  // and also begin(), end(), iterator, const_supports iterators.

.. doctest::

  >>> for row in table: print(row[-1], end=',')
  DG,DC,DG,DT,DA,DC,DG,DC,DSN,ALA,N2C,NCY,MVA,DSN,ALA,NCY,N2C,MVA,
  >>>
  >>> row = table[9]
  >>> for value in row: print(value, end=',')
  2,2,ALA,
  >>> row[2]
  'ALA'
  >>> row[-1]  # the same, because this table has two rows
  'ALA'
  >>> row.get(2) # the same, but would return None instead of IndexError
  'ALA'
  >>> row['mon_id']  # the same, but slightly slower than numeric index
  'ALA'
  >>> row['_entity_poly_seq.mon_id']  # the same
  'ALA'

and a few convenience functions, including::

  size_t Table::Row::size() const           // the width of the table
  std::string Table::Row::str(int n) const  // short-cut for cif::as_string(row.at(n))
  bool Table::Row::has(int n) const         // the same as Table::has_column(n)

.. doctest::

  >>> len(row)    # the same as table.width()
  3
  >>> row.str(2)  # if the value is in quotes, it gets un-quoted
  'ALA'
  >>> row.has(2)  # the same as table.has_column(2)
  True

and a property::

  int row_index

.. doctest::

  >>> row.row_index
  9

Individual values in a row can be directly modified.
This includes tags, which are accessible as a special row.
As an example, let us swap two names
(these two tend to have identical values, so no one will notice):

.. literalinclude:: code/cif_cc.cpp
   :language: cpp
   :lines: 37-39

.. doctest::

  >>> tags = block.find('_atom_site.', ['label_atom_id', 'auth_atom_id']).tags
  >>> tags[0], tags[1] = tags[1], tags[0]

----

.. _append_row:

Function `Table::append_row` appends a row.
It won't work if the table is composed of tag-value pairs (which isn't
the case here, but let's make this example as general as possible),
so we call `Table::ensure_loop()` first.
This function, in turn, won't work well if the table contains only a subset
of the category, thus `find_mmcif_category()`:

.. doctest::

  >>> block.find_mmcif_category('_entity_poly_seq.').ensure_loop()
  >>> table.append_row(['3', '4', 'new'])
  >>> table[-1]
  <gemmi.cif.Table.Row: 3 4 new>
  >>> _.row_index
  18

If new tags are added to this category in the future,
the corresponding values in appended rows will be filled with `.`.

We can also move a row to a different position:

.. doctest::

  >>> table.move_row(old_pos=-1, new_pos=0)  # move the last row to the front
  >>> table[0]
  <gemmi.cif.Table.Row: 3 4 new>

and remove a row from Table (function `Table::remove_row`):

.. doctest::

  >>> table.remove_row(18)  # or: del table[18]

To efficiently remove multiple consecutive rows use
C++ function `Table::remove_rows` or Python `__delitem__` with slice:

.. doctest::

  >>> del table[12:15]

As is usual with any containers (in both Python and C++),
if you want to remove items while iterating over them,
the iteration must be done backward.
Here is an example that removes atoms with zero occupancy:

.. doctest::

  >>> doc = cif.read('../tests/3dg1_final.cif')
  >>> atom_table = doc[0].find(['_atom_site.occupancy'])
  >>> for i in range(len(atom_table)-1, -1, -1):
  ...   if float(atom_table[i][0]) == 0:
  ...     del atom_table[i]
  ...
  >>> doc.write_file('out.cif')

Column-wise access
------------------

`Table` gives also access to columns, represented by the previously
introduced `Column`::

    Column Table::column(int index)

    // alternatively, specify tag name
    Column Table::find_column(const std::string& tag)

.. doctest::

  >>> table.column(0)
  <gemmi.cif.Column _entity_poly_seq.entity_id length 15>
  >>> table.find_column('_entity_poly_seq.mon_id')
  <gemmi.cif.Column _entity_poly_seq.mon_id length 15>

If the table is created in a function that uses prefix,
the prefix can be omitted in `find_column`::

  Table t = block.find("_entity_poly_seq.", {"entity_id", "num", "mon_id"});
  Column col = t.find_column(2);
  // is equivalent to
  Column col = t.find_column("_entity_poly_seq.mon_id");
  // is equivalent to
  Column col = t.find_column("mon_id");

.. doctest::

  >>> table.find_column('mon_id')
  <gemmi.cif.Column _entity_poly_seq.mon_id length 15>

C++ note:
both `Column` and `Table::Row` have functions `begin()` and `end()`
in const and non-const variants, returning `iterator` and
`const_iterator` types, respectively. These types satisfy requirements
of the BidirectionalIterator concept.
Conversely, the iterator over the rows of `Table` is a minimalistic
structure -- just enough get the range-for work.


mmCIF categories
----------------

mmCIF files group data into categories. All mmCIF tags have a dot
(e.g. `_entry.id`) and the category name is the part before the dot.

C++
~~~

We have two functions to work with categories.
One returns a list of all categories in the block::

  std::vector<std::string> Block::get_mmcif_category_names() const

The other returns a `Table` with all tags (and values) belonging to
the specified category::

  Table Block::find_mmcif_category(std::string cat)


Python
~~~~~~

Python bindings have the same two functions:

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

Additionally, two Python-specific functions: `get_mmcif_category`
and `set_mmcif_category` translate between an mmCIF category and
Python dictionary:

.. doctest::

  >>> block.get_mmcif_category('_entry.')
  {'id': ['1PFE']}
  >>> sorted(block.get_mmcif_category('_database_2.').keys())
  ['database_code', 'database_id']

`?` and `.` are translated to `None` and `False`.
To disable this translation and get "raw" strings, add argument `raw=True`.

.. doctest::

  >>> default = block.get_mmcif_category('_software')
  >>> raw = block.get_mmcif_category('_software', raw=True)
  >>> for name in ['name', 'version', 'citation_id']:
  ...     print(default[name][0], raw[name][0])
  ...
  EPMR EPMR
  False .
  None ?

.. _set_mmcif_category:

`set_mmcif_category()` takes a category name and a dictionary (that
maps tags to lists) and creates or replaces the corresponding category:

.. doctest::

  >>> ocean_data = {'id': [1, 2, 3], 'name': ['Atlantic', 'Pacific', 'Indian']}
  >>> block.set_mmcif_category('_ocean', ocean_data)

It also takes the optional `raw` parameter (default: `False`).
Specify `raw=True` if the values are already quoted, otherwise
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


Finally, we have a variant of `init_loop` for working with mmCIF categories:

.. doctest::

  >>> loop = block.init_mmcif_loop('_ocean.', ['id', 'name'])
  >>> loop.add_row(['1', 'Atlantic'])

The subtle difference from generic `init_loop` is that if the block
has other name-value pairs in the same category (say, `_ocean.depth 8.5`)
`init_loop` leaves them untouched, but `init_mmcif_loop` removes them.
Additionally, like in other `_mmcif_` functions, the trailing dot
in the category name may be omitted (but the leading underscore is required).

Validation
==========

A CIF document can conform to a dictionary (ontology, think DTD for XML
or JSON Schema for JSON). A dictionary, written in one of the versions
of DDL (Dictionary Definition Language), is itself a CIF document.
There are three versions of DDL:

* DDL1 is the simplest. It is used, for instance, for small molecule CIFs.
* DDL2 is used for PDBx/mmCIF, with activity in this area
  centered around the PDB.
* DDLm is a newer version (from around 2011) from the IUCr's COMCIFS
  (Committee for the Maintenance of the CIF Standard). It's not widely
  used yet and, like CIF2, is not supported by Gemmi.

Gemmi is primarily used in structural biology and is mostly exercised
with mmCIF and DDL2. DDL1 is supported to a limited extent (which could be
expanded if there was a good use case).

.. note::

   In most cases, it's simpler to use the command-line program
   :ref:`gemmi validate <gemmi-validate>` instead of the functions
   described below. If you use mmCIF-like files, make sure you read
   :ref:`notes about DDL2 <DDL2>`.

The validation capabilities are implemented in class `cif::Ddl`.
Let's start with a simple example, a pet weighting experiment:

.. doctest::

  >>> pet_example = '''\
  ...   data_pets
  ...   loop_
  ...    _pet_id
  ...    _pet_species
  ...    _pet_weight
  ...    1 parrot 2
  ...    2 dog    15
  ... '''

Now let's create a contrived DDL1 dictionary for it:

.. doctest::

  >>> pet_ddl = cif.read_string('''\
  ...  data_pet_index
  ...    _name               '_pet_id'
  ...    _category           pet
  ...    _type               numb
  ...
  ...  data_pet_species
  ...    _name               '_pet_species'
  ...    _category           pet
  ...    _type               char
  ...    loop_ _enumeration  parrot cat dog
  ...
  ...  data_pet_weight
  ...    _name               '_pet_weight'
  ...    _category           pet
  ...    _type               numb
  ...    _enumeration_range  0.0:100.0
  ...    _units              kg
  ... ''')

The `Ddl` class must be first supplied with a dictionary and can then validate
CIF files.

.. doctest::

  >>> validator = cif.Ddl(logger=sys.stdout)
  >>> validator.read_ddl(pet_ddl)
  >>> validator.validate_cif(cif.read_string(pet_example))
  True

Now let's append a line that will trigger errors:

.. doctest::

  >>> pet_example += '''
  ...    3 hippo  3000
  ... '''
  >>> validator.validate_cif(cif.read_string(pet_example))
  string:2 [pets] _pet_species: hippo is not one of the allowed values:
    parrot
    cat
    dog
  string:2 [pets] _pet_weight: value out of expected range: 3000
  False

Errors are sent to a logger as described in a :ref:`separate section <logger>`.
The logger is set in the constructor and can be changed at any point:

.. doctest::

  >>> validator.set_logger((None, 0))
  >>> validator.validate_cif(cif.read_string(pet_example))
  False

Calling `read_ddl()` moves the content of a `Document` to the `Ddl` class,
leaving the original object empty (it's slightly faster this way).
`read_ddl()` can be called multiple times to use multiple dictionaries
(or extensions) simultaneously.

`Ddl` has a few flags to enable or disable certain types of checks.
These correspond to the optional checks listed in the documentation
of the :ref:`gemmi validate <gemmi-validate>` subcommand.
In C++, these are member variables that can be set directly.
In Python, they are set through keyword arguments in the constructor
(except for `use_deposition_checks`, which is set directly).

The minimal example above used a contrived dictionary. Normally, you will
use a dictionary downloaded from the IUCr, wwPDB or another source --
perhaps with your own extensions. So you'll use `cif.read()` instead of
`cif.read_string()`.

.. _ddl2:

Notes on DDL2
-------------

The commonly used DDL2-based dictionaries are
`available from wwPDB <https://mmcif.wwpdb.org/dictionaries/downloads.html>`_.
To validate mmCIF files, use the current version of the PDBx/mmCIF
dictionary (`mmcif_pdbx_v50.dic` as of 2025).
The original IUCr mmCIF dictionary (`cif_mm.dic`) is now only of historical
interest. It was actively developed in the 1990s, but in the 2000s development
was taken over by the PDB under the PDBx/mmCIF name. Formally, PDBx/mmCIF
is an extension of the IUCr mmCIF, but for all practical purposes it can be
thought of as the current mmCIF dictionary.
When we talk about mmCIF files, it's shorthand for PDBx/mmCIF
or mmCIF-like files. No software targets the original mmCIF specification.

The mmCIF dictionary itself is massive---over 5MB of text---so it, too, can
use some validation. That's what the `mmcif_ddl.dic` (DDL2) dictionary is for.
This dictionary can also validate itself, closing the loop:

.. code-block:: console

  $ gemmi validate -d mmcif_ddl.dic mmcif_ddl.dic

The PDBx/mmCIF dictionary can be used to validate coordinate,
structure factor and chemical component files.
From a validation perspective, they are all the same thing.
That's why all categories in mmCIF files are, according to the dictionary,
optional (`_category.mandatory_code no`) -- we can't tell what the file
must contain. However, many items are marked as mandatory within
categories (`_item.mandatory_code`). If you've ever wondered
about the difference between null values `?` and `.` in mmCIF files:
the PDB's software writes `?` and `.` for optional and mandatory items,
respectively (an implementation detail that deviates from the CIF 1.1 spec).
If an item is *mandatory*, it only means that if its category is present,
the tag must also be present, but its value can be unknown or n/a.

DDL2 is missing a comprehensive specification. What is not covered in
International Tables for Crystallography, vol G (2006), has to be
inferred from studying dictionaries and asking around.
Parent-child relationships are particularly challenging.
Tags may have associated parent tags (e.g. `_entity.id` is the parent
of `_entity_poly.entity_id`), and groups of tags may have associated
parent groups (defined in the `pdbx_item_linked_group` category).
But it's unclear if every parent must exist.
The PDB's own validation software (CifCheck from
`mmcif-dict-suite <https://sw-tools.rcsb.org/apps/MMCIF-DICT-SUITE/>`_)
checks for the presence of parent tags but has a long list of arbitrary
exceptions hardcoded into the program, otherwise most of the files from
the PDB wouldn't validate. In gemmi, the relationships are not checked
by default, but there is an option for it (`-p` in `gemmi validate`).

In some cases, broken relationships are fixable. In others, there is a
fundamental mismatch between the design of the mmCIF schema and the capabilities
of DDL2. For example, some aspects of polymers and non-polymers are
described in different categories, but a residue can't be conditionally
linked to one or the other. So, it's linked only to polymeric categories,
leaving the schema partially incorrect.

In addition to `_item_type` (data type), the PDBx/mmCIF spec also features
`_pdbx_item_type` (described as an "alternate data type").
Similarly, `_item_range` (permissible range) has its twin
`_pdbx_item_range` (alternate permissible range), and
`_item_enumeration` has the alternate `_pdbx_item_enumeration`.
While the purpose of this dualism is not documented, the alternates are
interpreted as deposition constraints (typically, stricter checks).
Gemmi provides the `use_deposition_checks` option (`--depo` in the command-line
program) to apply the alternate criteria instead of the standard ones.

If you run gemmi validation in verbose mode, you might see warnings about
incorrect regular expressions in a dictionary. In general, regexes
come in various flavors. Over the years, some flavors have been formally
defined and standardized (POSIX BRE, ERE, RegExp in EcmaScript, etc.).
I think the regexes used in DDL2 are closest to POSIX ERE (Extended RegExp).
Gemmi has hacks for parsing all the regexes that have been in mmCIF dictionaries
for a long time, but sometimes new ones are added that are inconsistent with
the older ones, so full support can't be guaranteed.

Dictionaries allow only relatively simple checks.
When you deposit files to the PDB, they are primarily validated by other means.
Coordinate files are processed by a program called MAXIT,
and structure factor files -- by SF-CONVERT. These are part of a C++ codebase
that has been developed at RCSB since the late 1990s.
The files you deposit don't need to strictly conform to the dictionary,
but they must be able to pass through the processing programs, which rewrite
them anyway. Attempts to make an mmCIF file more conformant with the spec
sometimes backfire, choking the OneDep pipeline.

Validation helps spot certain types of mistakes
but shouldn't be overemphasized.
When generating an mmCIF file, the goal is to ensure that it can be
read and correctly interpreted by the software it will be used with next.
Validation against a dictionary is a guideline, not the goal.


JSON
====

`cif::Document` can be stored in a JSON format, and it can be read from
JSON file. This is in general true about CIF files - their content can
be converted to JSON and back. In gemmi we have a number of options to
customize the translation. In particular, both mmCIF and CIF-JSON flavours
are supported.
More details about the flavours are given in the description of
:ref:`gemmi cif2json <cif2json>`.

C++
---

Header `gemmi/to_json.hpp` provides code for serializing
`cif::Document` as JSON.

Such JSON files can be read back into the `cif::Document` structure
using function from `gemmi/json.hpp`.

Python
------

`Document.as_json()` returns the document serialized in JSON string.
To output mmJSON add argument `mmjson=True`.

`cif.read_mmjson_string()` reads mmJSON from a string.

mmJSON file (possibly gzipped) can be read using function `cif.read_mmjson()`.
In addition, the function `cif.read()` will also read mmJSON format
if the file name ends with `.json` or `.json.gz`.


Design choices
==============

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
but I decided to use PEGTL for this task (and I'm very happy with this choice).
I was convinced by the
`TAOC++ JSON <https://github.com/taocpp/json>`_
parser that is based on PEGTL and has a good balance of simplicity
and performance.

`PEGTL <https://github.com/taocpp/PEGTL/>`_ is a C++ library (not a generator)
for creating PEG parsers.
PEG stands for Parsing Expression Grammar -- a simpler approach than
the traditional Context Free Grammar.

As a result, our parser depends on a third-party (header-only) library,
but the parser itself is pretty simple.

Data structures
---------------

The next thing is how we store the data read from file.
We decided to rely on the C++ standard library where we can.

Generally, storage of such data involves (in C++ terms) some containers
and a union/variant type for storing values of mixed types.

We use primarily `std::vector` as a container.

A custom structure `Item` with (unrestricted) union is used as
a variant-like class.

Strings are stored in `std::string` and it is fast enough.
Mainstream C++ standard libraries have short string optimization (SSO)
for up to 15 or 22 characters, which covers most of the values in mmCIF files.

Performance
===========

Gemmi has `the fastest <https://github.com/project-gemmi/mmcif-benchmark>`_
open-source CIF parser (at least in the hands of the author).
While further improvement would be possible (some JSON parsers are
`much faster <https://github.com/project-gemmi/benchmarking-json>`_
and parsing CIF and JSON is not that different),
it is not a priority.

.. _cif_examples:

Examples
========

The examples here use C++ or Python.
Full working code can be found in the examples__ directory.

The examples below can be run on one or more PDBx/mmCIF files.
The ones that perform PDB-wide analysis are meant to be run on a
`local copy <https://www.wwpdb.org/download/downloads>`_ of the mmCIF
archive (30GB+ gzipped, don't uncompress!).

__ https://github.com/project-gemmi/gemmi/tree/master/examples

mmCIF to XYZ
------------

To start with something simple, let us convert mmCIF
to the `XYZ format <https://en.wikipedia.org/wiki/XYZ_file_format>`_:

.. literalinclude:: code/cif_cc.cpp
   :language: cpp
   :lines: 20-32


.. _auth_label_example:

auth vs label
-------------

When you look at the list of atoms (`_atom_site.*`) in mmCIF files
some columns seem to be completely redundant. Are they?

.. raw:: html

    <table style="font-size:11px; border-spacing:5px 0; margin-top:60px; width:100%; border-collapse:separate; overflow-wrap: normal;">
    <colgroup>
    <col span="3">
    <col class="orange-bg">
    <col>
    <col class="yellow-bg">
    <col class="cyan-bg">
    <col>
    <col class="greenish-bg">
    <col span="12">
    <col class="greenish-bg">
    <col class="yellow-bg">
    <col class="cyan-bg">
    <col class="orange-bg">
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
where `_atom_site.auth_atom_id` differs from `_atom_site.label_atom_id`, or
where `_atom_site.auth_comp_id` differs from `_atom_site.label_comp_id`.
So here is a small C++ program:

.. literalinclude:: ../examples/auth_label.cpp

We compile it, run it, and come back after an hour:

.. code-block:: none

    $ g++ -O2 -Iinclude examples/auth_label.cpp src/gz.cpp -lz
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

Update: the result above is from 2017.
In the meantime, the PDB removed the differences.
Tags auth_atom_id and auth_comp_id are now completely redundant
(they always have been, except for mistakes).

Amino acid frequency
--------------------

.. highlight:: python

After seeing a `paper <https://doi.org/10.1093/nar/gkw978>`_
in which amino acid frequency was averaged over 5000+ proteomes
(Ala 8.76%, Cys 1.38%, Asp 5.49%, etc)
we may want to compare it with the frequency in the PDB database.
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
we can recalculate `_chem_comp.formula_weight` from `_chem_comp.formula`.
The full code is in :file:`examples/weights.py`.
It includes a function that converts `'C2 H6 O2'` to `{C:2, H:6, O:2}`.
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

Now let us try to re-calculate `_entity.formula_weight` from the chem_comp
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
But to show something else, here we will use `gemmi grep`, a little
utility that is documented in a :ref:`separate section <grep>`.

First we try how to extract interesting data from a single entry:

.. code-block:: console

    $ gemmi grep --delimiter=' ' _struct_conn.conn_type_id \
    > -a _struct_conn.ptnr1_label_comp_id -a _struct_conn.ptnr1_label_atom_id \
    > -a _struct_conn.ptnr2_label_comp_id -a _struct_conn.ptnr2_label_atom_id \
    > 5CBL
    5CBL disulf CYS SG BME S2
    5CBL covale ILE CD1 BME C2

Then we pipe the output through Unix shell utilities.
`| grep disulf` limits the output to the disulfide bonds.
``| awk '{ print $3, $4 "\n" $5, $6 }'`` changes each line into two;
the first output line above becomes:

.. code-block:: none

    CYS SG
    BME S2

Then we run it on the whole PDB archive, sort, count and print the
statistics. The complete command is:

.. code-block:: console

    $ gemmi grep --delimiter=' ' _struct_conn.conn_type_id \
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

And what other bond types are annotated in `_struct_conn`?

.. code-block:: console

    $ gemmi grep -O -b _struct_conn.conn_type_id /hdd/mmCIF/ | sort | uniq -c
     501851 covale
     184231 disulf
    4035115 hydrog
    1394732 metalc


mmJSON-like data
----------------

Gemmi has a built-in support for mmJSON and comes with
converters :ref:`cif2json <cif2json>` and :ref:`json2cif <json2cif>`,
but just as an exercise let us convert mmJSON to mmCIF in Python:

.. literalinclude:: ../examples/from_json.py
   :lines: 5-

.. _ccd_example:

Chemical Component Dictionary
-----------------------------

For something a bit different, let us look at the data from
the :file:`components.cif` from `CCD <https://www.wwpdb.org/data/ccd>`_.
This file describes all the monomers (residues,
ligands, solvent molecules) from the PDB entries.

As an exercise, let us check heavy atoms in selenomethionine:

.. doctest::

    >>> from gemmi import cif
    >>> cif.read('components.cif')  # doctest: +ELLIPSIS
    <gemmi.cif.Document with ... blocks (000, 001, 002...)>
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

Or which one has the largest number of heavy atoms:

.. doctest::

    >>> el_tag = '_chem_comp_atom.type_symbol'
    >>> max(ccd, key=lambda b: sum(el != 'H' for el in b.find_values(el_tag)))
    <gemmi.cif.Block X12>
    >>> _.find_value('_chem_comp.formula')
    '"C96 H153 N31 O25"'

The :file:`components.cif` file is big, so we may want to split it into
multiple files. As an example, here we create a new document with only
one block, and we write it to a file:

.. doctest::

    >>> block = ccd['X12']
    >>> d = cif.Document()
    >>> d.add_copied_block(block)
    <gemmi.cif.Block X12>
    >>> d.write_file('X12.cif')

This example can be expanded to write multiple selected blocks
to a file, as in `examples/sub_ccd.py`.
To write a single block, we could simply call:

.. doctest::

    >>> ccd['X12'].write_file('X12.cif')

In the next example we delete things we do not need.
Let us write only components on letter A to a new file.
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
