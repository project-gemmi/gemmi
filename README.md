# Gemmi Cif Parser

C++11 header-only library for reading CIF and STAR files (ATM superset of CIF 1.1 but subset of STAR).

## What are STAR and CIF?

(in case someone comes across this page when looking for serialization format)

STAR is a human-readable data serialization format
(think XML or YAML, but also ASN.1)
widely used in molecular-structure sciences.
It is accompanied by a schema language called DDL.

CIF (Crystallographic Information File) -- a standard file format
in crystallography -- is basically a restricted variant of STAR.
It is restricted in features (to make implementation easier),
but also imposes arbitrary limits -- for example on the line length.

Both STAR and CIF (as well as DDL) have multiple versions.
We will be more specific in the following sections.

The syntax is relatively simple, but may be confusing at first.
(Note that the initial version of STAR was published by Sydney Hall in 1991 --
before XML and much before JSON and YAML, not to mention TOML).

Key-value pairs have the form:

    _year 2017  # means year = 2017

Blocks (sections) are marked with the `data_` keyword immediately
followed by a block name:

    data_tomato   # [tomato]
    _color red    # color = "red"

Importantly, unlike other mentioned here formats, STAR is designed
to concisely represent tabular data. This example from TOML docs:

    # TOML
    points = [ { x = 1, y = 2, z = 3 },
	       { x = 7, y = 8, z = 9 },
	       { x = 2, y = 4, z = 8 } ]

could be expressed as:

    # STAR/CIF
    loop_  # <- keyword that starts so-called loop
    _points.x _points.y _points.z
    1 2 3
    7 8 9
    2 4 8

Typically, long loops make most of the file content:

    1    N N   . LEU A 11  ? 0.5281 0.5618 0.5305 -0.0327 -0.0621 0.0104
    2    C CA  . LEU A 11  ? 0.5446 0.5722 0.5396 -0.0317 -0.0632 0.0080
    # many, many lines...
    5331 S SD  . MET C 238 ? 2.2952 2.3511 2.3275 -0.0895 0.0372  -0.0230
    5332 C CE  . MET C 238 ? 1.5699 1.6247 1.6108 -0.0907 0.0388  -0.0244

#### References

International Tables for Crystallography [Vol. G
(2006)](http://it.iucr.org/Ga/contents/)
describes all of the STAR, CIF 1.1, DDL1 and DDL2.
If you don't have access to it -- IUCr website has specs of
[CIF1.1](http://www.iucr.org/resources/cif/spec/version1.1)
and [DDLs](http://www.iucr.org/resources/cif/ddl).
As far as I can tell all versions of STAR spec are behind paywall
(except the expired [patent](https://patents.google.com/patent/WO1991016682A1),
but it is not useful as a spec).

Later papers describe later versions of the formats:
[STAR (2.0)](https://dx.doi.org/10.1021/ci300074v) (2012),
[DDLm](http://pubs.acs.org/doi/abs/10.1021/ci300075z) and
[dREL](http://pubs.acs.org/doi/abs/10.1021/ci300076w) (2012),
and [CIF 2.0](http://journals.iucr.org/j/issues/2016/01/00/aj5269/) (2016).


## What is parsed?

We test only:

* mmCIF files (i.e. CIF 1.1)
* DDL1 and DDL2 dictionaries (in particular `mmcif_pdbx_v40.dic`)
* and monomer library a.k.a. Refmac dictionary
* example files from COD

Later on we'll make sure that all correct CIF 1.1 files from COD
are handled properly.

The parser could be extended to handle the new features of CIF 2.0
or even the full STAR format, but we don't have a good reason to do this yet.
The same goes for DDLm/dREL.

The parser handles:

* all constructs of CIF 1.1 (including *save frames*),
* the `global_` and `stop_` keywords from STAR -- needed for Refmac
  monomer library and `mmcif_nmr-star.dic`, respectively.

Additionally:

* names and lines can have any length (like in STAR),
* unquoted strings cannot start with keywords (STAR spec is ambiguous
  about this -- see
  [StarTools doc](http://www.globalphasing.com/startools/) for details)
* obviously, there needs to be a whitespace between the last tag and the
  first value in a loop -- although it's omitted in the CIF 1.1 spec.

TODO:

* handling the EOL\;EOL long line convention

Considered:

* handling Greek letters (`\m` -> µ), accented letters (`\'o` -> ó),
  special alphabetic characters (`\%A` -> Å) and some other codes
  (`\\infty` -> ∞) as per CIF 1.1 markup conventions

CIF parsers in the small-molecules field often face incorrect syntax.
The papers about [iotbx.cif](https://doi.org/10.1107/S0021889811041161)
and [COD::CIF::Parser](http://dx.doi.org/10.1107/S1600576715022396)
enumerate 12 and 10 common error categories, respectively.
As the MX community embraced the CIF format later, the problem may be less
severe for us. So we start with relatively strict parser and will be
pragmatically relaxing it when needed.


## Technical consideration

This parser is a low-level layer used by other parts of the Gemmi projects,
such as our (to-be-done) mmCIF reader. For now there is no need to make
it accessible from other languages, so it has only a C++ API.

### Parser

Parsing formal languages is a well-researched topic in computer science.
The first versions of lex and yacc - popular tools that generate lexical
analyzers and parsers - were written in 1970's. Today many tools exist
to translate grammar rules into C/C++ code.
On the other hand many compilers and high-profile tools use hand-coded
parsers - as it is more flexible.

Looking at other STAR and CIF parsers -- some use parser generators
(COD::CIF::Parser and starlib2 use yacc/bison, iotbx.cif uses Antlr3),
others are hand-coded.

I had experience with flex/bison and Boost::Spirit, but I picked PEGTL
for this task.
(I was also considering Lemon and re2c but couldn't try everything).
PEGTL is a C++ library (not a generator) for creating PEG parsers.
PEG stands for Parsing Expression Grammar -- a simpler approach than
tradional Context Free Grammar.

As a result, our parser depends on a third-party (header-only) library,
but the parser itself is pretty simple.

### Data structures

The next thing is how we store the data read from file.
We decided rely on the C++ standard library where we can.

Generally, storage of such data involves (in C++ terms) some containers
and a union/variant type for storing values of mixed types.

As the container we use primarily `std::vector`.
In special cases we use `std::map` and `std::unordered_map`
for quicker access.
(Currently maps are used only as indexes for DDL dictionaries).

We use `std::string`, althouth it may not be optimal for performance.
`sizeof(std::string)` is typically 24 or 32 bytes, with SSO (short string
optimization) for up to at least 15 characters.
So at least the SSO covers most of the strings in the loops in mmCIF.

We do not use std::variant as it is new in C++17.
Instead, we use a custom structure with (unrestricted) union.
(But switching to a third-party C++11 implementation of a variant
could be considered.)

Many CIF readers do not store comments, as they are not part of the data.
In our case we may want to read, manipulate and save mmCIF file preserving
as much of the original file as possible. So we record the order of items,
comments and line breaks (but for now not any other whitespaces).

Such a data structure is far from optimal when reading DDL dictionaries.
On the other hand, even the biggest one -  PDBx/mmCIF (3.7MB) - gets parsed
and indexed in below 50ms. So we have no need to optimize it.

(data structure details are still evolving)


## How to use it in own program

It is a header-only library depending only on PEGTL (which is also header-only).
To use it, include the header `<gemmi/cif.hh>` in your C++ program,
make sure that Gemmi and PEGTL headers are in your project's include path,
and compile your program as C++11 or later.

The CIF parsing functionality can be used in two ways that roughly
correspond to DOM and SAX parsing:

1. Parsing file (or C++ stream or string) into a document:

        gemmi::cif::Document d;
        d.parse_file("my.cif")

   (more details later, examples will be provided with the code)

2. Go to a lower level and specify own actions corresponding to the grammar
   rules. This requires defining
   [PEGTL Actions](https://github.com/ColinH/PEGTL/blob/master/doc/Actions-and-States.md)
   for selected grammar rules from `cif.hh`.



## Utilities

This part of the Gemmi library is accompanied by two utilities:

* A validator that checks the syntax and, optionally, also ontology
  using a corresponding DDL1/DDL2 dictionary,

* a CIF-to-JSON converter -- useful for testing the parser.
  Similar converters are included in cod-tools and in JMol.
