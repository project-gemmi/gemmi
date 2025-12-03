Miscellaneous utils
###################

FASTA and PIR reader
====================

Gemmi provides a function to parse two sequence file formats, FASTA and PIR.
The function takes a string containing the file's content as an argument:

.. doctest::

  >>> with open('P0C805.fasta') as f:
  ...     fasta_str = f.read()
  >>> gemmi.read_pir_or_fasta(fasta_str)  #doctest: +ELLIPSIS
  [<gemmi.FastaSeq object at 0x...>]

The string must start with a header line that begins with `>`.
In the case of the PIR format, which starts with `>P1;` (or F1, DL, DC, RL, RC,
or XX instead of P1), the next line is also part of the header.
The sequence file may contain multiple sequences, each preceded by a header.
Whitespace in a sequence is ignored, except for blank lines,
which are only allowed between sequences.
A sequence can contain letters, dashes, and residue names in parentheses.
The latter is an extension inspired by the format used in mmCIF files,
in which non-standard residues are given in parentheses, e.g., `MA(MSE)GVN`.
The sequence may end with `*`.

`FastaSeq` objects, returned from `read_pir_or_fasta()`,
contain only two strings:

.. doctest::

  >>> (fasta_seq,) = _
  >>> fasta_seq.header
  'sp|P0C805|PSMA3_STAA8 Phenol-soluble modulin alpha 3 peptide OS=Staphylococcus aureus (strain NCTC 8325 / PS 47) OX=93061 GN=psmA3 PE=1 SV=1'
  >>> fasta_seq.seq
  'MEFVAKLFKFFKDLLGKFLGNN'

.. _logger:

Logger
======

Logger is a tiny helper class for passing messages from a gemmi function
to the calling function.
The messages being passed are usually info or warnings that a command-line
program would print to stdout or stderr.

The Logger has two member variables:

.. literalinclude:: ../include/gemmi/logger.hpp
   :language: cpp
   :start-at: ///
   :end-at: int threshold

and a few member functions for sending messages.

When a function takes a Logger argument, we can pass:

.. tab:: C++

 * `{&Logger::to_stderr}` to redirect messages to stderr
   (to_stderr() calls fprintf),
 * `{&Logger::to_stdout}` to redirect messages to stdout,
 * `{&Logger::to_stdout, 3}` to print only warnings (threshold=3),
 * `{nullptr, 0}` to disable all messages,
 * `{}` to throw exception on error and ignore other messages
   (the default, see Quirk above),
 * `{[](const std::string& s) { do_anything(s);}}` to do anything else.

.. tab:: Python

 * `sys.stderr` or `sys.stdout` or any other stream (an object with `write`
   and `flush` methods), to redirect messages to that stream,
 * `(sys.stdout, 3)` to print only warnings (threshold=3),
 * `(None, 0)` to disable all messages,
 * `None` to raise exception on error and ignore other messages
   (the default, see Quirk above),
 * a function that takes a message string as its only argument
   (e.g. `lambda s: print(s.upper())`).


.. _pdb_dir:

Copy of the PDB archive
=======================

Some examples in this documentation work with a local copy
of the Protein Data Bank archive. This section describes
our setup and functions for working with such a setup.

As in BioJava, we assume that the `$PDB_DIR` environment variable
points to a directory containing `structures/divided/mmCIF` -- the same
arrangement as in the
`PDB Archive <https://www.wwpdb.org/ftp/pdb-ftp-sites>`_.

.. code-block:: console

    $ cd $PDB_DIR
    $ du -sh structures/*/*  # as of Jun 2017
    34G    structures/divided/mmCIF
    25G    structures/divided/pdb
    101G   structures/divided/structure_factors
    2.6G   structures/obsolete/mmCIF

The PDB recommends using rsync for bulk file downloads.
We can keep the local copy up-to-date using, for instance, such a script:

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

Gemmi provides a helper function `expand_if_pdb_code()` for using the local
archive copy. This function is used by many subcommands of the gemmi
:ref:`program <program>`, so they can take a PDB code instead of a path.

When called with a PDB code (case insensitive), `expand_if_pdb_code()`
expands it to the path of the corresponding:

* coordinate mmCIF file (if the second argument is absent or 'M'),
* PDB file ('P'),
* or structure factor mmCIF file ('S').

.. doctest::

  >>> os.environ['PDB_DIR'] = '/copy'
  >>> gemmi.expand_if_pdb_code('1ABC', 'P') # PDB file
  '/copy/structures/divided/pdb/ab/pdb1abc.ent.gz'
  >>> gemmi.expand_if_pdb_code('1abc', 'M') # mmCIF file
  '/copy/structures/divided/mmCIF/ab/1abc.cif.gz'
  >>> gemmi.expand_if_pdb_code('1abc', 'S') # SF-mmCIF file
  '/copy/structures/divided/structure_factors/ab/r1abcsf.ent.gz'

If the first argument is not a PDB code,
the function returns it unchanged:

.. doctest::

  >>> arg = 'file.cif'
  >>> gemmi.is_pdb_code(arg)
  False
  >>> gemmi.expand_if_pdb_code(arg, 'M')
  'file.cif'


Directory walking
=================

.. note::

  This functionality was developed primarily for C++ before `std::filesystem`
  and is kept for backward compatibility.
  In new code consider using `std::filesystem::recursive_directory_iterator`
  in C++ and `os.walk` in Python.

Many of the utilities and examples developed for this project
work with archives of files, such as CIF files from wwPDB or COD.
To make it easier to iterate over all files of a particular type
in a directory tree, we provide classes that end with `Walk`,
such as `CifWalk`.

.. tab:: C++

 .. code-block:: cpp

  #include <gemmi/dirwalk.hpp>

  // ...
  // throws std::system_error if directory doesn't exist
  for (const std::string& cif_file : gemmi::CifWalk(top_dir)) {
    auto doc = gemmi::read_cif_gz(cif_file);
    // ...
  }

.. tab:: Python

 .. doctest::

  >>> import gemmi
  >>> # throws FileNotFoundError if directory doesn't exist
  >>> for path in gemmi.CifWalk('../tests/'):
  ...     doc = gemmi.cif.read(path)
  ...     # ...
  ...

Here are all the classes in this category.

* `DirWalk` (C++ only) -- returns all files. It's actually a template
  and the items below are instantiations of this template class,
  for instance: `using CifWalk = DirWalk<true, impl::IsCifFile>`.
  But it's simpler to just call them classes.
* `CifWalk` -- returns `*.cif` and `r*.ent` files (the latter matches
  SF-mmCIF files from the PDB archive). Here, and in all the classes below
  (except for `GlobWalk`), `.gz` at the end is ignored, and comparison
  is case insensitive. For example, `foo.CiF.gZ` counts as `*.cif`.
* `MmCifWalk` (C++ only) -- matches `*.cif` or `*.mmcif`.
* `PdbWalk` (C++ only) -- matches `*.pdb` or `*.ent`, but not `r*.ent`.
* `CoorFileWalk` -- combines the two above, except for `*-sf.mmcif`.
* `GlobWalk` (C++ only) -- can use any glob pattern to match file names.
  For example: `GlobWalk(top_level_path, "*.json")`.

If the user doesn't have permission to read one of the traversed directories,
the function will raise an error (`std::runtime_error` / `RuntimeError`).

All these directory walking functions are powered by the
`tinydir <https://github.com/cxong/tinydir>`_ library
(a single-header library copied into `include/gemmi/third_party`).

Input abstractions (decompressing on the fly)
=============================================

Gemmi can decompress `.gz` files on the fly.

File reading is decoupled from content parsing.
We use input classes that read content from a plain file,
a gzipped file (uncompressing it on the fly) or from a memory buffer.

This is done by using polymorphic input classes with virtual methods
for reading bytes or words or lines. Gemmi supports both text-based file formats
(PDB, CIF, XDS ACII, mmJSON) and binary ones (MRC/CCP4, MTZ).
As a side note: initially, for the first few years of development,
gemmi was a header-only library. As the library grew, this became
less and less practical, and gemmi started a gradual transition to a conventional
split into .hpp headers and .cpp implementation files.
The low-level way of reading a file, for example, a CIF file, is:

 .. code-block:: cpp

  cif::read(gemmi::BasicInput{mmcif_file});

and a gzipped file can be read with:

 .. code-block:: cpp

  cif::read(gemmi::MaybeGzipped{mmcif_gzipped_file});

Similarly with other formats, for example:

 .. code-block:: cpp

  xds.read_input(gemmi::MaybeGzipped(input_path));

BasicInput and MaybeGzipped are little helpers used to provide access to the actual
input stream class: one of FileStream, MemoryStream, GzStream, all of which are
derived from AnyStream (see `gemmi/input.hpp`).

The above approach used to be recommended and documented.
But after switching to a compiled library, this recommendation changed.
Gemmi now has more user-friendly wrapper functions, such as:

 .. code-block:: cpp

  read_structure_gz(mmcif_gzipped_file);

which are compiled as part of the library. Use these wrappers instead
of explicit MaybeGzipped (which still works)
to make your code easier to read, understand and faster to compile.
That's what will be documented through the rest of this manual,
but you may still find `MaybeGzipped` in older code that uses Gemmi,
so I wanted to mention it here.

Notes
-----

When the library was header-only, we avoided the overhead of virtual
functions by using static (compile-time) polymorphism in the input (FileStream,
MemoryStream, GzStream) classes.
At that time, linking with the Zlib library was necessary only for programs using
MaybeGzipped.
Now, the input classes (\*Stream) use virtual functions;
the overhead of the virtual calls turned out to be negligible.
As a side-effect, libgemmi (and  therefore any program using it)
must be linked with zlib or zlib-ng.

Compilation is much faster, though.
With static polymorphism, the parsing code had to be compiled in every
compilation unit that used `cif::read`, and functions with different
parameters (`BasicInput` vs `MaybeGzipped`) resulted in separate template
instantiations, which optimized performance at the cost of binary size bloat.

At first, we thought we might need more types of inputs (other compression
formats) later on, and exposing details such as `gemmi::MaybeGzipped` could help
the user grasp how to implement other inputs (without changing gemmi).
This concern turned out to be unimportant (no one used another compression
format during 8 years of Gemmi development), so using wrapper functions
is preferable, as they are easier to read and understand.

If you compile gemmi yourself, note that linking it with zlib-ng instead of zlib
(as done in official binaries on PyPI) makes reading gzipped files noticeably faster.
