Miscellaneous utils
###################

FASTA and PIR reader
--------------------

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

Gemmi Logger is a tiny helper class for passing messages from a gemmi function
to the calling function. It doesn't belong in this section, but it's
documented here because it's used in the previous subsection and I haven't found
a better spot for it.

The messages being passed are usually info or warnings that a command-line
program would print to stdout or stderr.

The Logger has two member variables:

.. literalinclude:: ../include/gemmi/logger.hpp
   :language: cpp
   :start-at: ///
   :end-at: int threshold

and a few member functions for sending messages.

When a function takes a Logger argument, we can pass:

**C++**

* `{&Logger::to_stderr}` to redirect messages to stderr
  (to_stderr() calls fprintf),
* `{&Logger::to_stdout}` to redirect messages to stdout,
* `{&Logger::to_stdout, 3}` to print only warnings (threshold=3),
* `{nullptr, 0}` to disable all messages,
* `{}` to throw errors and ignore other messages (the default, see Quirk above),
* `{[](const std::string& s) { do_anything(s);}}` to do anything else.

**Python**

* `sys.stderr` or `sys.stdout` or any other stream (an object with `write`
  and `flush` methods), to redirect messages to that stream,
* `(sys.stdout, 3)` to print only warnings (threshold=3),
* `(None, 0)` to disable all messages,
* `None` to throw errors and ignore other messages (the default, see Quirk above),
* a function that takes a message string as its only argument
  (e.g. `lambda s: print(s.upper())`).


.. _pdb_dir:

Copy of the PDB archive
=======================

Some of the examples in this documentation work with a local copy
of the Protein Data Bank archive. This subsection describes
the assumed setup and functions for working with this setup.

Like in BioJava, we assume that the `$PDB_DIR` environment variable
points to a directory that contains `structures/divided/mmCIF` -- the same
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

Gemmi has a helper function for using the local archive copy.
It takes a PDB code (case insensitive) and a symbol denoting what file
is requested: P for PDB, M for mmCIF, S for SF-mmCIF.

.. doctest::

  >>> os.environ['PDB_DIR'] = '/copy'
  >>> gemmi.expand_if_pdb_code('1ABC', 'P') # PDB file
  '/copy/structures/divided/pdb/ab/pdb1abc.ent.gz'
  >>> gemmi.expand_if_pdb_code('1abc', 'M') # mmCIF file
  '/copy/structures/divided/mmCIF/ab/1abc.cif.gz'
  >>> gemmi.expand_if_pdb_code('1abc', 'S') # SF-mmCIF file
  '/copy/structures/divided/structure_factors/ab/r1abcsf.ent.gz'

If the first argument is not in the PDB code format (4 characters for now)
the function returns the first argument.

.. doctest::

  >>> arg = 'file.cif'
  >>> gemmi.is_pdb_code(arg)
  False
  >>> gemmi.expand_if_pdb_code(arg, 'M')
  'file.cif'
