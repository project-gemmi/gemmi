"utilities that simplify examples"

import os
import sys
import argparse


def get_file_paths_from_args():
    """\
    Process arguments as filenames or directories with .cif(.gz) files,
    and yield the file paths.
    Normally we first test our scripts on a few files:
      ./myscript 1mru.cif another.cif
    and then do pdb-wide analysis:
      ./myscript $PDB_DIR/structures/divided/mmCIF
    If we may have a list of PDB codes to analyze (one code per line, the code
    must be the first word, but may be followed by others), we do:
      ./myscript -f my-list.txt $PDB_DIR/structures/divided/mmCIF
    Sometimes we want to check a specific PDB codes. If $PDB_DIR is set
    it is enough to do:
      ./myscript 1ABC 2def
    """
    parser = argparse.ArgumentParser(usage='%(prog)s [options] path [...]')
    parser.add_argument('path', nargs='+', help=argparse.SUPPRESS)
    parser.add_argument('--only', metavar='LIST',
                        help='Use only files that match names in this file')
    args = parser.parse_args()
    only = None
    if args.only:
        with open(args.only) as list_file:
            only = set(line.split()[0].lower() for line in list_file
                       if line.strip())
    for arg in args.path:
        if os.path.isdir(arg):
            for root, name, extlen in sorted_cif_search(arg):
                if not only or name[:-extlen].lower() in only:
                    yield os.path.join(root, name)
        elif len(arg) == 4 and arg.isalnum():
            pdb_dir = os.getenv('PDB_DIR')
            if not pdb_dir:
                sys.exit('Error: $PDB_DIR not set, where to look for ' + arg)
            yield os.path.join(pdb_dir, 'structures', 'divided', 'mmCIF',
                               arg[1:3].lower(), arg.lower() + '.cif.gz')
        else:
            yield arg


def sorted_cif_search(top_dir):
    for root, dirs, files in os.walk(top_dir):
        dirs.sort()
        for name in sorted(files):
            if name.endswith('.cif'):
                yield root, name, 4
            elif name.endswith('.cif.gz'):
                yield root, name, 7


def formula_to_dict(formula):
    '"O4 P -3" -> {O:4, P:1}'
    fdict = {}
    for elnum in formula.split():
        na = sum(e.isalpha() for e in elnum)
        if na == len(elnum):
            fdict[elnum] = 1
        elif na != 0:
            fdict[elnum[:na]] = int(elnum[na:])
    return fdict
