#!/usr/bin/env python

# Read an archive of SF-mmCIF data and for each block print information:
# - PDB ID
# - 0-based block index
# - main category: R=_refln, D=_diffrn_refln
# - number of columns
# - number of rows
# - number of unique Miller indices (all symmetry-equivalent reflections
#   and Friedel pairs are counted as one)
# - do Miller indices repeat? N=no, F=only Friedel pairs, Y=yes

import sys
import os
import gemmi

def process(path):
    doc = gemmi.cif.read(path)
    code = os.path.basename(path)[1:5]
    for n, rblock in enumerate(gemmi.as_refln_blocks(doc)):
        loop = rblock.default_loop
        assert loop is not None, 'missing tag index_h'
        cat = loop.tags[0][1].upper()
        (data_type, unique) = gemmi.check_data_type_under_symmetry(rblock)
        if data_type == gemmi.DataType.Mean:
            r = 'N'
        elif data_type == gemmi.DataType.Anomalous:
            r = 'F'
        elif data_type == gemmi.DataType.Unmerged:
            r = 'Y'
        else:  # gemmi.DataType.Unknown - when space group is missing
            r = 'X'
        width = loop.width()
        length = loop.length()
        print(f'{code}\t{n}\t{cat}\t{width}\t{length}\t{unique}\t{r}')

def main():
    for arg in sys.argv[1:]:
        for path in gemmi.CifWalk(arg, try_pdbid='S'):
            try:
                process(path)
            except Exception as e:
                print(path, 'was not processed:', e)

main()
