#!/usr/bin/env python
from __future__ import print_function
from collections import Counter
from gemmi import cif

# To keep this example small we moved handling of command-line args to util.py.
# When a directory is given as an argument get_file_paths_from_args()
# yields all the cif(.gz) paths under this directory.
from util import get_file_paths_from_args

# Check amino-acid frequency in the PDB database
totals = Counter()
for path in get_file_paths_from_args():
    # read file (uncompressing on the fly) and get the only block
    block = cif.read(path).sole_block()
    # find table with the sequence
    seq = block.find('_entity_poly_seq.', ['entity_id', 'mon_id'])
    # convert table with chain types (protein/DNA/RNA) to dict
    entity_types = dict(block.find('_entity_poly.', ['entity_id', 'type']))
    # and count these monomers that correspond to a protein chain
    aa_counter = Counter(row.str(1) for row in seq
                         if 'polypeptide' in entity_types[row.str(0)])
    totals += aa_counter
    # print residue counts for each file
    print(block.name, *('%s:%d' % c for c in aa_counter.most_common()))
# finally, print the total counts as percentages
f = 100.0 / sum(totals.values())
print('TOTAL', *('%s:%.2f%%' % (m, c*f) for (m, c) in totals.most_common(20)))
