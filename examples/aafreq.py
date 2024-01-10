#!/usr/bin/env python
# Check amino-acid frequency in the PDB database (or it's subset)
# by reading meta-data from mmCIF files.

import sys
from collections import Counter
from gemmi import cif, CifWalk

totals = Counter()
for arg in sys.argv[1:]:
    for path in CifWalk(arg, try_pdbid='M'):
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
        print(block.name, *(f'{m}:{c}' for (m, c) in aa_counter.most_common()))
# finally, print the total counts as percentages
f = 100.0 / sum(totals.values())
print('TOTAL', *(f'{m}:{c*f:.2f}%' for (m, c) in totals.most_common(20)))
