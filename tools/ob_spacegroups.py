#!/usr/bin/env python
# Check the space-groups.txt from Open Babel
# Usage: ob_spacegroups.py /path/to/space-groups.txt

from __future__ import print_function
import sys
import gemmi

space_groups_txt = sys.argv[1]
chunks = open(space_groups_txt).read().split('\n\n')
def parse_chunk(lines):
    return {'number': int(lines[0]),
            'hall': lines[1].strip(),
            'xhm': lines[2].strip(),
            'symops': [gemmi.Op(line) for line in lines[3:]]}
data = [parse_chunk(c.splitlines()) for c in chunks if c]
print(len(data), 'items')
for d in data:
    ops = gemmi.symops_from_hall(d['hall'])
    assert len(ops.sym_ops) * len(ops.cen_ops) == len(d['symops'])
    given = set(d['symops'])
    generated = set(ops)
    if given != generated:
        print(d['hall'])
        print('common:',  '  '.join(x.triplet() for x in given & generated))
        print('given:    ', '  '.join(x.triplet() for x in given - generated))
        print('generated:', '  '.join(x.triplet() for x in generated - given))
for d in data:
    names = d['xhm'].split(',')
    if all(gemmi.find_spacegroup_by_name(name) is None for name in names):
        print('not in gemmi:', names)
