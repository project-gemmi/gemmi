#!/usr/bin/env python3
# Show differences between two Refmac intermediate files.
# Usage: ./rstdiff.py first.rst second.rst

import sys
from gemmi import cif
from collections import namedtuple

COLUMNS = ['record', 'number', 'label', 'period', 'atom_id_1', 'atom_id_2',
           'atom_id_3', 'atom_id_4', 'value', 'dev', 'val_obs']

Restraint = namedtuple('Restraint', COLUMNS)

def read_rst(path):
    result = []
    for row in cif.read(path)['restraints'].find('_restr.', COLUMNS):
        if row[0] in ['MONO', 'LINK']:
            result.append((row[0].lower(), row[2].lower(), []))
        else:
            data = Restraint(*[x.lower() if x != '.' else None for x in row])
            result[-1][2].append(data)
    return result

def main():
    file1, file2 = sys.argv[1:]
    r1 = read_rst(file1)
    r2 = read_rst(file2)
    mono_count_1 = sum(t[0] == 'mono' for t in r1)
    mono_count_2 = sum(t[0] == 'mono' for t in r2)
    if mono_count_1 == mono_count_2:
        print('MONO count same:', mono_count_1)
    else:
        print('MONO count differs:', mono_count_1, mono_count_2)
    link_count_1 = sum(t[0] == 'link' for t in r1)
    link_count_2 = sum(t[0] == 'link' for t in r2)
    if link_count_1 == link_count_2:
        print('LINK count same:', link_count_1)
    else:
        print('LINK count differs:', link_count_1, link_count_2)
    # compare up to the first major difference
    def nums_differ(a, b, eps=0.002):
        return a != b and abs(float(a) - float(a)) > eps
    for n, (a, b) in enumerate(zip(r1, r2)):
        if a[0] != b[0] or a[1] != b[1]:
            print('item', n, 'differs:', a[:2], 'vs', b[:2])
            if a[0] != b[0] or a[1] != '.':
                print('stop')
                break
        if len(a[2]) != len(b[2]):
            print('item %d: %d vs %d restr.' % (n, len(a[2]), len(b[2])))
        else:
            for m, (restr1, restr2) in enumerate(zip(a[2], b[2])):
                if restr1[:8] != restr2[:8] or \
                   any(nums_differ(restr1[i], restr2[i]) for i in (8, 9, 10)):
                    print('Different restraint %d:%d:\n%s\nvs\n%s\n' %
                          (n, m, restr1, restr2))

main()
