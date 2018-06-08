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
    def same_nums(a, b, eps=0.003, mod360=False):
        if a == b:
            return True
        abs_diff = abs(float(a) - float(b))
        return abs_diff < eps or (mod360 and abs(abs_diff - 360) < eps)
    for n, (a, b) in enumerate(zip(r1, r2)):
        if a[:2] != b[:2]:
            print('item', n, 'differs:', a[:2], 'vs', b[:2])
            if a[0] != b[0] or (a[0] == 'mono' and a[1] != '.'):
                print('stop')
                break
        if len(a[2]) != len(b[2]):
            print('item %d: %d vs %d restr.' % (n, len(a[2]), len(b[2])))
        else:
            for m, (rst1, rst2) in enumerate(zip(a[2], b[2])):
                is_tors = (rst1.record == 'tors')
                if rst1[:8] != rst2[:8]:
                    print('Different restraint %d:%d:\n%s\nvs\n%s\n' %
                          (n, m, rst1, rst2))
                elif not same_nums(rst1.value, rst2.value):
                    print('Different value for %d:%d (%s vs %s) in:\n%s\n' %
                          (n, m, rst1.value, rst2.value, rst1))
                elif not same_nums(rst1.dev, rst2.dev):
                    print('Different dev for %d:%d (%s vs %s) in:\n%s\n' %
                          (n, m, rst1.dev, rst2.dev, rst1))
                elif not same_nums(rst1.val_obs, rst2.val_obs,
                                   eps=(0.1 if is_tors else 0.003),
                                   mod360=is_tors):
                    print('Different val_obs for %d:%d (%s vs %s) in:\n%s\n' %
                          (n, m, rst1.val_obs, rst2.val_obs, rst1))

main()
