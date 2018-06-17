#!/usr/bin/env python3
# Show differences between two Refmac intermediate files.
# Usage: ./rstdiff.py first.rst second.rst

import argparse
import difflib
from gemmi import cif
from collections import namedtuple

COLUMNS = ['record', 'number', 'label', 'period',
           'atom_id_1', 'atom_id_2', 'atom_id_3', 'atom_id_4',
           'value', 'dev', 'val_obs']

Restraint = namedtuple('Restraint', COLUMNS)
Crd = namedtuple('Crd', ['atoms', 'real_serial'])

def read_rst(path):
    result = []
    for row in cif.read(path)['restraints'].find('_restr.', COLUMNS):
        if row[0] in ['MONO', 'LINK']:
            result.append((row[0].lower(), row[2].lower(), []))
        else:
            data = Restraint(*[x.lower() if x != '.' else None for x in row])
            result[-1][2].append(data)
    return result

def read_crd(path):
    block = cif.read(path).sole_block()
    sites = block.find('_atom_site.', ['id', 'label_atom_id', 'label_alt_id',
                                       'label_comp_id', 'calc_flag'])
    atoms = [a for a in sites if a[-1] != 'M']
    real_serial = {None: None, '.': '.'}
    for a in atoms:
        real_serial[a[0]] = len(real_serial)
    return Crd(atoms, real_serial)

def main():
    parser = argparse.ArgumentParser()
    #parser.add_argument('--crd', action='store_true',
    #                    help='check also corresponding crd files')
    parser.add_argument('file1_rst', metavar='file1.rst')
    parser.add_argument('file2_rst', metavar='file2.rst')
    args = parser.parse_args()

    crd1 = read_crd(args.file1_rst[:-4] + '.crd')
    crd2 = read_crd(args.file2_rst[:-4] + '.crd')
    if len(crd1.atoms) != len(crd2.atoms):
        print('_atom_site count differs: %d vs %d' %
              (len(crd1.atoms), len(crd2.atoms)))
    for a1, a2 in zip(crd1.atoms, crd2.atoms):
        if any(a1.str(i) != a2.str(i) for i in range(1, 5)):
            print('First difference:')
            print('ATOM %s %s %s %s' % (a1[0], a1.str(1), a1[2], a1[3]))
            print('ATOM %s %s %s %s' % (a2[0], a2.str(1), a2[2], a2[3]))
            print()
            break

    r1 = read_rst(args.file1_rst)
    r2 = read_rst(args.file2_rst)
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
            print('Item %d has different restr count: %d vs %d' %
                  (n, len(a[2]), len(b[2])))
            for line in difflib.unified_diff(
                    [' '.join(s or '.' for s in x[:10]) for x in a[2]],
                    [' '.join(s or '.' for s in x[:10]) for x in b[2]]):
                print(line)
        else:
            for m, (rst1, rst2) in enumerate(zip(a[2], b[2])):
                is_tors = (rst1.record == 'tors')
                r_str = '%s/%s restraint %d:%d' % (b[0], b[1], n, m)
                if rst1[:4] != rst2[:4]:
                    print('Different %s:\n%s\nvs\n%s\n' %
                          (r_str, rst1, rst2))
                if not all(crd1.real_serial[u] == crd2.real_serial[v]
                           for u, v in zip(rst1[4:8],rst2[4:8])):
                    print('Different atom id in %s:\n%s\nvs\n%s\n' %
                          (r_str, rst1, rst2))
                elif not same_nums(rst1.value, rst2.value):
                    print('Different value for %s (%s vs %s) in:\n%s\n' %
                          (r_str, rst1.value, rst2.value, rst1))
                elif not same_nums(rst1.dev, rst2.dev):
                    print('Different dev for %s (%s vs %s) in:\n%s\n' %
                          (r_str, rst1.dev, rst2.dev, rst1))
                elif not same_nums(rst1.val_obs, rst2.val_obs,
                                   eps=(0.15 if is_tors else 0.003),
                                   mod360=is_tors):
                    print('Different val_obs for %s (%s vs %s) in:\n%s\n' %
                          (r_str, rst1.val_obs, rst2.val_obs, rst1))

main()
