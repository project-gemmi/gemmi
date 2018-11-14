#!/usr/bin/env python3
# Show differences between two Refmac intermediate files.
# Usage: ./rstdiff.py first.rst second.rst

import argparse
from collections import namedtuple
from itertools import zip_longest
from gemmi import cif

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
    return [r for r in result if r[2]]

def read_crd(path):
    block = cif.read(path).sole_block()
    sites = block.find('_atom_site.', ['id', 'label_atom_id', 'label_alt_id',
                                       'label_comp_id', 'occupancy',
                                       'calc_flag'])
    atoms = [a for a in sites if a[-1] != 'M']
    real_serial = {None: None, '.': '.'}
    for a in atoms:
        real_serial[a[0]] = len(real_serial) - 1
    return Crd(atoms, real_serial)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', action='store_true', help='verbose output')
    parser.add_argument('file1_rst', metavar='file1.rst')
    parser.add_argument('file2_rst', metavar='file2.rst')
    args = parser.parse_args()

    crd1 = read_crd(args.file1_rst[:-4] + '.crd')
    crd2 = read_crd(args.file2_rst[:-4] + '.crd')
    if len(crd1.atoms) != len(crd2.atoms):
        print('_atom_site count differs: %d vs %d' %
              (len(crd1.atoms), len(crd2.atoms)))
    for a1, a2 in zip(crd1.atoms, crd2.atoms):
        if any(a1.str(i) != a2.str(i) for i in [1, 2, 3, 5]):
            print('First difference:')
            print('ATOM %s %s %s %s' % (a1[0], a1.str(1), a1[2], a1[3]))
            print('ATOM %s %s %s %s' % (a2[0], a2.str(1), a2[2], a2[3]))
            break
        if float(a1[4]) != float(a2[4]):
            print('ATOM %s %s %s %s occupancy %s vs %s' % (
                  a1[0], a1.str(1), a1[2], a1[3], a1[4], a2[4]))

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

    def same_nums(a, b, eps=0.003, mod360=False):
        if a == b:
            return True
        abs_diff = abs(float(a) - float(b))
        return abs_diff < eps or (mod360 and abs(abs_diff - 360) < eps)

    def fmt_atom_id(atom_id, n):
        if atom_id is None:
            return '.'
        crd = (crd1 if n == 1 else crd2)
        idx = crd.real_serial[atom_id] - 1
        atom = crd.atoms[idx]
        alt = ''
        if atom[2] != '.':
            alt = '.' + atom[2]
        return '%s(%s %s%s)' % (atom_id, atom[3], atom[1], alt)

    def fmt(r, n=1):
        return '%s %s %s %s  %s : %s : %s : %s   %s ± %s  (%s)' % (
               r.record, r.number, r.label or '.', r.period or '.',
               fmt_atom_id(r.atom_id_1, n), fmt_atom_id(r.atom_id_2, n),
               fmt_atom_id(r.atom_id_3, n), fmt_atom_id(r.atom_id_4, n),
               r.value, r.dev, r.val_obs)

    for n, (a, b) in enumerate(zip_longest(r1, r2, fillvalue=[None, None, []])):
        if a[:2] != b[:2]:
            print('item', n, 'differs:', a[:2], 'vs', b[:2])
            if a[0] != b[0]:
                print('stop')
                break
        else:
            if args.v:
                print('---', n, a[0], a[1])
        if len(a[2]) != len(b[2]):
            print('Item %d has different restr count: %d vs %d' %
                  (n, len(a[2]), len(b[2])))
        fillvalue = Restraint(*[None]*11)
        for m, (rst1, rst2) in enumerate(zip_longest(a[2], b[2],
                                                     fillvalue=fillvalue)):
            r_str = '%s/%s restraint %d:%d' % (b[0], b[1], n, m)
            if rst1.record != rst2.record:
                print('Different records %s:\n%s\nvs\n%s\n' %
                      (r_str, fmt(rst1), fmt(rst2, 2)))
            elif (rst1.label != rst2.label if rst1.record != 'chir'
                  else rst1.label[0] != rst2.label[0]):
                print('Different labels in %s:\n%s\nvs\n%s\n' %
                      (r_str, fmt(rst1), fmt(rst2, 2)))
            elif rst1.period != rst2.period:
                print('Different period in %s:\n%s\nvs\n%s\n' %
                      (r_str, fmt(rst1), fmt(rst2, 2)))
            elif not all(crd1.real_serial[u] == crd2.real_serial[v]
                         for u, v in zip(rst1[4:8],rst2[4:8])):
                print('Different atom id in %s:\n%s\nvs\n%s\n' %
                      (r_str, fmt(rst1), fmt(rst2, 2)))
            elif not same_nums(rst1.value, rst2.value):
                print('Different value for %s (%s vs %s) in:\n%s\n' %
                      (r_str, rst1.value, rst2.value, fmt(rst1)))
            elif not same_nums(rst1.dev, rst2.dev):
                print('Different dev for %s (%s vs %s) in:\n%s\n' %
                      (r_str, rst1.dev, rst2.dev, fmt(rst1)))
            elif not same_nums(rst1.val_obs, rst2.val_obs,
                               eps=val_obs_eps(rst1.record),
                               mod360=(rst1.record == 'tors')):
                print('Different val_obs for %s (%s vs %s) in:\n%s\n' %
                      (r_str, rst1.val_obs, rst2.val_obs, fmt(rst1)))
            elif rst1.number != rst2.number:
                print('Serial number differs in %s: %s -> %s' %
                      (r_str, rst1[1], rst2[1]))

def val_obs_eps(record):
    if record == 'tors':
        return 0.15
    elif record == 'plan':
        return 0.015
    else:
        return 0.004

main()
