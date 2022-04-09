#!/usr/bin/env python3
# Show differences between two Refmac intermediate files.
# The intermediate files should have content of both crd and rst files.
# Usage: ./crd-diff.py first.crd second.crd

import argparse
from collections import namedtuple
from itertools import zip_longest
from gemmi import cif

COLUMNS = ['record', 'number', 'label', 'period',
           'atom_id_1', 'atom_id_2', 'atom_id_3', 'atom_id_4',
           'value', 'dev', 'val_obs']

Restraint = namedtuple('Restraint', COLUMNS)
Crd = namedtuple('Crd', ['atoms', 'real_serial'])
ignore_tors = False

def read_rst(doc):
    result = []
    for row in doc['restraints'].find('_restr.', COLUMNS):
        if row[0] in ['MONO', 'LINK']:
            result.append((row[0].lower(), row[2].lower(), []))
        else:
            data = Restraint(*[x.lower() if x != '.' else None for x in row])
            if ignore_tors and data.record == 'tors':
                continue
            result[-1][2].append(data)
    return [r for r in result if r[2]]

def read_crd(doc):
    block = doc[0]
    sites = block.find('_atom_site.', ['id', 'label_atom_id', 'label_alt_id',
                                       'label_comp_id', 'occupancy',
                                       'calc_flag'])
    atoms = [a for a in sites if a[-1] != 'M']
    len_diff = len(sites) - len(atoms)
    if len_diff != 0:
        print('Removed %d atoms with calc_flag=M' % len_diff)
    real_serial = {None: None, '.': '.'}
    for a in atoms:
        real_serial[a[0]] = len(real_serial) - 1
    return Crd(atoms, real_serial)

def get_crd_atom_row(crd, atom_id):
    idx = crd.real_serial[atom_id] - 1
    return crd.atoms[idx]

# hydrogen distances in Refmac may not be ideal
def can_have_wrong_val_obs(crd, rst1):
    def is_hydrogen(atom_id):
        atom = get_crd_atom_row(crd, atom_id)
        return cif.as_string(atom[1]).startswith('H')
    if rst1.record == 'bond':
        return is_hydrogen(rst1.atom_id_2)
    elif rst1.record == 'angl':
        return is_hydrogen(rst1.atom_id_3)
    elif rst1.record == 'tors':
        return is_hydrogen(rst1.atom_id_4)
    elif rst1.record == 'plan':
        # we can't easily check here if there are H's in the plane
        return True
        #return is_hydrogen(rst1.atom_id_1)
    return False

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', action='store_true', help='verbose output')
    parser.add_argument('--no-tors', action='store_true',
                        help='ignore torsion angles')
    parser.add_argument('file1_crd', metavar='file1.crd')
    parser.add_argument('file2_crd', metavar='file2.crd')
    args = parser.parse_args()
    if args.no_tors:
        global ignore_tors
        ignore_tors = True

    doc1 = cif.read(args.file1_crd)
    doc2 = cif.read(args.file2_crd)

    crd1 = read_crd(doc1)
    crd2 = read_crd(doc2)
    if len(crd1.atoms) != len(crd2.atoms):
        print('_atom_site count differs: %d vs %d' %
              (len(crd1.atoms), len(crd2.atoms)))
    for a1, a2 in zip(crd1.atoms, crd2.atoms):
        if any(a1.str(i) != a2.str(i) for i in [1, 2, 3, 5]):
            print('First difference:')
            print('ATOM %s %s %s %s' % (a1[0], a1.str(1), a1[2], a1[3]))
            print('ATOM %s %s %s %s' % (a2[0], a2.str(1), a2[2], a2[3]))
            return
        if float(a1[4]) != float(a2[4]):
            print('ATOM %s/%s %s %s %s occupancy %s vs %s' % (
                  a1[0], a2[0], a1.str(1), a1[2], a1[3], a1[4], a2[4]))

    r1 = read_rst(doc1)
    r2 = read_rst(doc2)
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

    def fmt_atom_id(atom_id, crd):
        if atom_id is None:
            return '.'
        atom = get_crd_atom_row(crd, atom_id)
        alt = ''
        if atom[2] != '.':
            alt = '.' + atom[2]
        return '%s(%s %s%s)' % (atom_id, atom[3], atom[1], alt)

    def fmt(r, crd):
        return '%s %s %s %s  %s : %s : %s : %s   %s ± %s  (%s)' % (
               r.record, r.number, r.label or '.', r.period or '.',
               fmt_atom_id(r.atom_id_1, crd), fmt_atom_id(r.atom_id_2, crd),
               fmt_atom_id(r.atom_id_3, crd), fmt_atom_id(r.atom_id_4, crd),
               r.value, r.dev, r.val_obs)

    def has_same_ids(rst1, rst2):
        return (rst1.record == rst2.record
                and all(crd1.real_serial[u] == crd2.real_serial[v]
                        for u, v in zip(rst1[4:8], rst2[4:8])))

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
            print('Item %d (%s) has different restr count: %d vs %d' %
                  (n, a[0], len(a[2]), len(b[2])))
        n_rst1 = 0
        n_rst2 = 0
        while n_rst1 < len(a[2]) or n_rst2 < len(b[2]):
            r_str = '%s/%s restraint %d:%d/%d' % (b[0], b[1], n, n_rst1, n_rst2)
            if n_rst1 >= len(a[2]):
                rst2 = b[2][n_rst2]
                print('ADDED record %s:\n%s\n' % (r_str, fmt(rst2, crd2)))
                n_rst2 += 1
                continue
            if n_rst2 >= len(b[2]):
                rst1 = a[2][n_rst1]
                print('REMOVED record %s:\n%s\n' % (r_str, fmt(rst1, crd1)))
                n_rst1 += 1
                continue
            rst1 = a[2][n_rst1]
            rst2 = b[2][n_rst2]
            same_ids = has_same_ids(rst1, rst2)
            if not same_ids:
                a2 = a[2]
                b2 = b[2]
                if any(n_rst2+k < len(b2) and has_same_ids(rst1, b2[n_rst2+k])
                       for k in range(1,10)):
                    print('ADDED record %s:\n%s\n' % (r_str, fmt(rst2, crd2)))
                    n_rst2 += 1
                    continue
                if any(n_rst1+k < len(a2) and has_same_ids(a2[n_rst1+k], rst2)
                       for k in range(1,10)):
                    print('REMOVED record %s:\n%s\n' % (r_str, fmt(rst1, crd1)))
                    n_rst1 += 1
                    continue
            if rst1.record != rst2.record:
                print('Different records %s:\n%s\nvs\n%s\n' %
                      (r_str, fmt(rst1, crd1), fmt(rst2, crd2)))
            elif (rst1.label != rst2.label if rst1.record != 'chir'
                  else rst1.label[0] != rst2.label[0]):
                print('Different labels in %s:\n%s\nvs\n%s\n' %
                      (r_str, fmt(rst1, crd1), fmt(rst2, crd2)))
            elif rst1.period != rst2.period:
                print('Different period in %s:\n%s\nvs\n%s\n' %
                      (r_str, fmt(rst1, crd1), fmt(rst2, crd2)))
            elif not same_ids:
                print('Different atom id in %s:\n%s\nvs\n%s\n' %
                      (r_str, fmt(rst1, crd1), fmt(rst2, crd2)))
            elif not same_nums(rst1.value, rst2.value):
                print('Different value for %s (%s vs %s) in:\n%s\n' %
                      (r_str, rst1.value, rst2.value, fmt(rst1, crd1)))
            elif not same_nums(rst1.dev, rst2.dev):
                print('Different dev for %s (%s vs %s) in:\n%s\n' %
                      (r_str, rst1.dev, rst2.dev, fmt(rst1, crd1)))
            elif not (same_nums(rst1.val_obs, rst2.val_obs,
                                eps=val_obs_eps(rst1.record),
                                mod360=(rst1.record == 'tors'))
                      or can_have_wrong_val_obs(crd1, rst1)):
                print('val_obs differs for %s (%s vs %s) in:\n%s\n' %
                      (r_str, rst1.val_obs, rst2.val_obs, fmt(rst1, crd1)))
            #elif rst1.number != rst2.number:
            #    print('Serial number differs in %s: %s -> %s' %
            #          (r_str, rst1[1], rst2[1]))
            n_rst1 += 1
            n_rst2 += 1

def val_obs_eps(record):
    if record == 'tors':
        return 0.15
    elif record == 'plan':
        return 0.015
    else:
        return 0.004

if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        pass
