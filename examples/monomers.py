#!/usr/bin/env python
from __future__ import print_function
import os
import argparse
from collections import Counter
from gemmi import cif
import util


def to_formula(cnt):
    '{C:12, O:8, N:1} -> "C12 N O8"'
    return ' '.join(k + (str(v) if v != 1 else '')
                    for k, v in sorted(cnt.items()))


def get_monomer_cifs(mon_path):
    'yield all files mon_path/?/*.cif in alphabetic order'
    for d in sorted(s for s in os.listdir(mon_path) if len(s) == 1):
        subdir_path = os.path.join(mon_path, d)
        if os.path.isdir(subdir_path):
            for name in sorted(os.listdir(subdir_path)):
                if name.endswith('.cif'):
                    yield os.path.join(subdir_path, name)


def check_formulas(ccd):
    '''\
    Print the chemical components in CCD (components.cif) that have
    formula not consistent with the list of atoms.
    '''
    for b in ccd:
        atoms = Counter(a.as_str(0).upper() for a in
                        b.find('_chem_comp_atom.type_symbol'))
        formula = cif.as_string(b.find_value('_chem_comp.formula')).upper()
        fdict = util.formula_to_dict(formula)
        if fdict != atoms:
            print('[%s]' % b.name, formula, '<>', to_formula(atoms))


def compare_monlib_with_ccd(mon_path, ccd):
    'compare monomers from monomer library and CCD that have the same names'
    PRINT_MISSING_ENTRIES = False
    for path in get_monomer_cifs(mon_path):
        mon = cif.read_any(path)
        for mb in mon:
            if mb.name in ('', 'comp_list'):
                continue
            assert mb.name.startswith('comp_')
            name = mb.name[5:]
            cb = ccd.find_block(name)
            if cb:
                compare_chem_comp(mb, cb, name)
            elif PRINT_MISSING_ENTRIES:
                print('Not in CCD:', name)


def compare_chem_comp(mb, cb, name):
    mon_atoms = Counter(a.as_str(0).upper() for a in
                        mb.find('_chem_comp_atom.type_symbol'))
    ccd_atoms = Counter(a.as_str(0).upper() for a in
                        cb.find('_chem_comp_atom.type_symbol'))
    del mon_atoms['H']
    del ccd_atoms['H']
    if mon_atoms != ccd_atoms and mon_atoms + Counter('O') != ccd_atoms:
        print(name, ':', to_formula(mon_atoms), '<>', to_formula(ccd_atoms))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('ccd_path', metavar='/path/to/components.cif[.gz]')
    args = parser.parse_args()
    ccd = cif.read_any(args.ccd_path)
    check_formulas(ccd)
    mon_path = os.getenv('CLIBD_MON')
    if mon_path:
        compare_monlib_with_ccd(mon_path, ccd)
    else:
        print('$CLIBD_MON not set!')


main()
