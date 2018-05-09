#!/usr/bin/env python3

# Read CCD and fill-in one-letter-codes in gemmi/resinfo.hpp
# Usage: ./tools/resinfo.py > resinfo.hpp-new

import re
from sys import stderr
from gemmi import cif

ccd = cif.read('components.cif')

STANDARD = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLN', 'GLU',
            'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
            'PRO', 'SER', 'THR', 'TRP', 'TYR', 'UNK', 'VAL',
            'SEC', 'PYL',
            'A', 'C', 'G', 'I', 'U',
            'DA', 'DC', 'DG', 'DI', 'DT', 'DU']

pattern = re.compile("""case ID\("([A-Z0-9]+)"\): return .*'(.)',( +\d+) """)
for line in open('include/gemmi/resinfo.hpp'):
    m = pattern.search(line)
    if m:
        name = m.group(1)
        block = ccd[name]
        olc = block.find_value('_chem_comp.one_letter_code')
        if olc != '?':
            code = olc[0]
            if name not in STANDARD:
                code = code.lower()
            pos = m.start(2)
            line = line[:pos] + code + line[pos+1:]
        if len(olc) > 1:
            print('One-letter-codes for %s: %s' % (name, olc), file=stderr)
        parent = block.find_value('_chem_comp.mon_nstd_parent_comp_id')
        if parent and parent != '?' and ',' not in parent:
            parent_olc = ccd[parent].find_value('_chem_comp.one_letter_code')
            if parent_olc != olc:
                print(name, olc, '<>', parent, parent_olc, file=stderr)
        old_nh = int(m.group(3))
        formula = cif.as_string(block.find_value('_chem_comp.formula'))
        obj = re.search(r'\bH(\d+)\b', formula)
        if obj:
            nh = int(obj.group(1))
            line = line[:m.end(2)] + ("', %3d" % nh) + line[m.end(3):]
    print(line, end='')
