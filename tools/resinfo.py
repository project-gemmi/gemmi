#!/usr/bin/env python3

# Read CCD and fill-in one-letter-codes in gemmi/resinfo.hpp
# Usage: ./tools/resinfo.py > resinfo.hpp-new

import re
from sys import stderr
from gemmi import cif

ccd = cif.read('components.cif')

pattern = re.compile("""case ID\("([A-Z0-9]+)"\): return .*'(.)',""")
for line in open('include/gemmi/resinfo.hpp'):
    m = pattern.search(line)
    if m:
        pos = m.start(2)
        name = m.group(1)
        block = ccd[name]
        olc = block.find_value('_chem_comp.one_letter_code')
        if olc != '?':
            line = line[:pos] + olc[0] + line[pos+1:]
        if len(olc) > 1:
            print('One-letter-codes for %s: %s' % (name, olc), file=stderr)
        parent = block.find_value('_chem_comp.mon_nstd_parent_comp_id')
        if parent and parent != '?' and ',' not in parent:
            parent_olc = ccd[parent].find_value('_chem_comp.one_letter_code')
            if parent_olc != olc:
                print(name, olc, '<>', parent, parent_olc, file=stderr)
    print(line, end='')
