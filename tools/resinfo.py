#!/usr/bin/env python3

# Read CCD and fill-in one-letter-codes in gemmi/resinfo.hpp
# Usage: ./tools/resinfo.py > resinfo.hpp-new

import re
from sys import stderr
from gemmi import cif

ccd = cif.read('components.cif.gz')

pattern = re.compile("""case ID\("([A-Z0-9]+)"\): return .*'(.)',""")
for line in open('include/gemmi/resinfo.hpp'):
    m = pattern.search(line)
    if m:
        pos = m.start(2)
        name = m.group(1)
        olc = ccd[name].find_value('_chem_comp.one_letter_code')
        if olc != '?':
            line = line[:pos] + olc[0] + line[pos+1:]
        if len(olc) > 1:
            print('One-letter-codes for %s: %s' % (name, olc), file=stderr)
    print(line, end='')
