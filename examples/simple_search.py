#!/usr/bin/env python
"Find PDB entries with more than 50,000 anisotropic B-factors."

from __future__ import print_function
import sys
from gemmi import cif, CifWalk

for path in CifWalk(sys.argv[1]):
    block = cif.read(path).sole_block()
    anis = block.find_values("_atom_site_anisotrop.id")
    if len(anis) > 50000:
        print(block.name, len(anis))
