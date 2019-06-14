#!/usr/bin/env python
"Find PDB entries with more than 50,000 anisotropic B-factors."

from __future__ import print_function
from gemmi import cif
from util import get_file_paths_from_args

for path in get_file_paths_from_args():
    block = cif.read(path).sole_block()
    anis = block.find_values("_atom_site_anisotrop.id")
    if len(anis) > 50000:
        print(block.name, len(anis))
