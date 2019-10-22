#!/usr/bin/env python3

import sys
import gemmi

CCD_PATH = 'components.cif.gz'
OUTPUT_PATH = 'sub_ccd.cif'
COORDINATE_FILES = sys.argv[1:]

# Make a list of residue names that we need.
mon_names = set()
for coordinate_file in COORDINATE_FILES:
    st = gemmi.read_structure(coordinate_file)
    mon_names.update(st[0].get_all_residue_names())

# Copy blocks corresponding to these residues to a new file.
sub = gemmi.cif.Document()
for block in gemmi.cif.read(CCD_PATH):
    if block.name in mon_names:
        sub.add_copied_block(block)
sub.write_file(OUTPUT_PATH)
