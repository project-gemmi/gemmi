#!/usr/bin/env python
# This example shows how to put pythonic data structure into an mmCIF file.
# Two paths are required as arguments: input (mmJSON file) and output.
from __future__ import print_function

import json
import sys
from gemmi import cif

file_in, file_out = sys.argv[1:]

with open(file_in) as f:
    json_data = json.load(f)
assert len(json_data) == 1  # data_1ABC
(block_name, block_data), = json_data.items()
assert block_name.startswith('data_')
# Now block_data is a dictionary that maps category names to dictionaries
# that in turn map column names to lists with values.
doc = cif.Document()
block = doc.add_new_block(block_name[5:])
for cat, data in block_data.items():
    block.set_mmcif_category('_'+cat, data)
doc.write_file(file_out)
