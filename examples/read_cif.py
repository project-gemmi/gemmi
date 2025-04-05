#!/usr/bin/env python
import sys
from gemmi import cif

for path in sys.argv[1:]:
    try:
        doc = cif.read_file(path)  # copy all the data from mmCIF file
        for block in doc:  # iterate over blocks
            for cc in block.find("_chem_comp.", ["id", "formula_weight"]):
                print(cc[0], "weights", cc[1])
    except Exception as e:
        print("Oops. %s" % e)
