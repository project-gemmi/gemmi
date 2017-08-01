#!/usr/bin/env python
import sys
from gemmi import cif

greeted = set()
for path in sys.argv[1:]:
    try:
        doc = cif.read_file(path)  # copy all the data from mmCIF file
        block = doc.sole_block()  # mmCIF has exactly one block
        for s in block.find_loop("_atom_site.type_symbol"):
            if s not in greeted:
                print("Hello " + s)
                greeted.add(s)
    except Exception as e:
        print("Oops. %s" % e)
        sys.exit(1)
