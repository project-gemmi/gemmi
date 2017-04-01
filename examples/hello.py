#!/usr/bin/env python
import sys
import gemmi as cif

greeted = set()
if len(sys.argv) != 2: sys.exit(1)
try:
  doc = cif.Document(sys.argv[1]) # open mmCIF file & copy all the data from it
  block = doc.sole_block() # mmCIF has exactly one block
  for s in block.find_loop("_atom_site.type_symbol"):
    if s not in greeted:
      print("Hello " + s)
      greeted.add(s)
except Exception as e:
  print("Oops. %s" % e)
  sys.exit(1)
