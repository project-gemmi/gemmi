#!/usr/bin/env python
# Check column presence and order in the _atom_site category.
# In mmCIF v5 esd _atom_site.*esd columns were removed.

from __future__ import print_function
import sys
from gemmi import cif, CifWalk

#ESD = 'Cartn_x_esd Cartn_y_esd Cartn_z_esd occupancy_esd B_iso_or_equiv_esd '
ESD = ''
USUAL_ORDER = ('group_PDB id type_symbol label_atom_id label_alt_id '
               'label_comp_id label_asym_id label_entity_id label_seq_id '
               'pdbx_PDB_ins_code Cartn_x Cartn_y Cartn_z occupancy '
               'B_iso_or_equiv ' + ESD + 'pdbx_formal_charge '
               'auth_seq_id auth_comp_id auth_asym_id auth_atom_id '
               'pdbx_PDB_model_num')
counts = {}
for arg in sys.argv[1:]:
    for path in CifWalk(arg):
        block = cif.read(path).sole_block()
        loop_tags = block.find_loop("_atom_site.id").get_loop().tags
        assert all(t.startswith("_atom_site.") for t in loop_tags)
        tags = ' '.join(t[11:] for t in loop_tags)
        if tags != USUAL_ORDER:
            print(tags)
            print(USUAL_ORDER)
            print(block.name, tags)
            counts[tags] = counts.get(tags, 0) + 1

for key, value in counts.items():
    print(value, key)

# Results: in v4 a few EM structures (5A9Z 5AA0 5FKI 4UDF)
# had different order, with ATOM/HETATM in the middle.
