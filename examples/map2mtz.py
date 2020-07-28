#!/usr/bin/env python
# Convert CCP4 map to map coefficients in MTZ

import sys
import gemmi

RESOLUTION_LIMIT = 1.5  # set 0 for no limit

if len(sys.argv) != 3:
    sys.exit('Usage: map2mtz.py input.ccp4 output.mtz')

m = gemmi.read_ccp4_map(sys.argv[1])
m.setup()
sf = gemmi.transform_map_to_f_phi(m.grid, half_l=True)
data = sf.prepare_asu_data(dmin=RESOLUTION_LIMIT)

mtz = gemmi.Mtz(with_base=True)
mtz.spacegroup = sf.spacegroup
mtz.set_cell_for_all(sf.unit_cell)
mtz.add_dataset('unknown')
mtz.add_column('FWT', 'F')
mtz.add_column('PHWT', 'P')
mtz.set_data(data)
mtz.write_to_file(sys.argv[2])
