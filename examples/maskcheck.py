#!/usr/bin/env python
# Reads a mask, reads a coordinate file, and checks how the mask compares
# with a bulk solvent mask for the coordinate file.

# To check if we produce identical maps as cctbx, do:
#   cctbx.python -m mmtbx.command_line.mask file.pdb
#   maskcheck.py mask.ccp4 file.pdb

import sys
import numpy
import gemmi

def maskcheck(mask_path, coor_path):
    # read mask
    mask = gemmi.read_ccp4_mask(mask_path)
    mask.setup()
    grid = mask.grid

    # read coordinates
    st = gemmi.read_structure(coor_path)

    # check if unit cell and symmetry are the same
    if grid.unit_cell.parameters != st.cell.parameters:
        print('Cell parameters differ:',
              grid.unit_cell.parameters, 'vs', st.cell.parameters)
    sg = st.find_spacegroup()
    if grid.spacegroup != sg:
        print('Space groups differ:', grid.spacegroup.xhm(), 'vs', sg.xhm())

    # prepare a mask with the same size
    masker = gemmi.SolventMasker(gemmi.AtomicRadiiSet.Cctbx)
    grid2 = gemmi.Int8Grid()
    grid2.copy_metadata_from(grid)
    masker.put_mask_on_int8_grid(grid2, st[0])

    compare_mask_arrays(grid, grid2)


def compare_mask_arrays(grid1, grid2):
    arr1 = numpy.array(grid1, copy=False)
    arr2 = numpy.array(grid2, copy=False)
    if arr1.shape != arr2.shape:
        sys.exit('Different grid sizes %s and %s. Exiting.' %
                 (arr1.shape, arr2.shape))
    else:
        print('Size: %d x %d x %d,' % arr1.shape, 'total', arr1.size, 'points')
    t = 2 * (arr1 != 0) + (arr2 != 0)
    for (a, b) in [(0, 0), (1, 1), (0, 1), (1, 0)]:
        n = numpy.count_nonzero(t == 2*a+b)
        print('%d-%d %12d %6.2f%%' % (a, b, n, 100.*n/arr1.size))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit("Usage: maskcheck.py COORDINATE_FILE MAP_FILE")
    maskcheck(sys.argv[1], sys.argv[2])
