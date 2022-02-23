#!/usr/bin/env python
# Reads a mask, reads a coordinate file, and checks how the mask compares
# with a bulk solvent mask for the coordinate file.

# To check if we produce identical maps as cctbx, do:
#   cctbx.python -m mmtbx.command_line.mask file.pdb
#   maskcheck.py mask.ccp4 file.pdb

import os
import sys
import numpy
import gemmi

def maskcheck(mask_path, coor_path, output_diff_map=None):
    # read mask
    mask = gemmi.read_ccp4_mask(mask_path, setup=True)
    grid = mask.grid

    # read coordinates
    st = gemmi.read_structure(coor_path)

    # check if unit cell and symmetry are the same
    if not grid.unit_cell.approx(st.cell, epsilon=1e-4):
        print('Cell parameters differ. PDB:', st.cell.parameters)
        print('                       Mask:', grid.unit_cell.parameters)
    sg = st.find_spacegroup()
    if grid.spacegroup != sg:
        print('Space groups differ:', grid.spacegroup.xhm(), 'vs', sg.xhm())

    # prepare a mask with the same size
    if 'REFMAC' in os.environ:
        print('Generating REFMAC-compatible mask ...')
        masker = gemmi.SolventMasker(gemmi.AtomicRadiiSet.Refmac)
        masker.island_min_volume = 8 * 2.8**3
    else:
        print('Generating CCTBX-compatible mask ...')
        masker = gemmi.SolventMasker(gemmi.AtomicRadiiSet.Cctbx)
    grid2 = gemmi.Int8Grid()
    grid2.copy_metadata_from(grid)
    masker.put_mask_on_int8_grid(grid2, st[0])

    compare_mask_arrays(grid, grid2)
    #print_nearby_atoms(st, grid, grid2)
    if output_diff_map:
        write_diff_map(grid, grid2, output_diff_map)


def compare_mask_arrays(grid1, grid2):
    arr1 = numpy.array(grid1, copy=False)
    arr2 = numpy.array(grid2, copy=False)
    if arr1.shape != arr2.shape:
        sys.exit('Different grid sizes %s and %s. Exiting.' %
                 (arr1.shape, arr2.shape))
    print('Size: %d x %d x %d,' % arr1.shape, 'total', arr1.size, 'points')
    t = 2 * (arr1 != 0) + (arr2 != 0)
    print('File-Gemmi Count Fraction')
    for (a, b) in [(0, 0), (1, 1), (0, 1), (1, 0)]:
        n = numpy.count_nonzero(t == 2*a+b)
        print('%d-%d %12d %6.2f%%' % (a, b, n, 100.*n/arr1.size))

# Print nearby atom for each differing point
def print_nearby_atoms(st, grid1, grid2):
    ns = gemmi.NeighborSearch(st[0], st.cell, 4).populate()
    for (p1, p2) in zip(grid1, grid2):
        if p1.value != p2.value:
            pos = grid2.point_to_position(p2)
            mark = ns.find_nearest_atom(pos)
            cra = mark.to_cra(st[0])
            print('%d-%d near %s' % (p1.value, p2.value, cra))

def write_diff_map(grid1, grid2, output_diff_map):
    arr_out = numpy.array(grid1, copy=False) - numpy.array(grid2, copy=False)
    map_out = gemmi.Ccp4Map()
    map_out.grid = gemmi.FloatGrid(arr_out.astype(dtype=numpy.float32))
    map_out.grid.copy_metadata_from(grid1)
    map_out.update_ccp4_header()
    map_out.write_ccp4_map(output_diff_map)


if __name__ == '__main__':
    if len(sys.argv) not in (3, 4):
        sys.exit("Usage: maskcheck.py MAP_FILE COORDINATE_FILE [OUTPUT_MAP]")
    maskcheck(*sys.argv[1:])
