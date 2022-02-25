#!/usr/bin/env python
# Reads a mask, reads a coordinate file, and checks how the mask compares
# with a bulk solvent mask for the coordinate file.

# To check if we produce identical maps as cctbx, do:
#   cctbx.python -m mmtbx.command_line.mask file.pdb
#   maskcheck.py mask.ccp4 file.pdb

'Usage: maskcheck.py [-v] MAP_FILE COORDINATE_FILE [OUTPUT_MAP]'

import os
import sys
import numpy
import gemmi

def maskcheck(mask_path, coor_path, output_diff_map=None, verbose=False):
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
    if verbose:
        print_nearby_atoms(st, grid, grid2)
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
    diff_grid = get_diff_grid(grid1, grid2)
    for negate in (True, False):
        blobs = gemmi.find_blobs_by_flood_fill(diff_grid, cutoff=0,
                                               min_volume=0, min_score=0,
                                               negate=negate)
        blobs.sort(key=lambda a: a.volume, reverse=True)
        print('\n%s %d blobs' % ((negate and '0-1' or '1-0'), len(blobs)))
        for blob in blobs:
            cra = ns.find_nearest_atom(blob.centroid).to_cra(st[0])
            print('    %.1f A^3 near %s' % (blob.volume, cra))

def get_diff_grid(grid1, grid2):
    arr = (grid1.array - grid2.array).astype(dtype=numpy.float32)
    return gemmi.FloatGrid(arr, grid1.unit_cell, grid1.spacegroup)

def write_diff_map(grid1, grid2, output_diff_map):
    map_out = gemmi.Ccp4Map()
    map_out.grid = get_diff_grid(grid1, grid2)
    map_out.update_ccp4_header()
    map_out.write_ccp4_map(output_diff_map)

def main():
    args = sys.argv[1:]
    verbose = ('-v' in args)
    if verbose:
        args.remove('-v')
    if len(args) not in (2, 3):
        sys.exit(__doc__)
    maskcheck(*args, verbose=verbose)

if __name__ == '__main__':
    main()
