#!/usr/bin/env python

import sys
import numpy
import gemmi

def maskdiff(path1, path2):
    mask1 = gemmi.read_ccp4_mask(path1)
    mask1.setup()
    arr1 = numpy.array(mask1.grid, copy=False)
    mask2 = gemmi.read_ccp4_mask(path2)
    mask2.setup()
    arr2 = numpy.array(mask2.grid, copy=False)
    print("Size: %d x %d x %d  and  %d x %d x %d" % (arr1.shape + arr2.shape))
    if arr1.shape != arr2.shape:
        sys.exit("Different sizes. Exiting.")
    t = 2 * (arr1 != 0) + (arr2 != 0)
    for (a, b) in [(0, 0), (1, 1), (0, 1), (1, 0)]:
        n = numpy.count_nonzero(t == 2*a+b)
        print('%d-%d %12d %6.2f%%' % (a, b, n, 100.*n/arr1.size))
    print('total %10d' % arr1.size)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit("Usage: maskdiff.py map1 map2")
    maskdiff(sys.argv[1], sys.argv[2])
