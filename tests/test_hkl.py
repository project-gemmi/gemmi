#!/usr/bin/env python

import os
import unittest
import gemmi
from common import full_path, get_path_for_tempfile
try:
    import numpy
except ImportError:
    numpy = None

class TestMtz(unittest.TestCase):
    def test_read_write(self):
        path = full_path('5e5z.mtz')
        mtz = gemmi.read_mtz_file(path)
        self.assertEqual(mtz.spacegroup.hm, 'P 1 21 1')
        out_name = get_path_for_tempfile()
        mtz.write_to_file(out_name)
        mtz2 = gemmi.read_mtz_file(out_name)
        os.remove(out_name)
        self.assertEqual(mtz2.spacegroup.hm, 'P 1 21 1')

    def fft_test(self, mtz, size):
        if numpy is None:
            return
        grid_full = mtz.get_f_phi_on_grid('FWT', 'PHWT', size, half_l=False)
        array_full = numpy.array(grid_full, copy=False)
        map1 = gemmi.transform_f_phi_grid_to_map(grid_full)
        map2 = numpy.fft.ifftn(array_full.conj())
        map2 = numpy.real(map2) * (map2.size / grid_full.unit_cell.volume)
        self.assertTrue(numpy.allclose(map1, map2, atol=5e-7, rtol=0))
        map3 = mtz.transform_f_phi_to_map('FWT', 'PHWT', size)
        self.assertTrue(numpy.allclose(map1, map3, atol=6e-7, rtol=0))

        grid2 = gemmi.transform_map_to_f_phi(map1, half_l=False)
        self.assertTrue(numpy.allclose(grid2, array_full, atol=1e-4, rtol=0))

        grid_half = mtz.get_f_phi_on_grid('FWT', 'PHWT', size, half_l=True)
        grid3 = gemmi.transform_map_to_f_phi(map1, half_l=True)
        self.assertTrue(numpy.allclose(grid3, grid_half, atol=1e-4, rtol=0))

    def test_f_phi_grid(self):
        path = full_path('5wkd_phases.mtz.gz')
        mtz = gemmi.read_mtz_file(path)
        size = mtz.get_size_for_hkl()
        for half_l in (False, True):
            grid1 = mtz.get_f_phi_on_grid('FWT', 'PHWT', size, half_l=half_l)
            grid2 = mtz.get_f_phi_on_grid('FWT', 'PHWT', size, half_l=half_l,
                                          hkl_orient=gemmi.HklOrient.LKH)
            if numpy is None:
                continue
            array1 = numpy.array(grid1, copy=False)
            array2 = numpy.array(grid2, copy=False)
            self.assertTrue((array2 == array1.transpose(2,1,0)).all())

        self.fft_test(mtz, size)


class TestSfMmcif(unittest.TestCase):
    def test_reading(self):
        doc = gemmi.cif.read(full_path('r5wkdsf.ent'))
        rblock = gemmi.as_refln_blocks(doc)[0]
        self.assertEqual(rblock.spacegroup.hm, 'C 1 2 1')

if __name__ == '__main__':
    unittest.main()
