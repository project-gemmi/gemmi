#!/usr/bin/env python

import math
import unittest
import gemmi
from common import full_path
try:
    import numpy
except ImportError:
    numpy = None

class TestFloatGrid(unittest.TestCase):
    def test_reading(self):
        m = gemmi.read_ccp4_map(full_path('5i55_tiny.ccp4'))
        self.assertEqual(m.grid.nu, 8)
        self.assertEqual(m.grid.nv, 6)
        self.assertEqual(m.grid.nw, 10)
        self.assertEqual(m.header_i32(28), 0)
        m.set_header_i32(28, 20140)  # set NVERSION
        self.assertEqual(m.header_i32(28), 20140)
        dmax = m.header_float(21)
        self.assertEqual(dmax, max(p.value for p in m.grid))
        self.assertNotEqual(m.grid.axis_order, gemmi.AxisOrder.XYZ)
        m.setup()
        self.assertEqual(m.grid.axis_order, gemmi.AxisOrder.XYZ)
        self.assertEqual(m.grid.nu, 60)
        self.assertEqual(m.grid.nv, 24)
        self.assertEqual(m.grid.nw, 60)
        self.assertEqual(m.grid.point_count, 60 * 24 * 60)
        self.assertEqual(m.header_float(14), 90.0)  # 14 - alpha angle
        self.assertEqual(m.grid.unit_cell.alpha, 90.0)
        self.assertEqual(m.grid.spacegroup.ccp4, 4)  # P21

        pos = gemmi.Position(19.4, 3., 21.)
        frac = m.grid.unit_cell.fractionalize(pos)
        pos_value = 2.1543798446655273
        self.assertAlmostEqual(m.grid.interpolate_value(pos), pos_value)
        self.assertAlmostEqual(m.grid.interpolate_value(frac), pos_value)

        # this spacegroup has symop -x, y+1/2, -z
        m.grid.set_value(60-3, 24//2+4, 60-5, 100)  # image of (3, 4, 5)
        self.assertEqual(m.grid.get_value(60-3, 24//2+4, 60-5), 100)
        self.assertTrue(math.isnan(m.grid.get_value(3, 4, 5)))
        m.grid.symmetrize_max()
        self.assertEqual(m.grid.get_value(3, 4, 5), 100)
        m.grid.set_value(3, 4, 5, float('nan'))
        self.assertTrue(math.isnan(m.grid.get_value(3, 4, 5)))
        m.grid.symmetrize_min()
        self.assertEqual(m.grid.get_value(3, 4, 5), 100)
        m.grid.set_value(60-3, 24//2+4, 60-5, float('nan'))
        m.grid.symmetrize_max()
        self.assertEqual(m.grid.get_value(60-3, 24//2+4, 60-5), 100)
        if numpy:
            arr = numpy.array(m.grid, copy=False)
            self.assertEqual(arr.shape, (60, 24, 60))
            self.assertEqual(arr[3][4][5], 100)
            grid2 = gemmi.FloatGrid(arr)
            self.assertTrue(numpy.allclose(m.grid, grid2, atol=0.0, rtol=0,
                                           equal_nan=True))

    def test_new(self):
        N = 24
        m = gemmi.FloatGrid(N, N, N)
        self.assertEqual(m.nu, N)
        self.assertEqual(m.nv, N)
        self.assertEqual(m.nw, N)
        m.set_value(1,2,3, 1.0)
        self.assertEqual(m.sum(), 1.0)
        m.spacegroup = gemmi.find_spacegroup_by_name('C2')
        self.assertEqual(m.spacegroup.number, 5)
        m.symmetrize_max()
        self.assertEqual(m.sum(), 4.0)
        m.get_point(0, 0, 0).value += 1
        self.assertEqual(m.sum(), 5.0)
        m.fill(2.0)
        m.spacegroup = gemmi.find_spacegroup_by_name('P 62 2 2')
        self.assertEqual(len(m.spacegroup.operations()), 12)
        m.set_value(1, 2, 3, 0.0)
        m.symmetrize_min()
        self.assertEqual(m.sum(), 2 * N * N * N - 2 * 12)


if __name__ == '__main__':
    unittest.main()
