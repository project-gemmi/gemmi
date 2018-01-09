#!/usr/bin/env python

import os
import unittest
import gemmi

class TestFloatGrid(unittest.TestCase):
    def test_reading(self):
        path = os.path.join(os.path.dirname(__file__), '5i55_tiny.ccp4')
        m = gemmi.read_ccp4_map(path)
        self.assertEqual(m.nu, 8)
        self.assertEqual(m.nv, 6)
        self.assertEqual(m.nw, 10)
        self.assertEqual(m.header_i32(28), 0)
        m.set_header_i32(28, 20140)  # set NVERSION
        self.assertEqual(m.header_i32(28), 20140)
        dmax = m.header_float(21)
        self.assertEqual(dmax, max(m))
        m.setup()
        self.assertEqual(m.nu, 60)
        self.assertEqual(m.nv, 24)
        self.assertEqual(m.nw, 60)
        self.assertEqual(m.header_float(14), 90.0)  # 14 - alpha angle
        self.assertEqual(m.unit_cell.alpha, 90.0)

    def test_new(self):
        N = 24
        m = gemmi.FloatGrid(N, N, N)
        self.assertEqual(m.nu, N)
        self.assertEqual(m.nv, N)
        self.assertEqual(m.nw, N)
        m.set_value(1,2,3, 1.0)
        self.assertEqual(sum(m), 1.0)
        m.space_group = gemmi.find_spacegroup_by_name('C2')
        self.assertEqual(m.space_group.number, 5)
        m.symmetrize_max()
        self.assertEqual(sum(m), 4.0)
        m.fill(2.0)
        m.space_group = gemmi.find_spacegroup_by_name('P 62 2 2')
        self.assertEqual(len(m.space_group.operations()), 12)
        m.set_value(1, 2, 3, 0.0)
        m.symmetrize_min()
        self.assertEqual(sum(m), 2 * N * N * N - 2 * 12)


if __name__ == '__main__':
    unittest.main()
