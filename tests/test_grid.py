#!/usr/bin/env python

import math
import os
import unittest
import gemmi
try:
    import numpy
except ImportError:
    numpy = None

class TestFloatGrid(unittest.TestCase):
    def test_reading(self):
        path = os.path.join(os.path.dirname(__file__), '5i55_tiny.ccp4')
        m = gemmi.read_ccp4_map(path)
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


# In 5a11 applying NCS causes atom clashing
FRAGMENT_5A11 = """\
CRYST1   48.367   89.613   83.842  90.00 101.08  90.00 P 1 21 1      4          
MTRIX1   1 -0.999980  0.006670  0.000050       13.74419                         
MTRIX2   1  0.006260  0.941760 -0.336220       35.56956                         
MTRIX3   1 -0.002290 -0.336210 -0.941780      205.34927                         
ATOM    256  SG  CYS A  37      -1.002 -31.125  88.394  1.00 14.38           S  
ATOM   2969  SG  CYS B  37      14.582 -23.455 132.554  1.00 18.14           S  
"""  # noqa: W291 - trailing whitespace

# In 1gtv two different chains (with partial occupancy) are exactly in the
# same place
FRAGMENT_1GTV = """\
CRYST1   76.353   76.353  134.815  90.00  90.00 120.00 P 65 2 2     24
ATOM    635  SG  CYS A  85      42.948   6.483  17.913  0.48 23.86           S
ATOM   2293  SG  CYS B  85      42.948   6.483  17.913  0.52 23.86           S
"""

class TestNeighborSearch(unittest.TestCase):
    def test_5a11(self, use_populate=True):
        st = gemmi.read_pdb_string(FRAGMENT_5A11)
        a1 = st[0].sole_residue('A', gemmi.SeqId(37, ' '))[0]
        ns = gemmi.NeighborSearch(st[0], st.cell, 5)
        if use_populate:
            ns.populate()
        else:
            for n_ch, chain in enumerate(st[0]):
                for n_res, res in enumerate(chain):
                    for n_atom, atom in enumerate(res):
                        ns.add_atom(atom, n_ch, n_res, n_atom)
        marks = ns.find_atoms(a1.pos, a1.altloc, 3)
        m1, m2 = sorted(marks, key=lambda m: ns.dist(a1.pos, m.pos()))
        self.assertAlmostEqual(ns.dist(a1.pos, m1.pos()), 0, delta=5e-6)
        self.assertAlmostEqual(ns.dist(a1.pos, m2.pos()), 0.13, delta=5e-3)
        cra2 = m2.to_cra(st[0])
        self.assertEqual(cra2.chain.name, 'B')
        self.assertEqual(str(cra2.residue.seqid), '37')
        self.assertEqual(cra2.atom.name, 'SG')
        marks2 = ns.find_neighbors(a1, 0.1, 3)
        self.assertEqual(len(marks2), 1)
        self.assertEqual(marks2[0], m2)

    def test_5a11_using_add_atom(self):
        self.test_5a11(use_populate=False)

    def test_1gtv(self):
        st = gemmi.read_pdb_string(FRAGMENT_1GTV)
        a1 = st[0].sole_residue('A', gemmi.SeqId(85, ' '))[0]
        ns = gemmi.NeighborSearch(st[0], st.cell, 5)
        ns.populate()
        marks = ns.find_atoms(a1.pos, a1.altloc, 3)
        self.assertEqual(len(marks), 2)
        for mark in marks:
            d = ns.dist(a1.pos, mark.pos())
            self.assertAlmostEqual(d, 0, delta=5e-6)
        marks2 = ns.find_neighbors(a1, 0.1, 3)
        self.assertEqual(len(marks2), 0)

class TestContactSearch(unittest.TestCase):
    def test_contact_search(self):
        cs = gemmi.ContactSearch(4.0)
        hg = gemmi.Element('Hg')
        self.assertEqual(cs.get_radius(hg), 0)
        cs.setup_atomic_radii(1, 0)
        cs.set_radius(hg, 1.5)
        self.assertEqual(cs.get_radius(hg), 1.5)


if __name__ == '__main__':
    unittest.main()
