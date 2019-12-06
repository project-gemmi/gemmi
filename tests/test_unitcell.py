#!/usr/bin/env python

import unittest
from math import pi  # , isnan
from random import random
from gemmi import Position, Fractional, UnitCell, calculate_dihedral

class TestUnitCell(unittest.TestCase):
    def test_dummy_cell(self):
        cell = UnitCell()
        self.assertEqual([cell.a, cell.b, cell.c], [1, 1, 1])
        self.assertEqual([cell.alpha, cell.beta, cell.gamma], [90, 90, 90])
        self.assertEqual(cell.volume, 1.0)

    def test_ortho_cell(self):
        cell = UnitCell(25.14, 39.50, 45.07, 90, 90, 90)
        pos = Position(5, -6, 7)
        frac = cell.fractionalize(pos)
        self.assertAlmostEqual(frac.x, 0.198886, delta=1e-6)
        self.assertAlmostEqual(frac.y, -0.151899, delta=1e-6)
        self.assertAlmostEqual(frac.z, 0.155314, delta=1e-6)
        pos2 = cell.orthogonalize(frac)
        self.assertAlmostEqual(pos.x, pos2.x, delta=1e-12)
        self.assertAlmostEqual(pos.y, pos2.y, delta=1e-12)
        self.assertAlmostEqual(pos.z, pos2.z, delta=1e-12)
        corner = cell.orthogonalize(Fractional(1, 1, 1))
        self.assertAlmostEqual(corner.x, cell.a, delta=1e-12)
        self.assertAlmostEqual(corner.y, cell.b, delta=1e-12)
        self.assertAlmostEqual(corner.z, cell.c, delta=1e-12)

    def test_triclinic_cell(self):
        cell = UnitCell(35.996, 41.601, 45.756, 67.40, 66.90, 74.85)
        pos = Position(-15, -17, 190)
        frac = cell.fractionalize(pos)
        pos2 = cell.orthogonalize(frac)
        self.assertAlmostEqual(pos.x, pos2.x, delta=1e-12)
        self.assertAlmostEqual(pos.y, pos2.y, delta=1e-12)
        self.assertAlmostEqual(pos.z, pos2.z, delta=1e-12)
        # tested against values from uctbx:
        #  from cctbx import uctbx
        #  uc = uctbx.unit_cell((35.996, 41.601, 45.756, 67.40, 66.90, 74.85))
        #  uc.d_star_sq((-3, -2, 1))
        #  uc.d((3, 4, 5))
        self.assertAlmostEqual(cell.calculate_1_d2([-3, -2, 1]),
                               0.0128229081865688, delta=1e-17)
        self.assertAlmostEqual(cell.calculate_d([3, 4, 5]),
                               7.7319559244298, delta=1e-13)


class TestAngles(unittest.TestCase):
    def test_dihedral_special_cases(self):
        a = Position(random(), random(), random())
        # not sure what it should be in such undefined cases
        #self.assertTrue(isnan(calculate_dihedral(a, a, a, a)))
        self.assertEqual(calculate_dihedral(a, a, a, a), 0.0)
        # Special cases from scitbx tst_math.py
        # atan2 is guaranteed to give exact values (I think)
        p000 = Position(0, 0, 0)
        p100 = Position(1, 0, 0)
        p010 = Position(0, 1, 0)
        def xy_dihedral(last_point):
            return calculate_dihedral(p100, p000, p010, last_point)
        self.assertEqual(xy_dihedral(Position(1, 1, 0)), 0.0)
        self.assertEqual(xy_dihedral(Position(-1, 1, 0)), pi)
        p01_ = Position(0, 1, -1)
        self.assertEqual(xy_dihedral(p01_), pi/2)
        p01_.z = 1
        self.assertEqual(xy_dihedral(p01_), -pi/2)

    def test_dihedral(self):
        # based on from https://stackoverflow.com/questions/20305272/
        p0 = Position(24.969, 13.428, 30.692)  # N
        p1 = Position(24.044, 12.661, 29.808)  # CA
        p2 = Position(22.785, 13.482, 29.543)  # C
        p3 = Position(21.951, 13.670, 30.431)  # O
        p4 = Position(23.672, 11.328, 30.466)  # CB
        p5 = Position(22.881, 10.326, 29.620)  # CG
        p6 = Position(23.691,  9.935, 28.389)  # CD1
        p7 = Position(22.557,  9.096, 30.459)  # CD2
        def check_dihedral(a, b, c, d, angle):
            deg = calculate_dihedral(a, b, c, d) * 180 / pi
            self.assertAlmostEqual(deg, angle, places=4)
        check_dihedral(p0, p1, p2, p3, -71.21515)
        check_dihedral(p0, p1, p4, p5, -171.94319)
        check_dihedral(p1, p4, p5, p6, 60.82226)
        check_dihedral(p1, p4, p5, p7, -177.63641)

if __name__ == '__main__':
    unittest.main()
