#!/usr/bin/env python

import unittest
from math import pi, isnan
from random import random
from gemmi import Position, calculate_dihedral


class TestAngles(unittest.TestCase):
    def test_dihedral_special_cases(self):
        a = Position(random(), random(), random())
        self.assertTrue(isnan(calculate_dihedral(a, a, a, a)))
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
        p0 = Position(24.969, 13.428, 30.692) # N
        p1 = Position(24.044, 12.661, 29.808) # CA
        p2 = Position(22.785, 13.482, 29.543) # C
        p3 = Position(21.951, 13.670, 30.431) # O
        p4 = Position(23.672, 11.328, 30.466) # CB
        p5 = Position(22.881, 10.326, 29.620) # CG
        p6 = Position(23.691,  9.935, 28.389) # CD1
        p7 = Position(22.557,  9.096, 30.459) # CD2
        def check_dihedral(a, b, c, d, angle):
            deg = calculate_dihedral(a, b, c, d) * 180 / pi
            self.assertAlmostEqual(deg, angle, places=4)
        check_dihedral(p0, p1, p2, p3, -71.21515)
        check_dihedral(p0, p1, p4, p5, -171.94319)
        check_dihedral(p1, p4, p5, p6, 60.82226)
        check_dihedral(p1, p4, p5, p7, -177.63641)

if __name__ == '__main__':
    unittest.main()
