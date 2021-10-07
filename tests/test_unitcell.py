#!/usr/bin/env python

import unittest
from math import pi  # , isnan
from random import random
import gemmi
from gemmi import Position, UnitCell

def assert_almost_equal_seq(self, a, b, delta=None):
    for x,y in zip(a, b):
        self.assertAlmostEqual(x, y, delta=delta)

class TestMath(unittest.TestCase):
    def test_Vec3(self):
        a = gemmi.Vec3(-3, 3, 13)
        b = gemmi.Vec3(2, 0, -2)
        self.assertEqual((-a).tolist(), [3, -3, -13])
        self.assertEqual((a+b).tolist(), [-1, 3, 11])
        self.assertEqual((a-b).tolist(), [-5, 3, 15])
        a -= b
        self.assertEqual(a.tolist(), [-5, 3, 15])
        a += b
        self.assertEqual(a.tolist(), [-3, 3, 13])

    def test_SMat33_transformed_by(self):
        tensor = gemmi.SMat33f(random(), random(), random(),
                               random(), random(), random())
        mat = gemmi.Mat33()
        mat.fromlist([[random() for _ in range(3)] for _ in range(3)])
        t1 = tensor.transformed_by(mat).as_mat33().tolist()
        t2 = mat.multiply(tensor.as_mat33()).multiply(mat.transpose()).tolist()
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(t1[i][j], t2[i][j])


class TestUnitCell(unittest.TestCase):
    def test_dummy_cell(self):
        cell = UnitCell()
        self.assertEqual([cell.a, cell.b, cell.c], [1, 1, 1])
        self.assertEqual([cell.alpha, cell.beta, cell.gamma], [90, 90, 90])
        self.assertEqual(cell.volume, 1.0)

    def test_ortho_cell(self):
        cell = UnitCell(25.14, 39.50, 45.07, 90, 90, 90)
        self.assertTrue(cell.orthogonalization_matrix.row_copy(1)
                        .approx(gemmi.Vec3(0, 39.5, 0), 0))
        self.assertTrue(cell.orthogonalization_matrix.column_copy(0)
                        .approx(gemmi.Vec3(25.14, 0, 0), 0))
        pos = Position(5, -6, 7)
        frac = cell.fractionalize(pos)
        self.assertAlmostEqual(frac.x, 0.198886, delta=1e-6)
        self.assertAlmostEqual(frac.y, -0.151899, delta=1e-6)
        self.assertAlmostEqual(frac.z, 0.155314, delta=1e-6)
        pos2 = cell.orthogonalize(frac)
        self.assertAlmostEqual(pos.x, pos2.x, delta=1e-12)
        self.assertAlmostEqual(pos.y, pos2.y, delta=1e-12)
        self.assertAlmostEqual(pos.z, pos2.z, delta=1e-12)
        corner = cell.orthogonalize(gemmi.Fractional(1, 1, 1))
        self.assertAlmostEqual(corner.x, cell.a, delta=1e-12)
        self.assertAlmostEqual(corner.y, cell.b, delta=1e-12)
        self.assertAlmostEqual(corner.z, cell.c, delta=1e-12)
        rec = cell.reciprocal()
        self.assertEqual([rec.alpha, rec.beta, rec.gamma], [90, 90, 90])
        self.assertAlmostEqual(rec.a, 1 / cell.a, delta=1e-17)

    def test_triclinic_cell(self):
        cell = UnitCell(35.996, 41.601, 45.756, 67.40, 66.90, 74.85)
        # this would give syntax error with Python < 3.5
        #o_f = cell.orthogonalization_matrix @ cell.fractionalization_matrix
        o_f = cell.orthogonalization_matrix.multiply(
            cell.fractionalization_matrix)
        self.assertTrue(o_f.approx(gemmi.Mat33(), 1e-15))
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
        #  uc.metrical_matrix()
        cctbx_mm = [1295.712016, 1730.643201, 2093.611536,
                    391.3591013825865, 646.1921687548228, 731.5043620154578]
        mt = cell.metric_tensor()
        assert_almost_equal_seq(self, mt.elements_pdb(), cctbx_mm, delta=1e-12)
        #  uc.reciprocal_metrical_matrix()
        cctbx_rmm = [0.00092792089082916, 0.000689632633981, 0.0006277651322979,
                     -0.000104162588996, -0.000250008091601, -0.000208806754807]
        rmt = cell.reciprocal_metric_tensor()
        assert_almost_equal_seq(self, rmt.elements_pdb(), cctbx_rmm,
                                delta=1e-15)

    def test_change_of_basis(self):
        uc = gemmi.UnitCell(20, 30, 39, 73, 93, 99)
        op = gemmi.Op('y-x/2,-2/3*z+2/3*y,3*x')
        uc2 = uc.changed_basis_backward(op, set_images=False)
        # compare with result from cctbx:
        #  from cctbx import sgtbx, uctbx
        #  u = uctbx.unit_cell((20,30,39, 73,93,99))
        #  op = sgtbx.change_of_basis_op('y-x/2,-2/3*z+2/3*y,3*x').inverse()
        #  print(u.change_basis(cb_op=op).parameters())
        expected = (117.9468784563987, 25.977921933207348, 20.0,
                    130.5, 107.65517573180257, 82.63132106791868)
        assert_almost_equal_seq(self, uc2.parameters, expected)
        uc3 = uc2.changed_basis_forward(op, set_images=True)
        assert_almost_equal_seq(self, uc3.parameters, uc.parameters)

    def test_atom_to_site(self):
        cell = UnitCell(35.996, 41.601, 45.756, 67.40, 66.90, 74.85)
        atom = gemmi.Atom()
        atom.aniso = gemmi.SMat33f(13.1, 20.1, 11.1, -3.5, 5.5, -0.4)
        site = gemmi.SmallStructure.Site(atom, cell)
        # tested against values from cctbx:
        # from cctbx import uctbx, adptbx
        # uc = uctbx.unit_cell((35.996, 41.601, 45.756, 67.40, 66.90, 74.85))
        # aniso = (13.1, 20.1, 11.1, -3.5, 5.5, -0.4)
        # ucif = adptbx.u_cart_as_u_cif(uc, aniso)
        ucif = [11.537759976524049, 19.43436271641311, 11.1,
                -8.078683096677723, 1.4787260755519491, -3.9018967241279157]
        assert_almost_equal_seq(self, site.aniso.elements_pdb(), ucif,
                                delta=1e-6)

    def test_pickling(self):
        try:
            import cPickle as pickle  # Use cPickle on Python 2.7
        except ImportError:
            import pickle

        cell = UnitCell(35.996, 41.601, 45.756, 67.40, 66.90, 74.85)
        pkl_string = pickle.dumps(cell, protocol=pickle.HIGHEST_PROTOCOL)
        result = pickle.loads(pkl_string)
        self.assertTrue(isinstance(result, UnitCell))
        self.assertEqual(cell.parameters, result.parameters)

class TestAngles(unittest.TestCase):
    def test_dihedral_special_cases(self):
        a = Position(random(), random(), random())
        # not sure what it should be in such undefined cases
        #self.assertTrue(isnan(gemmi.calculate_dihedral(a, a, a, a)))
        self.assertEqual(gemmi.calculate_dihedral(a, a, a, a), 0.0)
        # Special cases from scitbx tst_math.py
        # atan2 is guaranteed to give exact values (I think)
        p000 = Position(0, 0, 0)
        p100 = Position(1, 0, 0)
        p010 = Position(0, 1, 0)
        def xy_dihedral(last_point):
            return gemmi.calculate_dihedral(p100, p000, p010, last_point)
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
            deg = gemmi.calculate_dihedral(a, b, c, d) * 180 / pi
            self.assertAlmostEqual(deg, angle, places=4)
        check_dihedral(p0, p1, p2, p3, -71.21515)
        check_dihedral(p0, p1, p4, p5, -171.94319)
        check_dihedral(p1, p4, p5, p6, 60.82226)
        check_dihedral(p1, p4, p5, p7, -177.63641)

class TestGruber(unittest.TestCase):
    def test_reduction(self):
        cell = gemmi.UnitCell(687.9, 687.9, 1933.3, 90.0, 90.0, 90.0)
        sg = gemmi.SpaceGroup('I 4 2 2')
        gv = gemmi.GruberVector(cell, sg)
        gv.niggli_reduce()
        self.assertTrue(gv.is_niggli())
        self.assertTrue(gv.is_buerger())
        self.assertTrue(gv.is_normalized())
        p = cell.a
        q = 1082.134662368783
        t = 108.5325886
        par = (p*p, p*p, q*q, -p*p, -p*p, 0)
        assert_almost_equal_seq(self, gv.parameters, par)
        print(gv.cell_parameters())
        self.assertTrue(gv.cell_parameters(), (p, p, q, t, t, 0))

    def test_near_degenerate(self):
        cell = gemmi.UnitCell(15.53, 91.94, 4.35, 110.326, 7.337, 103.014)
        gv = gemmi.GruberVector(cell, None)
        self.assertFalse(gv.is_normalized())
        gv.normalize()
        self.assertTrue(gv.is_normalized())
        self.assertFalse(gv.is_buerger())
        self.assertFalse(gv.is_niggli())
        n = gv.niggli_reduce(iteration_limit=100)
        self.assertEqual(n, 100)
        self.assertFalse(gv.is_niggli())
        n = gv.niggli_reduce(iteration_limit=100)
        self.assertTrue(n < 100)
        self.assertTrue(gv.is_niggli())
        expected = (2.814597242, 3.077205425, 7.408935896,
                    100.42421409, 94.02885284, 95.07179187)
        assert_almost_equal_seq(self, gv.cell_parameters(), expected)

if __name__ == '__main__':
    unittest.main()
