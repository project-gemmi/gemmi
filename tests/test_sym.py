#!/usr/bin/env python

import unittest
import random
from gemmi import sym
try:
    from cctbx import sgtbx
    print('(w/ sgtbx)')
except ImportError:
    sgtbx = None

TDEN = sym.Op.TDEN

CANONICAL_SINGLES = {
 "x"     : [ 1, 0, 0, 0],
 "z"     : [ 0, 0, 1, 0],
 "-y"    : [ 0,-1, 0, 0],
 "-z"    : [ 0, 0,-1, 0],
 "x-y"   : [ 1,-1, 0, 0],
 "-x+y"  : [-1, 1, 0, 0],
 "x+1/2" : [ 1, 0, 0, TDEN//2],
 "y+1/4" : [ 0, 1, 0, TDEN//4],
 "z+3/4" : [ 0, 0, 1, TDEN*3//4],
 "z+1/3" : [ 0, 0, 1, TDEN*1//3],
 "z+1/6" : [ 0, 0, 1, TDEN//6],
 "z+2/3" : [ 0, 0, 1, TDEN*2//3],
 "z+5/6" : [ 0, 0, 1, TDEN*5//6],
 "-x+1/4": [-1, 0, 0, TDEN//4],
 "-y+1/2": [ 0,-1, 0, TDEN//2],
 "-y+3/4": [ 0,-1, 0, TDEN*3//4],
 "-z+1/3": [ 0, 0,-1, TDEN//3],
 "-z+1/6": [ 0, 0,-1, TDEN//6],
 "-z+2/3": [ 0, 0,-1, TDEN*2//3],
 "-z+5/6": [ 0, 0,-1, TDEN*5//6],
}

OTHER_SINGLES = {
 # order and letter case may vary
 "Y-x"   : [-1, 1, 0,  0],
 "-X"    : [-1, 0, 0,  0],
 "-1/2+Y": [ 0, 1, 0, -TDEN//2],
 # we want to handle non-crystallographic translations
 "x+3"   : [ 1, 0, 0,  TDEN*3],
 "1+Y"   : [ 0, 1, 0,  TDEN],
 "-2+Y"  : [ 0, 1, 0, -TDEN*2],
 "-z-5/6": [ 0, 0,-1, -TDEN*5//6],
}


class TestSymmetry(unittest.TestCase):
    def test_parse_triplet_part(self):
        for single, row in CANONICAL_SINGLES.items():
            calculated = sym.parse_triplet_part(single)
            self.assertEqual(calculated, row)
        for single, row in OTHER_SINGLES.items():
            calculated = sym.parse_triplet_part(single)
            self.assertEqual(calculated, row)

    def test_make_triplet_part(self):
        self.assertEqual(sym.make_triplet_part(0, 0, 0, 1),
                         '1/%d' % sym.Op.TDEN)
        for single, row in CANONICAL_SINGLES.items():
            calculated = sym.make_triplet_part(*row)
            self.assertEqual(calculated, single)

    def test_triplet_roundtrip(self):
        singles = list(CANONICAL_SINGLES.keys())
        for i in range(4):
            items = [random.choice(singles) for j in range(3)]
            triplet = ','.join(items)
            op = sym.parse_triplet(triplet)
            self.assertEqual(op.triplet(), triplet)

    def test_combine(self):
        a = sym.Op('x+1/3,z,-y')
        self.assertEqual(sym.combine(a, a).triplet(), 'x+2/3,-y,-z')
        self.assertEqual('x,-y,z' * sym.Op('-x,-y,z'), '-x,y,z')
        a = sym.Op('-y+1/4,x+3/4,z+1/4')
        b = sym.Op('-x+1/2,y,-z')
        self.assertEqual((a * b).triplet(), '-y+1/4,-x+1/4,-z+1/4')
        c = '-y,-z,-x'
        self.assertNotEqual(b * c, c * b)
        self.assertEqual((a * c).triplet(), 'z+1/4,-y+3/4,-x+1/4')
        self.assertEqual(b * c, sym.Op('y+1/2,-z,x'))
        self.assertEqual(c * b, '-y,z,x+1/2')

    def test_invert(self):
        for xyz in ['-y,-x,-z+1/4', 'y,-x,z+3/4', 'y,x,-z', 'y+1/2,x,-z+1/3']:
            op = sym.Op(xyz)
            self.assertEqual(op * op.inverted(), 'x,y,z')
            self.assertEqual(op.inverted().inverted(), op)
        op = sym.Op("-y+z,x+z,-x+y+z") # det=3
        self.assertRaises(RuntimeError, op.inverted)

    def test_generators_from_hall(self):
        # first test on example matrices from
        # http://cci.lbl.gov/sginfo/hall_symbols.html
        self.assertEqual(sym.generators_from_hall('p -2xc').sym_ops,
                         ['x,y,z', '-x,y,z+1/2'])
        self.assertEqual(sym.generators_from_hall('p 3*').sym_ops,
                         ['x,y,z', 'z,x,y'])
        self.assertEqual(sym.generators_from_hall('p 4vw').sym_ops,
                         ['x,y,z', '-y,x+1/4,z+1/4'])
        self.assertEqual(sym.generators_from_hall('p 61 2 (0 0 -1)').sym_ops,
                         ['x,y,z', 'x-y,x,z+1/6', '-y,-x,-z+5/6'])
        # then on examples from the 530 settings
        self.assertEqual(sym.generators_from_hall('P -2 -2').sym_ops,
                         ['x,y,z', 'x,y,-z', '-x,y,z'])
        # the same operations in different notation
        a = sym.generators_from_hall('P 3*')
        b = sym.generators_from_hall('R 3 (-y+z,x+z,-x+y+z)')
        self.assertEqual(a.sym_ops, b.sym_ops)
        self.assertEqual(a.cen_ops, b.cen_ops)

    def compare_hall_symops_with_sgtbx(self, hall):
        cctbx_sg = sgtbx.space_group(hall)
        cctbx_ops = set(m.as_xyz() for m in cctbx_sg.all_ops(mod=1))
        gemmi_sg = sym.symops_from_hall(hall)
        self.assertEqual(len(gemmi_sg.sym_ops), cctbx_sg.order_p())
        self.assertEqual(len(gemmi_sg.cen_ops), cctbx_sg.n_ltr())
        gemmi_ops = set(m.triplet() for m in gemmi_sg)
        self.assertEqual(cctbx_ops, gemmi_ops)

    def test_with_sgtbx(self):
        if sgtbx is None:
            return
        for s in sgtbx.space_group_symbol_iterator():
            self.compare_hall_symops_with_sgtbx(s.hall())
        self.compare_hall_symops_with_sgtbx('C -4 -2b')

    def test_find_spacegroup(self):
        self.assertEqual(sym.SpaceGroup('P21212').hm, 'P 21 21 2')
        self.assertEqual(sym.find_spacegroup_by_name('P21').hm, 'P 1 21 1')
        self.assertEqual(sym.find_spacegroup_by_name('P 2').hm, 'P 1 2 1')
        def check_xhm(name, xhm):
            self.assertEqual(sym.SpaceGroup(name).xhm(), xhm)
        check_xhm('R 3 2', 'R 3 2:H')
        check_xhm('R 3 2', 'R 3 2:H')
        check_xhm('R32:H', 'R 3 2:H')
        check_xhm('R 3 2:R', 'R 3 2:R')
        check_xhm('P6', 'P 6')
        check_xhm('P 6', 'P 6')
        check_xhm('P65', 'P 65')
        check_xhm('I1211', 'I 1 21 1')
        check_xhm('Aem2', 'A b m 2')
        check_xhm('C c c e', 'C c c a:1')
        self.assertRaises(ValueError, sym.SpaceGroup, 'i1')
        check_xhm('i2', 'I 1 2 1')
        self.assertEqual(sym.find_spacegroup_by_number(5).hm, 'C 1 2 1')
        self.assertEqual(sym.SpaceGroup(4005).hm, 'I 1 2 1')
        self.assertIsNone(sym.find_spacegroup_by_name('abc'))

    def test_operations(self):
        gops = sym.symops_from_hall('-P 2a 2ac (z,x,y)')
        self.assertEqual(set(sym.SpaceGroup('Pbaa').operations()), set(gops))
        self.assertEqual(sym.find_spacegroup_by_ops(gops).hm, 'P b a a')

if __name__ == '__main__':
    unittest.main()
