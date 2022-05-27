#!/usr/bin/env python

import unittest
import copy
import math
import gemmi

# Example cell parameters and space groups with twinning operators from cctbx.
# The following cctbx script was used to get twinning operators.
#    xs = cctbx.crystal.symmetry(unit_cell=uctbx.unit_cell(params),
#                                space_group_symbol=sg_symbol,
#                                correct_rhombohedral_setting_if_necessary=True)
#    ms = cctbx.miller.set(xs, cctbx.array_family.flex.miller_index())
#    ma = cctbx.miller.array(ms, data=cctbx.array_family.flex.double(()))
#    ops = [op.operator.as_xyz() for op in
#           mmtbx.scaling.twin_analyses.twin_laws(ma).operators]
#    print((params, sg_symbol, ops))
TWINNING_DATA = [
    ([96.6, 96.6, 91.1, 90.0, 90.0, 120.0], 'P 1',
     ['x-y,x,z', 'y,-x+y,z', '-y,x-y,z', '-x+y,-x,z', '-x,-y,z', 'x-y,-y,-z',
      '-x,-x+y,-z', 'y,x,-z', '-y,-x,-z', '-x+y,y,-z', 'x,x-y,-z']),
    ([52.304, 52.255, 120.1, 89.11, 88.24, 60.12], 'P 1',
     ['x+y,-x,z', '-y,x+y,z', 'y,-x-y,z', '-x-y,x,z', '-x,x+y,-z', 'x+y,-y,-z',
      '-x,-y,z', 'y,x,-z', '-y,-x,-z', 'x,-x-y,-z', '-x-y,y,-z']),
    ([107.61, 80.6, 110.1, 90.0, 88.59, 90.0], 'I 1 2 1',
     ['z,y,-x', '-z,-y,-x', 'x,-y,-z']),
    ([152.089, 30.767, 105.94, 90.0, 133.15, 90.0], 'C 1 2 1',
     ['x-z,y,2*x-z', 'x-z,-y,-z', '-x,-y,-2*x+z']),
    ([283.5, 401.8, 284.0, 90.0, 89.4, 90.0], 'P 1 21 1',
     ['z,y,-x', 'x,-y,-z', '-z,-y,-x']),
    ([114.6, 117.9, 114.0, 90.0, 90.0, 90.0], 'P 21 21 21',
     ['y,x,-z', '-x,-z,-y', 'z,-y,x', 'z,x,y', 'y,z,x']),
    ([17.75, 17.76, 42.77, 90.0, 90.0, 120.05], 'P 1 1 21',
     ['-x+y,-x,z', '-y,x-y,z', 'y,x,-z', 'x-y,-y,-z', '-x,-x+y,-z']),
    ([186.1, 185.9, 186.9, 90.0, 90.0, 90.0], 'F 2 2 2',
     ['z,-y,x', '-x,z,y', '-y,-x,-z', '-y,z,-x', '-z,-x,y']),
    ([83.122, 83.4, 83.917, 90.0, 90.0, 90.0], 'I 21 21 21',
     ['-x,-z,-y', 'z,-y,x', '-y,-x,-z', '-z,x,-y', 'y,-z,-x']),
    ([60.7, 60.7, 97.0, 90.0, 90.0, 120.0], 'P 32 2 1',
     ['-x,-y,z']),
    ([91.93, 168.02, 137.77, 90.0, 90.0, 90.0], 'C 2 2 2',
     ['1/2*x-3/2*y,-1/2*x-1/2*y,-z', '1/2*x+3/2*y,1/2*x-1/2*y,-z']),
    ([133.02, 133.02, 133.02, 90.0, 90.0, 90.0], 'I 21 3',
     ['-z,-y,-x']),
    ([174.22, 53.12, 75.17, 90.0, 115.29, 90.0], 'C 1 2 1',
     ['-x,-y,-2*x+z']),
    ([157.13, 86.69, 79.87, 90.0, 90.0, 90.0], 'C 2 2 21',
     ['-1/2*x-1/2*y,-3/2*x+1/2*y,-z', '-1/2*x+1/2*y,3/2*x+1/2*y,-z']),
    ([75.58, 81.06, 90.53, 86.23, 81.86, 63.92], 'P 1',
     ['-x-y,y,-z']),
    ([93.33, 93.33, 157.14, 90.0, 90.0, 90.0], 'P 43',
     ['x,-y,-z']),
    ([78.2, 78.2, 87.4, 90.0, 90.0, 120.0], 'P 64',
     ['x-y,-y,-z']),
    ([77.17, 77.17, 90.46, 90.0, 90.0, 120.0], 'H 3',
     ['-2/3*x-1/3*y+2/3*z,-1/3*x-2/3*y-2/3*z,2/3*x-2/3*y+1/3*z',
      '-x+1/3*y+2/3*z,-1/3*y+4/3*z,2/3*y+1/3*z',
      '-1/3*x-4/3*z,1/3*x-y-2/3*z,-2/3*x+1/3*z',
      '-x+2/3*y-2/3*z,1/3*y-4/3*z,-2/3*y-1/3*z',
      '1/3*x+4/3*z,2/3*x-y+2/3*z,2/3*x-1/3*z',
      '-1/3*x-2/3*y-2/3*z,-2/3*x-1/3*y+2/3*z,-2/3*x+2/3*y-1/3*z',
      'x-y,-y,-z']),
    ([162.9, 162.9, 100.7, 90.0, 90.0, 120.0], 'H 3 2',
     ['-1/3*x+2/3*z,1/3*x-y+1/3*z,4/3*x+1/3*z',
      '-2/3*x-1/3*y-1/3*z,-1/3*x-2/3*y+1/3*z,-4/3*x+4/3*y+1/3*z',
      '-x+1/3*y-1/3*z,-1/3*y-2/3*z,-4/3*y+1/3*z']),
    ([112.0, 112.0, 402.7, 90.0, 90.0, 90.0], 'P 43 21 2', []),
    ([95.829, 46.431, 65.186, 90.0, 115.39, 90.0], 'C 1 2 1', []),
    ([86.89, 86.89, 99.01, 90.0, 90.0, 120.0], 'H 3 2', []),
    ([72.2, 146.8, 88.9, 90.0, 90.0, 90.0], 'C 2 2 21', []),
    ([76.7, 72.6, 57.0, 111.5, 82.8, 62.6], 'P 1', []),
    ([106.06, 106.06, 294.13, 90.0, 90.0, 120.0], 'P 61 2 2', []),
    # small molecule examples
    ([4.1537, 4.1537, 6.862, 90.0, 90.0, 120.0], 'P -3 m 1', ['-x,-y,z']),
    ([4.358, 4.358, 4.358, 90.0, 90.0, 90.0], 'F -4 3 m', []),
]

class TestTwinning(unittest.TestCase):
    def test_lattice_2fold_ops(self):
        # from Andrey's example
        cell = gemmi.UnitCell(102.053, 46.612, 74.904, 95.00, 71.29, 95.00)
        op_scores = gemmi.find_lattice_2fold_ops(cell, max_obliq=5)
        self.assertEqual(len(op_scores), 1)
        op, score = op_scores[0]
        self.assertEqual(op.triplet(), 'x,-y,-x-z')
        self.assertAlmostEqual(score, math.degrees(0.07946439), delta=1e-6)

    def test_lattice_symmetry_r(self):
        cell = gemmi.UnitCell(30.0, 30.0, 219.0, 90.0, 90.0, 90.0)
        gops = gemmi.find_lattice_symmetry_r(cell, 5)
        sg = gemmi.find_spacegroup_by_ops(gops)
        self.assertTrue(sg and sg.hm == 'P 4 2 2')

    def test_lattice_symmetry(self):
        cell = gemmi.UnitCell(119.353, 119.47, 184.789, 89.76, 89.9, 90.22)
        # original space group: 'I 41 2 2'
        gops = gemmi.find_lattice_symmetry(cell, centring='I', max_obliq=5)
        lattice_sg = gemmi.find_spacegroup_by_ops(gops)
        self.assertTrue(lattice_sg and lattice_sg.hm == 'I 4 2 2')

    def test_potential_twinning(self):
        for (cell_params, sg_symbol, cctbx_ops) in TWINNING_DATA:
            sg = gemmi.SpaceGroup(sg_symbol)
            cell = gemmi.UnitCell(*cell_params)
            twin_ops = gemmi.find_twin_laws(cell, sg, 3.0, all_ops=False)
            self.assertEqual(len(twin_ops), len(cctbx_ops))
            # We should get the same cosets wrt. the point group,
            # but the coset representatives can differ.
            # To compare them, we combine them with the point group ops.
            gops1 = sg.operations()
            gops2 = copy.deepcopy(gops1)
            pg_symops = gops1.derive_symmorphic().sym_ops
            gops1.sym_ops = pg_symops + [gemmi.Op(o) for o in cctbx_ops]
            gops1.add_missing_elements()
            gops2.sym_ops = pg_symops + twin_ops
            gops2.add_missing_elements()
            self.assertEqual(set(gops1.sym_ops), set(gops2.sym_ops))


if __name__ == '__main__':
    unittest.main()
