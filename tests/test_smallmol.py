#!/usr/bin/env python

import os
import unittest
import gemmi

def full_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)

# tests for gemmi/smcif.hpp
class TestRealCif(unittest.TestCase):
    def test_cod_sic(self):
        path = full_path('1011031.cif')
        sic = gemmi.read_small_structure(path)
        self.assertEqual(sic.name, '1011031')
        self.assertEqual(sic.cell.a, 4.358)
        self.assertEqual(sic.cell.alpha, 90.0)
        self.assertEqual(len(sic.sites), 2)
        self.assertEqual(len(sic.get_all_unit_cell_sites()), 8)
        self.assertEqual(4 * sic.sites[1].orth(sic.cell).x, sic.cell.a)

    def test_fen4(self):
        path = full_path('2242624.cif')
        small = gemmi.read_small_structure(path)
        types = small.atom_types
        self.assertEqual(len(types), 2)
        self.assertEqual(types[0].symbol, 'Fe')
        self.assertEqual(types[1].element, gemmi.Element('N'))

    def test_perovskite(self):
        small = gemmi.read_small_structure(full_path('4003024.cif'))
        self.assertEqual(small.spacegroup_hm, 'P m -3 m')
        disorder_groups = [site.disorder_group for site in small.sites]
        self.assertEqual(disorder_groups, [0, 1, 0, 2])

    def test_mgi2(self):
        small = gemmi.read_small_structure(full_path('2013551.cif'))
        for site in small.sites:
            u_eq = small.cell.calculate_u_eq(site.aniso)
            self.assertAlmostEqual(u_eq, site.u_iso, delta=0.00012)

if __name__ == '__main__':
    unittest.main()
