#!/usr/bin/env python

import os
import unittest
import gemmi


class TestRealCif(unittest.TestCase):
    def test_cod_sic(self):
        path = os.path.join(os.path.dirname(__file__), '1011031.cif')
        sic = gemmi.read_atomic_structure(path)
        self.assertEqual(sic.name, '1011031')
        self.assertEqual(sic.cell.a, 4.358)
        self.assertEqual(sic.cell.alpha, 90.0)

if __name__ == '__main__':
    unittest.main()
