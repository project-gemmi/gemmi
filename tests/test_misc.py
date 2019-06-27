#!/usr/bin/env python

import unittest
import os
import gemmi

class TestMisc(unittest.TestCase):
    def test_pdb_code(self):
        self.assertTrue(gemmi.is_pdb_code('9ALA'))
        self.assertFalse(gemmi.is_pdb_code('BAD1'))
        os.environ['PDB_DIR'] = '/data'
        cif_4xyz = '/data/structures/divided/mmCIF/xy/4xyz.cif.gz'
        self.assertEqual(gemmi.expand_pdb_code_to_path('4XYZ', 'M'), cif_4xyz)
        self.assertEqual(gemmi.expand_if_pdb_code('4XYZ'), cif_4xyz)
        self.assertEqual(gemmi.expand_if_pdb_code('4XYZ', filetype='P'),
                         '/data/structures/divided/pdb/xy/pdb4xyz.ent.gz')

if __name__ == '__main__':
    unittest.main()
