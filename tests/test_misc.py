#!/usr/bin/env python

import unittest
import os
import gemmi

class TestMisc(unittest.TestCase):
    def test_pdb_code(self):
        self.assertTrue(gemmi.is_pdb_code('9ALA'))
        self.assertFalse(gemmi.is_pdb_code('BAD1'))
        pdb_dir = os.getenv('PDB_DIR')
        if not pdb_dir:
            return
        cif_4xyz = pdb_dir + '/structures/divided/mmCIF/xy/4xyz.cif.gz'
        self.assertEqual(gemmi.expand_pdb_code_to_path('4XYZ', 'M'), cif_4xyz)
        self.assertEqual(gemmi.expand_if_pdb_code('4XYZ'), cif_4xyz)
        self.assertEqual(gemmi.expand_if_pdb_code('4XYZ', filetype='P'),
                         pdb_dir + '/structures/divided/pdb/xy/pdb4xyz.ent.gz')

    def test_selections(self):
        def selstr(s):
            s2 = gemmi.Selection(s).str()
            s3 = gemmi.Selection(s2).str()
            self.assertEqual(s2, s3)
            return s2
        self.assertEqual(selstr('/1'), '/1/*//')
        self.assertEqual(selstr('CA:*'), '//*//CA')
        self.assertEqual(selstr('A//CA'), '//A//CA')
        self.assertEqual(selstr('///'), '////')
        self.assertEqual(selstr('[Cu]'), '//*//[Cu]')
        self.assertEqual(selstr('[Mg,O,X]'), '//*//[X,O,Mg]')
        self.assertEqual(selstr('[!Xe]'), '//*//[!Xe]')
        self.assertEqual(selstr('[!H,D]'), '//*//[!H,D]')
        self.assertEqual(selstr(';q=0'), '//*//;q=0')
        self.assertEqual(selstr('A/33.-120.A/[C]'), '//A/33.-120.A/[C]')

if __name__ == '__main__':
    unittest.main()
