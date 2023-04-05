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
        self.assertEqual(selstr('CA[C]:'), '//*//CA[C]:')
        self.assertEqual(selstr('CA[C]:A'), '//*//CA[C]:A')
        self.assertEqual(selstr('A;b>20'), '//A//;b>20')
        self.assertEqual(selstr('/1;q>0.5'), '/1/*//;q>0.5')
        self.assertEqual(selstr('[C];b<20'), '//*//[C];b<20')
        self.assertEqual(selstr(':;b=20'), '//*//:;b=20')
        self.assertEqual(selstr(':B'), '//*//:B')
        self.assertEqual(selstr(':B;q=0'), '//*//:B;q=0')
        self.assertEqual(selstr('A//CA'), '//A//CA')
        self.assertEqual(selstr('(ALA)'), '//*/(ALA)/')
        self.assertEqual(selstr('(ALA,GLY)'), '//*/(ALA,GLY)/')
        self.assertEqual(selstr('15-55'), '//*/15.-55./')
        self.assertEqual(selstr('15C-55B'), '//*/15.C-55.B/')
        self.assertEqual(selstr('15.C-55.B'), '//*/15.C-55.B/')
        self.assertEqual(selstr('///'), '////')
        self.assertEqual(selstr('[Cu]'), '//*//[Cu]')
        self.assertEqual(selstr('[Mg,O,X]'), '//*//[X,O,Mg]')
        self.assertEqual(selstr('[!Xe]'), '//*//[!Xe]')
        self.assertEqual(selstr('[!H,D]'), '//*//[!H,D]')
        self.assertEqual(selstr(';q=0'), '//*//;q=0')
        self.assertEqual(selstr('(ALA);polymer;b<30;q>0'),
                         '//*/(ALA)/;polymer;b<30;q>0')
        self.assertEqual(selstr('[!H,D];!polymer,solvent'),
                         '//*//[!H,D];!polymer,solvent')
        self.assertEqual(selstr('A/33.-120.A/[C]'), '//A/33.-120.A/[C]')
        self.assertEqual(selstr('B/7'), '//B/7./')

if __name__ == '__main__':
    unittest.main()
