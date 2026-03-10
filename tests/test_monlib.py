#!/usr/bin/env python

import os
import unittest
import gemmi

def full_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)

class TestMonLib(unittest.TestCase):
    def test_path(self):
        path = full_path('')
        monlib = gemmi.read_monomer_lib(path, [])
        self.assertEqual(monlib.path('ALA'), path + "a/ALA.cif")
        self.assertEqual(monlib.path('CON'), path + "c/CON_CON.cif")

    @unittest.skipIf(os.getenv('CLIBD_MON') is None, "$CLIBD_MON not defined.")
    def test_read_monomer_lib(self):
        st = gemmi.read_structure(full_path('4oz7.pdb'))
        resnames = st[0].get_all_residue_names()
        monlib = gemmi.MonLib()
        ok = monlib.read_monomer_lib(os.environ['CLIBD_MON'], resnames)
        self.assertTrue(ok)
        monlib.update_old_atom_names(st)
        topo = gemmi.prepare_topology(st, monlib, model_index=0)
        topo = gemmi.prepare_topology(st, monlib, model_index=0,
                                      h_change=gemmi.HydrogenChange.Shift)
        self.assertIsNotNone(topo)
        # topo = gemmi.prepare_topology(st, monlib, model_index=0,
        #                              h_change=gemmi.HydrogenChange.ReAdd)
        # RuntimeError: Placing of hydrogen bonded to A/22W 6/N failed:
        # Missing angle restraint HN-N-C.


if __name__ == '__main__':
    unittest.main()
