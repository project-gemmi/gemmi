#!/usr/bin/env python

import os
import unittest
import gemmi

SO2_FROM_MONOMER = """\
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1                      
HETATM    1  S   SO3     1      -5.979   0.717  17.353  1.00 20.00           S  
HETATM    2  O1  SO3     1      -5.035   1.876  17.325  1.00 20.00           O  
HETATM    3  O2  SO3     1      -7.003   1.053  16.315  1.00 20.00           O1-
HETATM    4  O3  SO3     1      -5.199  -0.407  16.748  1.00 20.00           O1-
"""  # noqa: W291 - trailing whitespace

def full_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)

# tests for gemmi/chemcomp_xyz.hpp
class TestChemCompXyz(unittest.TestCase):
    def test_reading_monomer_SO3_coordinates(self):
        path = full_path('SO3.cif')
        block = gemmi.cif.read(path)[-1]
        st = gemmi.make_structure_from_chemcomp_block(block)
        out = st.make_pdb_string(gemmi.PdbWriteOptions(minimal=True))
        self.assertEqual(out.splitlines(), SO2_FROM_MONOMER.splitlines())

    # comparing HEM.cif from PDB CCD with HEM.pdb from PDBe
    def test_reading_HEM(self):
        cif_path = full_path('HEM.cif')
        cif_block = gemmi.cif.read(cif_path).sole_block()
        cif_st = gemmi.make_structure_from_chemcomp_block(cif_block)
        self.assertEqual(len(cif_st), 2)
        # we compare not-ideal model only
        del cif_st['example_xyz']
        # PDBe files have residue number 0 and ATOM instead of HETATM
        residue = cif_st[0][0][0]
        residue.seqid.num = 0
        residue.het_flag = 'A'
        for atom in residue:
            atom.b_iso = 20
        minimal_opt = gemmi.PdbWriteOptions(minimal=True)
        cif_out = cif_st.make_pdb_string(minimal_opt)
        pdb_path = full_path('HEM.pdb')
        pdb_st = gemmi.read_structure(pdb_path)
        pdb_out = pdb_st.make_pdb_string(minimal_opt)
        self.assertEqual(cif_out.splitlines(), pdb_out.splitlines())

    # HEN.cif from CCD does not provide ideal coordinates
    def test_reading_HEN(self):
        path = full_path('HEN.cif')
        block = gemmi.cif.read(path).sole_block()
        st = gemmi.make_structure_from_chemcomp_block(block)
        self.assertEqual(len(st), 1)

# tests for gemmi/chemcomp.hpp
class TestChemComp(unittest.TestCase):
    def test_hem_ccd(self):
        path = os.path.join(os.path.dirname(__file__), 'HEM.cif')
        block = gemmi.cif.read(path).sole_block()
        cc = gemmi.make_chemcomp_from_block(block)
        self.assertEqual(len(cc.atoms), 75)
        self.assertEqual(len(cc.rt.bonds), 82)
        bond = cc.rt.get_bond("CGA", "O1A")
        self.assertEqual(bond.type, gemmi.BondType.Double)
        self.assertEqual(bond.aromatic, False)
        bond = cc.rt.get_bond("NA", "C1A")
        self.assertEqual(bond.type, gemmi.BondType.Single)
        self.assertEqual(bond.aromatic, True)

    def test_hen_ccd(self):
        path = os.path.join(os.path.dirname(__file__), 'HEN.cif')
        block = gemmi.cif.read(path).sole_block()
        cc = gemmi.make_chemcomp_from_block(block)
        self.assertEqual(len(cc.atoms), 45)
        self.assertEqual(len(cc.rt.bonds), 45)

    def test_so2_ccp4(self):
        path = os.path.join(os.path.dirname(__file__), 'SO3.cif')
        block = gemmi.cif.read(path)[-1]
        cc = gemmi.make_chemcomp_from_block(block)
        self.assertEqual(len(cc.atoms), 4)
        self.assertEqual(len(cc.rt.bonds), 3)
        self.assertAlmostEqual(cc.rt.get_bond("S", "O2").value, 1.496)
        self.assertAlmostEqual(cc.rt.get_bond("O3", "S").esd, 0.019)
        with self.assertRaises(RuntimeError):
            cc.rt.get_bond("O2", "O3")

    def test_shortest_path(self):
        path = full_path('HEM.cif')
        block = gemmi.cif.read(path)[-1]
        cc = gemmi.make_chemcomp_from_block(block)
        A = gemmi.Restraints.AtomId
        result = cc.rt.find_shortest_path(A('FE'), A('CMD'), [])
        self.assertEqual([aid.atom for aid in result],
                         ['FE', 'ND', 'C1D', 'C2D', 'CMD'])
        result = cc.rt.find_shortest_path(A('CMD'), A('FE'), [])
        self.assertEqual([aid.atom for aid in result],
                         ['CMD', 'C2D', 'C1D', 'ND', 'FE'])
        result = cc.rt.find_shortest_path(A('FE'), A('CMD'), [A('C1D')])
        self.assertEqual([aid.atom for aid in result],
                         ['FE', 'ND', 'C4D', 'C3D', 'C2D', 'CMD'])

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
