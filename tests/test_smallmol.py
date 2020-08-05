#!/usr/bin/env python

import os
import unittest
import gemmi

SO2_FROM_MONOMER = """\
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1                      
HETATM    1  S   SO3  -999      -5.979   0.717  17.353  1.00 20.00           S  
HETATM    2  O1  SO3  -999      -5.035   1.876  17.325  1.00 20.00           O  
HETATM    3  O2  SO3  -999      -7.003   1.053  16.315  1.00 20.00           O1-
HETATM    4  O3  SO3  -999      -5.199  -0.407  16.748  1.00 20.00           O1-
"""  # noqa: W291 - trailing whitespace

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

# tests for gemmi/chemcomp_xyz.hpp
class TestChemCompXyz(unittest.TestCase):
    def test_reading_monomer_SO3_coordinates(self):
        path = full_path('SO3.cif')
        block = gemmi.cif.read(path)[-1]
        st = gemmi.make_structure_from_chemcomp_block(block)
        out = st.make_minimal_pdb()
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
        cif_out = cif_st.make_minimal_pdb()
        pdb_path = full_path('HEM.pdb')
        pdb_st = gemmi.read_structure(pdb_path)
        pdb_out = pdb_st.make_minimal_pdb()
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

if __name__ == '__main__':
    unittest.main()
