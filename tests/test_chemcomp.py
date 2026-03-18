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

# tests for make_structure_from_chemcomp_block()
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
        self.assertEqual(len(cif_st), 2)  # two models: Ideal and Example

        which = gemmi.ChemCompModel.Xyz
        cif_st = gemmi.make_structure_from_chemcomp_block(cif_block, which)
        self.assertEqual(len(cif_st), 0)

        which = gemmi.ChemCompModel.Xyz | gemmi.ChemCompModel.Ideal
        cif_st = gemmi.make_structure_from_chemcomp_block(cif_block, which)
        self.assertEqual(len(cif_st), 1)

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

    def test_assign_chemcomp_ccp4_types_without_tables(self):
        path = os.path.join(os.path.dirname(__file__), 'ccd', 'ALA.cif')
        block = gemmi.cif.read(path).sole_block()
        cc = gemmi.make_chemcomp_from_block(block)

        gemmi.assign_chemcomp_ccp4_types(cc)

        atom_types = {atom.id: atom.chem_type for atom in cc.atoms}
        self.assertEqual(atom_types['N'], 'N32')
        self.assertEqual(atom_types['CA'], 'CH1')
        self.assertEqual(atom_types['CB'], 'CH3')
        self.assertEqual(atom_types['C'], 'C')
        self.assertEqual(atom_types['O'], 'O')

class TestSmarts(unittest.TestCase):
    def test_smarts_benzene(self):
        cif_text = """
data_comp_BEN
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
BEN C1 C 0
BEN C2 C 0
BEN C3 C 0
BEN C4 C 0
BEN C5 C 0
BEN C6 C 0
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
BEN C1 C2 arom
BEN C2 C3 arom
BEN C3 C4 arom
BEN C4 C5 arom
BEN C5 C6 arom
BEN C6 C1 arom
"""
        doc = gemmi.cif.read_string(cif_text)
        cc = gemmi.make_chemcomp_from_block(doc.sole_block())

        # [C] should NOT match aromatic carbons
        self.assertEqual(len(cc.match_smarts("[C]")), 0)

        # [c] should match aromatic carbons
        matches = cc.match_smarts("[c]")
        self.assertEqual(len(matches), 6)

        # cc (two aromatic carbons with implicit bond) should match
        self.assertGreater(len(cc.match_smarts("cc")), 0)

    def test_smarts_ester(self):
        cif_text = """
data_comp_MAE
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
C1 C
C2 C
O1 O
O2 O
C3 C
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
C1 C2 sing
C2 O1 doub
C2 O2 sing
O2 C3 sing
"""
        doc = gemmi.cif.read_string(cif_text)
        cc = gemmi.make_chemcomp_from_block(doc.sole_block())

        # Match carbonyl C=O
        matches = cc.match_smarts("C=O")
        self.assertEqual(len(matches), 1)
        self.assertEqual(cc.atoms[matches[0][0]].id, "C2")
        self.assertEqual(cc.atoms[matches[0][1]].id, "O1")

        # Match ester group CC(=O)OC
        matches_ester = cc.match_smarts("CC(=O)OC")
        self.assertEqual(len(matches_ester), 1)

        # Match oxygen with degree 2
        matches_o2 = cc.match_smarts("[OX2]")
        self.assertEqual(len(matches_o2), 1)
        self.assertEqual(cc.atoms[matches_o2[0][0]].id, "O2")


class TestAceRings(unittest.TestCase):
    def test_find_ace_rings_benzene(self):
        cif_text = """
data_comp_BEN
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
C1 C
C2 C
C3 C
C4 C
C5 C
C6 C
H1 H
H2 H
H3 H
H4 H
H5 H
H6 H
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
C1 C2 arom
C2 C3 arom
C3 C4 arom
C4 C5 arom
C5 C6 arom
C6 C1 arom
C1 H1 sing
C2 H2 sing
C3 H3 sing
C4 H4 sing
C5 H5 sing
C6 H6 sing
"""
        cc = gemmi.make_chemcomp_from_block(gemmi.cif.read_string(cif_text).sole_block())
        rings = gemmi.find_ace_rings(cc)

        self.assertEqual(len(rings), 1)
        self.assertEqual(len(rings[0].atoms), 6)
        self.assertTrue(rings[0].is_aromatic)
        self.assertTrue(rings[0].is_aromatic_permissive)


if __name__ == '__main__':
    unittest.main()
