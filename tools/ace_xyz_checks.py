#!/usr/bin/env python3

import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / 'build' / 'py'))

import gemmi  # type: ignore


def resolve_acedrg_tables_from_ccp4() -> str:
    ccp4 = os.environ.get('CCP4')
    if not ccp4:
        raise unittest.SkipTest('CCP4 environment variable is not set')
    tables_dir = Path(ccp4) / 'share' / 'acedrg' / 'tables'
    if not tables_dir.is_dir():
        raise unittest.SkipTest(f'AceDRG tables directory is not available: {tables_dir}')
    return str(tables_dir)


class TestChemCompCoordinateGeneration(unittest.TestCase):
    def test_generate_chemcomp_xyz_from_restraints_exact_chain(self):
        cc = gemmi.ChemComp()
        cc.name = 'TXYZ'
        cc.group = gemmi.ChemComp.Group.NonPolymer
        for atom_id, el, chem_type in [
            ('A1', 'C', 'C'),
            ('A2', 'C', 'C'),
            ('A3', 'C', 'C'),
            ('A4', 'O', 'O'),
        ]:
            atom = gemmi.ChemComp.Atom()
            atom.id = atom_id
            atom.el = gemmi.Element(el)
            atom.chem_type = chem_type
            cc.atoms.append(atom)

        def atom_id(name):
            return gemmi.Restraints.AtomId(name)

        def add_bond(a1, a2, value):
            bond = gemmi.Restraints.Bond()
            bond.id1 = atom_id(a1)
            bond.id2 = atom_id(a2)
            bond.type = gemmi.BondType.Single
            bond.value = value
            bond.esd = 0.02
            cc.rt.bonds.append(bond)

        def add_angle(a1, a2, a3, value):
            angle = gemmi.Restraints.Angle()
            angle.id1 = atom_id(a1)
            angle.id2 = atom_id(a2)
            angle.id3 = atom_id(a3)
            angle.value = value
            angle.esd = 2.0
            cc.rt.angles.append(angle)

        def add_torsion(a1, a2, a3, a4, value, period=3):
            tors = gemmi.Restraints.Torsion()
            tors.label = 'test'
            tors.id1 = atom_id(a1)
            tors.id2 = atom_id(a2)
            tors.id3 = atom_id(a3)
            tors.id4 = atom_id(a4)
            tors.value = value
            tors.esd = 5.0
            tors.period = period
            cc.rt.torsions.append(tors)

        add_bond('A1', 'A2', 1.50)
        add_bond('A2', 'A3', 1.40)
        add_bond('A3', 'A4', 1.30)
        add_angle('A1', 'A2', 'A3', 112.0)
        add_angle('A2', 'A3', 'A4', 121.0)
        add_torsion('A1', 'A2', 'A3', 'A4', 60.0)

        placed = gemmi.generate_chemcomp_xyz_from_restraints(cc)
        self.assertEqual(placed, 4)
        self.assertTrue(cc.has_coordinates)

        positions = {atom.id: atom.xyz for atom in cc.atoms}
        self.assertAlmostEqual(positions['A1'].dist(positions['A2']), 1.50, places=5)
        self.assertAlmostEqual(positions['A2'].dist(positions['A3']), 1.40, places=5)
        self.assertAlmostEqual(positions['A3'].dist(positions['A4']), 1.30, places=5)
        self.assertAlmostEqual(math.degrees(gemmi.calculate_angle(positions['A1'], positions['A2'], positions['A3'])),
                               112.0, places=4)
        self.assertAlmostEqual(math.degrees(gemmi.calculate_angle(positions['A2'], positions['A3'], positions['A4'])),
                               121.0, places=4)
        self.assertAlmostEqual(math.degrees(gemmi.calculate_dihedral(positions['A1'], positions['A2'], positions['A3'], positions['A4'])),
                               60.0, places=4)

    def test_generate_chemcomp_xyz_enforces_plane_restraint(self):
        cc = gemmi.ChemComp()
        cc.name = 'TPLN'
        cc.group = gemmi.ChemComp.Group.NonPolymer
        for atom_id, el, chem_type in [
            ('C1', 'C', 'C'),
            ('C2', 'C', 'C'),
            ('O1', 'O', 'O'),
            ('N1', 'N', 'N'),
        ]:
            atom = gemmi.ChemComp.Atom()
            atom.id = atom_id
            atom.el = gemmi.Element(el)
            atom.chem_type = chem_type
            cc.atoms.append(atom)

        def atom_id(name):
            return gemmi.Restraints.AtomId(name)

        def add_bond(a1, a2, value):
            bond = gemmi.Restraints.Bond()
            bond.id1 = atom_id(a1)
            bond.id2 = atom_id(a2)
            bond.type = gemmi.BondType.Single
            bond.value = value
            bond.esd = 0.02
            cc.rt.bonds.append(bond)

        def add_angle(a1, a2, a3, value):
            angle = gemmi.Restraints.Angle()
            angle.id1 = atom_id(a1)
            angle.id2 = atom_id(a2)
            angle.id3 = atom_id(a3)
            angle.value = value
            angle.esd = 2.0
            cc.rt.angles.append(angle)

        add_bond('C1', 'C2', 1.40)
        add_bond('C2', 'O1', 1.25)
        add_bond('C2', 'N1', 1.35)
        add_angle('C1', 'C2', 'O1', 120.0)
        add_angle('C1', 'C2', 'N1', 120.0)
        add_angle('O1', 'C2', 'N1', 120.0)

        plane = gemmi.Restraints.Plane()
        plane.label = 'pl1'
        plane.ids = [atom_id('C1'), atom_id('C2'), atom_id('O1'), atom_id('N1')]
        plane.esd = 0.02
        cc.rt.planes.append(plane)

        placed = gemmi.generate_chemcomp_xyz_from_restraints(cc)
        self.assertEqual(placed, 4)

        atoms = [a for a in cc.atoms if a.id in ('C1', 'C2', 'O1', 'N1')]
        base = atoms[1].xyz - atoms[0].xyz
        other = atoms[2].xyz - atoms[0].xyz
        normal = base.cross(other)
        self.assertGreater(normal.length(), 1e-4)
        normal = normal.normalized()
        for atom in atoms[3:]:
            dist = abs((atom.xyz - atoms[0].xyz).dot(normal))
            self.assertLess(dist, 1e-6)

    def test_generate_chemcomp_xyz_enforces_chirality(self):
        cc = gemmi.ChemComp()
        cc.name = 'TCHR'
        cc.group = gemmi.ChemComp.Group.NonPolymer
        for atom_id, el, chem_type in [
            ('CTR', 'C', 'C'),
            ('A1', 'O', 'O'),
            ('A2', 'N', 'N'),
            ('A3', 'C', 'C'),
            ('H1', 'H', 'H'),
        ]:
            atom = gemmi.ChemComp.Atom()
            atom.id = atom_id
            atom.el = gemmi.Element(el)
            atom.chem_type = chem_type
            cc.atoms.append(atom)

        def atom_id(name):
            return gemmi.Restraints.AtomId(name)

        def add_bond(a1, a2, value):
            bond = gemmi.Restraints.Bond()
            bond.id1 = atom_id(a1)
            bond.id2 = atom_id(a2)
            bond.type = gemmi.BondType.Single
            bond.value = value
            bond.esd = 0.02
            cc.rt.bonds.append(bond)

        def add_angle(a1, a2, a3, value):
            angle = gemmi.Restraints.Angle()
            angle.id1 = atom_id(a1)
            angle.id2 = atom_id(a2)
            angle.id3 = atom_id(a3)
            angle.value = value
            angle.esd = 2.0
            cc.rt.angles.append(angle)

        for other, dist in [('A1', 1.43), ('A2', 1.47), ('A3', 1.53), ('H1', 1.00)]:
            add_bond('CTR', other, dist)
        add_angle('A1', 'CTR', 'A2', 109.5)
        add_angle('A1', 'CTR', 'A3', 109.5)
        add_angle('A2', 'CTR', 'A3', 109.5)
        add_angle('A1', 'CTR', 'H1', 109.5)
        add_angle('A2', 'CTR', 'H1', 109.5)
        add_angle('A3', 'CTR', 'H1', 109.5)

        chir = gemmi.Restraints.Chirality()
        chir.id_ctr = atom_id('CTR')
        chir.id1 = atom_id('A1')
        chir.id2 = atom_id('A2')
        chir.id3 = atom_id('A3')
        chir.sign = gemmi.ChiralityType.Positive
        cc.rt.chirs.append(chir)

        placed = gemmi.generate_chemcomp_xyz_from_restraints(cc)
        self.assertEqual(placed, 5)

        pos = {atom.id: atom.xyz for atom in cc.atoms}
        v1 = pos['A1'] - pos['CTR']
        v2 = pos['A2'] - pos['CTR']
        v3 = pos['A3'] - pos['CTR']
        vol = v1.dot(v2.cross(v3))
        self.assertGreater(vol, 0.0)

    def test_generate_chemcomp_xyz_regularizes_small_ring(self):
        cc = gemmi.ChemComp()
        cc.name = 'TRNG'
        cc.group = gemmi.ChemComp.Group.NonPolymer
        for i in range(5):
            atom = gemmi.ChemComp.Atom()
            atom.id = f'C{i+1}'
            atom.el = gemmi.Element('C')
            atom.chem_type = 'C'
            cc.atoms.append(atom)

        def atom_id(name):
            return gemmi.Restraints.AtomId(name)

        def add_bond(a1, a2, value):
            bond = gemmi.Restraints.Bond()
            bond.id1 = atom_id(a1)
            bond.id2 = atom_id(a2)
            bond.type = gemmi.BondType.Single
            bond.value = value
            bond.esd = 0.02
            cc.rt.bonds.append(bond)

        def add_angle(a1, a2, a3, value):
            angle = gemmi.Restraints.Angle()
            angle.id1 = atom_id(a1)
            angle.id2 = atom_id(a2)
            angle.id3 = atom_id(a3)
            angle.value = value
            angle.esd = 2.0
            cc.rt.angles.append(angle)

        ring = [f'C{i+1}' for i in range(5)]
        for i in range(5):
            add_bond(ring[i], ring[(i + 1) % 5], 1.45)
        for i in range(5):
            add_angle(ring[i - 1], ring[i], ring[(i + 1) % 5], 108.0)

        placed = gemmi.generate_chemcomp_xyz_from_restraints(cc)
        self.assertEqual(placed, 5)
        pos = {atom.id: atom.xyz for atom in cc.atoms}
        for i in range(5):
            self.assertAlmostEqual(
                pos[ring[i]].dist(pos[ring[(i + 1) % 5]]), 1.45, places=5)

    def test_drg_only_xyz_on_prepared_file(self):
        gemmi_bin = REPO_ROOT / 'build' / 'gemmi'
        source = REPO_ROOT / 'ccd' / 'gemmi' / 'a' / 'ALA.cif'
        self.assertTrue(gemmi_bin.is_file(), msg=f'missing gemmi executable: {gemmi_bin}')
        self.assertTrue(source.is_file(), msg=f'missing prepared CIF: {source}')

        with tempfile.TemporaryDirectory(prefix='drg_only_xyz_') as tmpdir:
            out = Path(tmpdir) / 'ALA_only_xyz.cif'
            proc = subprocess.run(
                [str(gemmi_bin), 'drg', '--only-xyz', str(source), str(out)],
                text=True, capture_output=True, check=False)
            self.assertEqual(
                proc.returncode, 0,
                msg=f'drg --only-xyz failed:\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}')

            validate = subprocess.run(
                [str(gemmi_bin), 'validate', '-m', '-v', '--z-score=2', str(out)],
                text=True, capture_output=True, check=False)
            combined = (validate.stdout + validate.stderr).strip()
            self.assertEqual(
                validate.returncode, 0,
                msg=f'validate -m failed:\nSTDOUT:\n{validate.stdout}\nSTDERR:\n{validate.stderr}')
            self.assertNotIn('[atom.xyz]', combined, msg=combined)
            self.assertIn('OK', combined, msg=combined)

    def test_prepare_xyz_regression_harness(self):
        raise unittest.SkipTest(
            'prepare-mode raw embedding is still exploratory; use ace_xyz_regression.py manually')


if __name__ == '__main__':
    unittest.main()
