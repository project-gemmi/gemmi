#!/usr/bin/env python3

import gemmi  # type: ignore
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / 'build' / 'py'))


def resolve_acedrg_tables_from_ccp4() -> str:
    ccp4 = os.environ.get('CCP4')
    if not ccp4:
        raise unittest.SkipTest('CCP4 environment variable is not set')
    tables_dir = Path(ccp4) / 'share' / 'acedrg' / 'tables'
    if not tables_dir.is_dir():
        raise unittest.SkipTest(
            f'AceDRG tables directory is not available: {tables_dir}')
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
        self.assertAlmostEqual(
            positions['A1'].dist(
                positions['A2']), 1.50, places=5)
        self.assertAlmostEqual(
            positions['A2'].dist(
                positions['A3']), 1.40, places=5)
        self.assertAlmostEqual(
            positions['A3'].dist(
                positions['A4']), 1.30, places=5)
        self.assertAlmostEqual(
            math.degrees(
                gemmi.calculate_angle(
                    positions['A1'],
                    positions['A2'],
                    positions['A3'])),
            112.0,
            places=4)
        self.assertAlmostEqual(
            math.degrees(
                gemmi.calculate_angle(
                    positions['A2'],
                    positions['A3'],
                    positions['A4'])),
            121.0,
            places=4)
        self.assertAlmostEqual(
            math.degrees(
                gemmi.calculate_dihedral(
                    positions['A1'],
                    positions['A2'],
                    positions['A3'],
                    positions['A4'])),
            60.0,
            places=4)

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

        for other, dist in [('A1', 1.43), ('A2', 1.47),
                            ('A3', 1.53), ('H1', 1.00)]:
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

    def test_generate_chemcomp_xyz_preserves_abp_sugar_branch(self):
        path = REPO_ROOT / 'ccd' / 'gemmi' / 'a' / 'ABP.cif'
        cc = gemmi.make_chemcomp_from_block(
            gemmi.cif.read(str(path)).sole_block())
        for atom in cc.atoms:
            atom.xyz = gemmi.Position(float('nan'), float('nan'), float('nan'))
        placed = gemmi.generate_chemcomp_xyz_from_restraints(cc)
        self.assertEqual(placed, len(cc.atoms))
        pos = {atom.id: atom.xyz for atom in cc.atoms}
        self.assertLess(pos["C2'"].dist(pos["O2'"]), 2.0)
        self.assertLess(pos["C3'"].dist(pos["C2'"]), 2.0)
        self.assertLess(pos["C2'"].dist(pos["C1'"]), 2.0)

    def test_generate_chemcomp_xyz_preserves_alb_amide_bridges(self):
        path = REPO_ROOT / 'ccd' / 'gemmi' / 'a' / 'ALB.cif'
        cc = gemmi.make_chemcomp_from_block(
            gemmi.cif.read(str(path)).sole_block())
        for atom in cc.atoms:
            atom.xyz = gemmi.Position(float('nan'), float('nan'), float('nan'))
        placed = gemmi.generate_chemcomp_xyz_from_restraints(cc)
        self.assertEqual(placed, len(cc.atoms))
        pos = {atom.id: atom.xyz for atom in cc.atoms}
        self.assertAlmostEqual(pos['C21'].dist(pos['N6']), 1.3376, delta=0.15)
        self.assertAlmostEqual(pos['C22'].dist(pos['N6']), 1.4549, delta=0.15)
        self.assertAlmostEqual(pos['C24'].dist(pos['N7']), 1.3376, delta=0.15)
        self.assertAlmostEqual(pos['C25'].dist(pos['N7']), 1.4553, delta=0.15)

    def test_generate_chemcomp_xyz_preserves_a6e_planar_articulation(self):
        path = REPO_ROOT / 'ccd' / 'gemmi' / 'a' / 'A6E.cif'
        cc = gemmi.make_chemcomp_from_block(
            gemmi.cif.read(str(path)).sole_block())
        for atom in cc.atoms:
            atom.xyz = gemmi.Position(float('nan'), float('nan'), float('nan'))
        placed = gemmi.generate_chemcomp_xyz_from_restraints(cc)
        self.assertEqual(placed, len(cc.atoms))
        pos = {atom.id: atom.xyz for atom in cc.atoms}
        self.assertAlmostEqual(pos['C6'].dist(pos['N']), 1.3537, delta=0.10)
        self.assertLess(pos['C6'].dist(pos['N']), 2.0)

    def test_generate_chemcomp_xyz_preserves_a6x_planar_bridge(self):
        path = REPO_ROOT / 'ccd' / 'gemmi' / 'a' / 'A6X.cif'
        cc = gemmi.make_chemcomp_from_block(
            gemmi.cif.read(str(path)).sole_block())
        for atom in cc.atoms:
            atom.xyz = gemmi.Position(float('nan'), float('nan'), float('nan'))
        placed = gemmi.generate_chemcomp_xyz_from_restraints(cc)
        self.assertEqual(placed, len(cc.atoms))
        pos = {atom.id: atom.xyz for atom in cc.atoms}
        self.assertAlmostEqual(pos['C23'].dist(pos['N10']), 1.4608, delta=0.10)
        self.assertAlmostEqual(pos['N10'].dist(pos['C01']), 1.3211, delta=0.10)

    def test_generate_chemcomp_xyz_preserves_a1v_planar_chain(self):
        path = REPO_ROOT / 'ccd' / 'gemmi' / 'a' / 'A1V.cif'
        cc = gemmi.make_chemcomp_from_block(
            gemmi.cif.read(str(path)).sole_block())
        for atom in cc.atoms:
            atom.xyz = gemmi.Position(float('nan'), float('nan'), float('nan'))
        placed = gemmi.generate_chemcomp_xyz_from_restraints(cc)
        self.assertEqual(placed, len(cc.atoms))
        pos = {atom.id: atom.xyz for atom in cc.atoms}
        self.assertAlmostEqual(pos['C10'].dist(pos['N09']), 1.3858, delta=0.10)
        self.assertAlmostEqual(pos['N09'].dist(pos['C08']), 1.3072, delta=0.10)
        self.assertAlmostEqual(pos['C11'].dist(pos['S12']), 1.7059, delta=0.12)
        self.assertAlmostEqual(pos['S12'].dist(pos['C08']), 1.7235, delta=0.12)

    def test_generate_chemcomp_xyz_preserves_a4w_planar_linker(self):
        path = REPO_ROOT / 'ccd' / 'gemmi' / 'a' / 'A4W.cif'
        cc = gemmi.make_chemcomp_from_block(
            gemmi.cif.read(str(path)).sole_block())
        for atom in cc.atoms:
            atom.xyz = gemmi.Position(float('nan'), float('nan'), float('nan'))
        placed = gemmi.generate_chemcomp_xyz_from_restraints(cc)
        self.assertEqual(placed, len(cc.atoms))
        pos = {atom.id: atom.xyz for atom in cc.atoms}
        self.assertAlmostEqual(pos['C03'].dist(pos['N04']), 1.3951, delta=0.05)
        self.assertAlmostEqual(pos['N04'].dist(pos['C05']), 1.3524, delta=0.05)

    def test_generate_chemcomp_xyz_rescues_a7y_fragment_attachment(self):
        path = REPO_ROOT / 'ccd' / 'gemmi' / 'a' / 'A7Y.cif'
        cc = gemmi.make_chemcomp_from_block(
            gemmi.cif.read(str(path)).sole_block())
        for atom in cc.atoms:
            atom.xyz = gemmi.Position(float('nan'), float('nan'), float('nan'))
        placed = gemmi.generate_chemcomp_xyz_from_restraints(cc)
        self.assertEqual(placed, len(cc.atoms))
        pos = {atom.id: atom.xyz for atom in cc.atoms}
        self.assertAlmostEqual(pos['C41'].dist(pos['C38']), 1.4865, delta=0.05)
        self.assertAlmostEqual(pos['C38'].dist(pos['C39']), 1.4179, delta=0.05)
        self.assertAlmostEqual(pos['C38'].dist(pos['N37']), 1.3222, delta=0.05)

    def test_generate_chemcomp_xyz_rescues_awi_fragment_attachment(self):
        path = REPO_ROOT / 'ccd' / 'gemmi' / 'a' / 'AWI.cif'
        cc = gemmi.make_chemcomp_from_block(
            gemmi.cif.read(str(path)).sole_block())
        for atom in cc.atoms:
            atom.xyz = gemmi.Position(float('nan'), float('nan'), float('nan'))
        placed = gemmi.generate_chemcomp_xyz_from_restraints(cc)
        self.assertEqual(placed, len(cc.atoms))
        pos = {atom.id: atom.xyz for atom in cc.atoms}
        self.assertAlmostEqual(pos['C1'].dist(pos['O22']), 1.3841, delta=0.05)
        self.assertAlmostEqual(pos['C1'].dist(pos['C6']), 1.3874, delta=0.05)
        self.assertAlmostEqual(pos['O22'].dist(pos['C23']), 1.4239, delta=0.05)

    def test_generate_chemcomp_xyz_rescues_a4u_macrocycle_closure(self):
        path = REPO_ROOT / 'ccd' / 'gemmi' / 'a' / 'A4U.cif'
        cc = gemmi.make_chemcomp_from_block(
            gemmi.cif.read(str(path)).sole_block())
        for atom in cc.atoms:
            atom.xyz = gemmi.Position(float('nan'), float('nan'), float('nan'))
        placed = gemmi.generate_chemcomp_xyz_from_restraints(cc)
        self.assertEqual(placed, len(cc.atoms))
        pos = {atom.id: atom.xyz for atom in cc.atoms}
        self.assertAlmostEqual(pos['C47'].dist(pos['C48']), 1.1979, delta=0.10)
        self.assertAlmostEqual(pos['C48'].dist(pos['C03']), 1.4347, delta=0.10)
        self.assertAlmostEqual(pos['C46'].dist(pos['C47']), 1.4229, delta=0.10)

    def test_generate_chemcomp_xyz_places_methyl_hydrogens_symmetrically(self):
        cc = gemmi.ChemComp()
        cc.name = 'TMET'
        cc.group = gemmi.ChemComp.Group.NonPolymer
        for atom_id, el in [('C0', 'C'), ('C1', 'C'),
                            ('H1', 'H'), ('H2', 'H'), ('H3', 'H')]:
            atom = gemmi.ChemComp.Atom()
            atom.id = atom_id
            atom.el = gemmi.Element(el)
            atom.chem_type = el
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

        add_bond('C0', 'C1', 1.53)
        for h in ('H1', 'H2', 'H3'):
            add_bond('C1', h, 1.09)
            add_angle('C0', 'C1', h, 109.5)
        add_angle('H1', 'C1', 'H2', 109.5)
        add_angle('H1', 'C1', 'H3', 109.5)
        add_angle('H2', 'C1', 'H3', 109.5)

        placed = gemmi.generate_chemcomp_xyz_from_restraints(cc)
        self.assertEqual(placed, 5)
        pos = {atom.id: atom.xyz for atom in cc.atoms}
        self.assertAlmostEqual(
            math.degrees(
                gemmi.calculate_angle(
                    pos['H1'],
                    pos['C1'],
                    pos['H2'])),
            109.5,
            places=3)
        self.assertAlmostEqual(
            math.degrees(
                gemmi.calculate_angle(
                    pos['H1'],
                    pos['C1'],
                    pos['H3'])),
            109.5,
            places=3)
        self.assertAlmostEqual(
            math.degrees(
                gemmi.calculate_angle(
                    pos['H2'],
                    pos['C1'],
                    pos['H3'])),
            109.5,
            places=3)

    def test_drg_only_xyz_on_prepared_file(self):
        gemmi_bin = REPO_ROOT / 'build' / 'gemmi'
        source = REPO_ROOT / 'ccd' / 'gemmi' / 'a' / 'ALA.cif'
        self.assertTrue(
            gemmi_bin.is_file(),
            msg=f'missing gemmi executable: {gemmi_bin}')
        self.assertTrue(source.is_file(), msg=f'missing prepared CIF: {source}')

        with tempfile.TemporaryDirectory(prefix='drg_only_xyz_') as tmpdir:
            out = Path(tmpdir) / 'ALA_only_xyz.cif'
            proc = subprocess.run(
                [str(gemmi_bin), 'drg', '--only-xyz', str(source), str(out)],
                text=True, capture_output=True, check=False)
            proc_msg = (
                'drg --only-xyz failed:\n'
                f'STDOUT:\n{proc.stdout}\n'
                f'STDERR:\n{proc.stderr}'
            )
            self.assertEqual(
                proc.returncode, 0,
                msg=proc_msg)

            validate = subprocess.run(
                [
                    str(gemmi_bin), 'validate', '-m', '-v',
                    '--z-score=2', str(out),
                ],
                text=True, capture_output=True, check=False)
            combined = (validate.stdout + validate.stderr).strip()
            validate_msg = (
                'validate -m failed:\n'
                f'STDOUT:\n{validate.stdout}\n'
                f'STDERR:\n{validate.stderr}'
            )
            self.assertEqual(
                validate.returncode, 0,
                msg=validate_msg)
            self.assertNotIn('[atom.xyz]', combined, msg=combined)
            self.assertIn('OK', combined, msg=combined)

    def test_prepare_xyz_regression_harness(self):
        raise unittest.SkipTest(
            'prepare-mode raw embedding is still exploratory; '
            'use ace_xyz_regression.py manually')


if __name__ == '__main__':
    unittest.main()
