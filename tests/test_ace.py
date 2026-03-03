#!/usr/bin/env python

import contextlib
import os
import subprocess
import tempfile
import unittest
import gemmi

TOP_DIR = os.path.join(os.path.dirname(__file__), "..")


@contextlib.contextmanager
def temp_env(var, value):
    had = var in os.environ
    old = os.environ.get(var)
    if value is None:
        os.environ.pop(var, None)
    else:
        os.environ[var] = value
    try:
        yield
    finally:
        if had:
            os.environ[var] = old
        else:
            os.environ.pop(var, None)


class TestAcePrepareChemComp(unittest.TestCase):
    @staticmethod
    def make_atom(atom_id, el, chem_type, charge=0.0):
        atom = gemmi.ChemComp.Atom()
        atom.id = atom_id
        atom.el = gemmi.Element(el)
        atom.charge = charge
        atom.chem_type = chem_type
        return atom

    @staticmethod
    def make_bond(a1, a2, bond_type, value, esd, aromatic=False):
        b = gemmi.Restraints.Bond()
        b.id1 = gemmi.Restraints.AtomId(a1)
        b.id2 = gemmi.Restraints.AtomId(a2)
        b.type = bond_type
        b.aromatic = aromatic
        b.value = value
        b.esd = esd
        b.value_nucleus = float('nan')
        b.esd_nucleus = float('nan')
        return b

    @staticmethod
    def prepare(cc):
        tables = gemmi.AcedrgTables()
        gemmi.prepare_chemcomp(cc, tables, {}, True)

    def test_prepare_chemcomp_cleans_invalid_restraint_refs(self):
        cc = gemmi.ChemComp()
        cc.name = 'TINV'
        cc.group = gemmi.ChemComp.Group.NonPolymer
        cc.atoms.append(self.make_atom('C1', 'C', 'C'))
        cc.atoms.append(self.make_atom('C2', 'C', 'C'))
        cc.atoms.append(self.make_atom('C3', 'C', 'C'))

        cc.rt.bonds.append(self.make_bond('C1', 'C2', gemmi.BondType.Single, 1.50, 0.02))

        a = gemmi.Restraints.Angle()
        a.id1 = gemmi.Restraints.AtomId('C1')
        a.id2 = gemmi.Restraints.AtomId('C2')
        a.id3 = gemmi.Restraints.AtomId('X')
        a.value = 120.0
        a.esd = 3.0
        cc.rt.angles.append(a)

        t = gemmi.Restraints.Torsion()
        t.label = 'auto'
        t.id1 = gemmi.Restraints.AtomId('C1')
        t.id2 = gemmi.Restraints.AtomId('C2')
        t.id3 = gemmi.Restraints.AtomId('X')
        t.id4 = gemmi.Restraints.AtomId('C3')
        t.value = 180.0
        t.esd = 10.0
        t.period = 3
        cc.rt.torsions.append(t)

        ch = gemmi.Restraints.Chirality()
        ch.id_ctr = gemmi.Restraints.AtomId('X')
        ch.id1 = gemmi.Restraints.AtomId('C1')
        ch.id2 = gemmi.Restraints.AtomId('C2')
        ch.id3 = gemmi.Restraints.AtomId('C3')
        ch.sign = gemmi.ChiralityType.Positive
        cc.rt.chirs.append(ch)

        p = gemmi.Restraints.Plane()
        p.label = 'p'
        p.ids = [
            gemmi.Restraints.AtomId('C1'),
            gemmi.Restraints.AtomId('X'),
            gemmi.Restraints.AtomId('C2'),
        ]
        p.esd = 0.02
        cc.rt.planes.append(p)

        self.prepare(cc)
        existing = {a.id for a in cc.atoms}
        for b in cc.rt.bonds:
            self.assertIn(b.id1.atom, existing)
            self.assertIn(b.id2.atom, existing)
        for ang in cc.rt.angles:
            self.assertIn(ang.id1.atom, existing)
            self.assertIn(ang.id2.atom, existing)
            self.assertIn(ang.id3.atom, existing)
        for tor in cc.rt.torsions:
            self.assertIn(tor.id1.atom, existing)
            self.assertIn(tor.id2.atom, existing)
            self.assertIn(tor.id3.atom, existing)
            self.assertIn(tor.id4.atom, existing)
        for chir in cc.rt.chirs:
            self.assertIn(chir.id_ctr.atom, existing)
            self.assertIn(chir.id1.atom, existing)
            self.assertIn(chir.id2.atom, existing)
            self.assertIn(chir.id3.atom, existing)
        for plane in cc.rt.planes:
            for atom_id in plane.ids:
                self.assertIn(atom_id.atom, existing)

    def test_prepare_chemcomp_deduplicates_restraints_centrally(self):
        cc = gemmi.ChemComp()
        cc.name = 'TDED'
        cc.group = gemmi.ChemComp.Group.NonPolymer
        for atom_id in ('C1', 'C2', 'C3', 'C4'):
            cc.atoms.append(self.make_atom(atom_id, 'C', 'C'))

        cc.rt.bonds.append(self.make_bond('C1', 'C2', gemmi.BondType.Single, 1.50, 0.02))
        cc.rt.bonds.append(self.make_bond('C2', 'C1', gemmi.BondType.Single, 1.50, 0.02))
        cc.rt.bonds.append(self.make_bond('C2', 'C3', gemmi.BondType.Single, 1.50, 0.02))
        cc.rt.bonds.append(self.make_bond('C3', 'C4', gemmi.BondType.Single, 1.50, 0.02))

        a1 = gemmi.Restraints.Angle()
        a1.id1 = gemmi.Restraints.AtomId('C1')
        a1.id2 = gemmi.Restraints.AtomId('C2')
        a1.id3 = gemmi.Restraints.AtomId('C3')
        a1.value = 120.0
        a1.esd = 3.0
        cc.rt.angles.append(a1)
        a2 = gemmi.Restraints.Angle()
        a2.id1 = gemmi.Restraints.AtomId('C3')
        a2.id2 = gemmi.Restraints.AtomId('C2')
        a2.id3 = gemmi.Restraints.AtomId('C1')
        a2.value = 120.0
        a2.esd = 3.0
        cc.rt.angles.append(a2)

        t1 = gemmi.Restraints.Torsion()
        t1.label = 'dup'
        t1.id1 = gemmi.Restraints.AtomId('C1')
        t1.id2 = gemmi.Restraints.AtomId('C2')
        t1.id3 = gemmi.Restraints.AtomId('C3')
        t1.id4 = gemmi.Restraints.AtomId('C4')
        t1.value = 180.0
        t1.esd = 10.0
        t1.period = 3
        cc.rt.torsions.append(t1)
        t2 = gemmi.Restraints.Torsion()
        t2.label = 'dup'
        t2.id1 = gemmi.Restraints.AtomId('C4')
        t2.id2 = gemmi.Restraints.AtomId('C3')
        t2.id3 = gemmi.Restraints.AtomId('C2')
        t2.id4 = gemmi.Restraints.AtomId('C1')
        t2.value = 180.0
        t2.esd = 10.0
        t2.period = 3
        cc.rt.torsions.append(t2)

        ch1 = gemmi.Restraints.Chirality()
        ch1.id_ctr = gemmi.Restraints.AtomId('C2')
        ch1.id1 = gemmi.Restraints.AtomId('C1')
        ch1.id2 = gemmi.Restraints.AtomId('C3')
        ch1.id3 = gemmi.Restraints.AtomId('C4')
        ch1.sign = gemmi.ChiralityType.Positive
        cc.rt.chirs.append(ch1)
        ch2 = gemmi.Restraints.Chirality()
        ch2.id_ctr = gemmi.Restraints.AtomId('C2')
        ch2.id1 = gemmi.Restraints.AtomId('C3')
        ch2.id2 = gemmi.Restraints.AtomId('C1')
        ch2.id3 = gemmi.Restraints.AtomId('C4')
        ch2.sign = gemmi.ChiralityType.Positive
        cc.rt.chirs.append(ch2)

        p1 = gemmi.Restraints.Plane()
        p1.label = 'pl'
        p1.ids = [
            gemmi.Restraints.AtomId('C1'),
            gemmi.Restraints.AtomId('C2'),
            gemmi.Restraints.AtomId('C3'),
        ]
        p1.esd = 0.02
        cc.rt.planes.append(p1)
        p2 = gemmi.Restraints.Plane()
        p2.label = 'pl'
        p2.ids = [
            gemmi.Restraints.AtomId('C3'),
            gemmi.Restraints.AtomId('C2'),
            gemmi.Restraints.AtomId('C1'),
        ]
        p2.esd = 0.02
        cc.rt.planes.append(p2)

        self.prepare(cc)
        self.assertEqual(len(cc.rt.bonds), 3)
        self.assertEqual(len(cc.rt.angles), 1)
        self.assertEqual(len(cc.rt.torsions), 1)
        self.assertEqual(len(cc.rt.chirs), 1)
        self.assertEqual(len(cc.rt.planes), 1)

    def test_prepare_chemcomp_normalizes_nitro_groups(self):
        cc = gemmi.ChemComp()
        cc.name = 'TNIT'
        cc.group = gemmi.ChemComp.Group.NonPolymer
        cc.atoms.append(self.make_atom('N1', 'N', 'N'))
        cc.atoms.append(self.make_atom('C1', 'C', 'C'))
        cc.atoms.append(self.make_atom('O1', 'O', 'O'))
        cc.atoms.append(self.make_atom('O2', 'O', 'O'))

        cc.rt.bonds.append(self.make_bond('N1', 'C1', gemmi.BondType.Single, 1.45, 0.02))
        cc.rt.bonds.append(self.make_bond('N1', 'O1', gemmi.BondType.Double, 1.22, 0.02))
        cc.rt.bonds.append(self.make_bond('N1', 'O2', gemmi.BondType.Double, 1.22, 0.02))

        self.prepare(cc)
        atoms = {a.id: a for a in cc.atoms}
        self.assertAlmostEqual(atoms['N1'].charge, 1.0)
        o_neg = int(atoms['O1'].charge < -0.5) + int(atoms['O2'].charge < -0.5)
        self.assertEqual(o_neg, 1)

        n_o_single = 0
        n_o_double = 0
        for b in cc.rt.bonds:
            n_o = (
                (b.id1.atom == 'N1' and b.id2.atom in ('O1', 'O2'))
                or (b.id2.atom == 'N1' and b.id1.atom in ('O1', 'O2'))
            )
            if not n_o:
                continue
            if b.type == gemmi.BondType.Single:
                n_o_single += 1
            if b.type in (gemmi.BondType.Double, gemmi.BondType.Deloc):
                n_o_double += 1
        self.assertEqual(n_o_single, 1)
        self.assertEqual(n_o_double, 1)

    def test_prepare_chemcomp_protonates_terminal_amine(self):
        cc = gemmi.ChemComp()
        cc.name = 'TAMN'
        cc.group = gemmi.ChemComp.Group.NonPolymer
        cc.atoms.append(self.make_atom('N', 'N', 'N'))
        cc.atoms.append(self.make_atom('H', 'H', 'H'))
        cc.atoms.append(self.make_atom('H2', 'H', 'H'))
        cc.atoms.append(self.make_atom('CA', 'C', 'C'))
        cc.atoms.append(self.make_atom('C', 'C', 'C'))
        cc.atoms.append(self.make_atom('O', 'O', 'O'))
        cc.atoms.append(self.make_atom('OXT', 'O', 'O'))
        cc.atoms.append(self.make_atom('HXT', 'H', 'H'))

        cc.rt.bonds.append(self.make_bond('N', 'H', gemmi.BondType.Single, 1.0, 0.02))
        cc.rt.bonds.append(self.make_bond('N', 'H2', gemmi.BondType.Single, 1.0, 0.02))
        cc.rt.bonds.append(self.make_bond('N', 'CA', gemmi.BondType.Single, 1.46, 0.02))
        cc.rt.bonds.append(self.make_bond('CA', 'C', gemmi.BondType.Single, 1.53, 0.02))
        cc.rt.bonds.append(self.make_bond('C', 'O', gemmi.BondType.Double, 1.24, 0.02))
        cc.rt.bonds.append(self.make_bond('C', 'OXT', gemmi.BondType.Single, 1.32, 0.02))
        cc.rt.bonds.append(self.make_bond('OXT', 'HXT', gemmi.BondType.Single, 1.0, 0.02))

        self.prepare(cc)
        atoms = {a.id: a for a in cc.atoms}
        self.assertAlmostEqual(atoms['N'].charge, 1.0)
        self.assertIn('H3', atoms)
        self.assertNotIn('HXT', atoms)

    def test_prepare_chemcomp_applies_metal_neighbor_charge_correction(self):
        cc = gemmi.ChemComp()
        cc.name = 'TMET'
        cc.group = gemmi.ChemComp.Group.NonPolymer
        cc.atoms.append(self.make_atom('O1', 'O', 'O'))
        cc.atoms.append(self.make_atom('C1', 'C', 'C'))
        cc.atoms.append(self.make_atom('ZN1', 'Zn', 'ZN'))

        cc.rt.bonds.append(self.make_bond('O1', 'C1', gemmi.BondType.Single, 1.34, 0.02))
        cc.rt.bonds.append(self.make_bond('O1', 'ZN1', gemmi.BondType.Single, 2.0, 0.04))

        self.prepare(cc)
        atoms = {a.id: a for a in cc.atoms}
        self.assertAlmostEqual(atoms['O1'].charge, -1.0)

    def test_prepare_chemcomp_strict_fails_on_duplicate_plane_atom_ids(self):
        cc = gemmi.ChemComp()
        cc.name = 'TPLN'
        cc.group = gemmi.ChemComp.Group.NonPolymer
        cc.atoms.append(self.make_atom('C1', 'C', 'C'))
        cc.atoms.append(self.make_atom('C2', 'C', 'C'))
        cc.rt.bonds.append(self.make_bond('C1', 'C2', gemmi.BondType.Single, 1.50, 0.02))

        p = gemmi.Restraints.Plane()
        p.label = 'p'
        p.ids = [
            gemmi.Restraints.AtomId('C1'),
            gemmi.Restraints.AtomId('C1'),
            gemmi.Restraints.AtomId('C2'),
        ]
        p.esd = 0.02
        cc.rt.planes.append(p)

        with temp_env('GEMMI_ACE_STRICT', '1'):
            with self.assertRaisesRegex(RuntimeError, 'ACE strict validation failed at post-angle-seed'):
                self.prepare(cc)

    def test_prepare_chemcomp_strict_fails_on_final_nan_bond_value(self):
        cc = gemmi.ChemComp()
        cc.name = 'TNAN'
        cc.group = gemmi.ChemComp.Group.NonPolymer
        cc.atoms.append(self.make_atom('X1', 'X', ''))
        cc.atoms.append(self.make_atom('X2', 'X', ''))
        cc.rt.bonds.append(self.make_bond('X1', 'X2', gemmi.BondType.Single, float('nan'), 0.02))

        with temp_env('GEMMI_ACE_STRICT', '1'):
            with self.assertRaisesRegex(RuntimeError, 'ACE strict validation failed at final'):
                self.prepare(cc)

    def test_prepare_chemcomp_trace_mode_smoke(self):
        cc = gemmi.ChemComp()
        cc.name = 'TTRC'
        cc.group = gemmi.ChemComp.Group.NonPolymer
        cc.atoms.append(self.make_atom('C1', 'C', 'C'))
        cc.atoms.append(self.make_atom('C2', 'C', 'C'))
        cc.rt.bonds.append(self.make_bond('C1', 'C2', gemmi.BondType.Single, 1.50, 0.02))

        with temp_env('GEMMI_ACE_TRACE', '1'):
            self.prepare(cc)
        self.assertEqual(len(cc.rt.bonds), 1)


class TestAceDrgBatch(unittest.TestCase):
    def test_drg_batch_strict_regression_pack(self):
        gemmi_bin = os.path.join(TOP_DIR, 'build', 'gemmi')
        tables_dir = os.path.join(TOP_DIR, 'acedrg', 'tables')
        if not os.path.isfile(gemmi_bin):
            self.skipTest('build/gemmi is not available')
        if not os.path.isdir(tables_dir):
            self.skipTest('acedrg tables directory is not available')

        rel_inputs = [
            'ccd/orig/ALA.cif',
            'ccd/orig/ATP.cif',
            'ccd/orig/CYS.cif',
            'ccd/orig/HEM.cif',
            'ccd/orig/HIS.cif',
            'ccd/orig/SEC.cif',
            'ccd/orig/TRP.cif',
            'ccd/orig/TYR.cif',
        ]
        inputs = [os.path.join(TOP_DIR, p) for p in rel_inputs]
        for path in inputs:
            if not os.path.isfile(path):
                self.skipTest(f'missing test input: {path}')

        with tempfile.TemporaryDirectory(prefix='gemmi-ace-batch-') as out_dir:
            env = os.environ.copy()
            env['GEMMI_ACE_STRICT'] = '1'
            env['ACEDRG_TABLES'] = tables_dir
            subprocess.run(
                [gemmi_bin, 'drg', '--output-dir', out_dir] + inputs,
                cwd=TOP_DIR,
                env=env,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            for input_path in inputs:
                out_path = os.path.join(out_dir, os.path.basename(input_path))
                self.assertTrue(os.path.isfile(out_path), out_path)
                with open(out_path, encoding='utf-8') as f:
                    text = f.read()
                self.assertIn('_chem_comp_bond.atom_id_1', text, out_path)
                self.assertIn('_chem_comp_angle.atom_id_1', text, out_path)


if __name__ == '__main__':
    unittest.main()
