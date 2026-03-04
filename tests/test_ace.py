#!/usr/bin/env python

import contextlib
import math
import os
import unittest
import gemmi

TOP_DIR = os.path.join(os.path.dirname(__file__), "..")
TESTS_DIR = os.path.dirname(__file__)
CCD_TEST_DIR = os.path.join(TESTS_DIR, 'ccd')


def resolve_acedrg_tables_from_ccp4():
    ccp4 = os.environ.get('CCP4')
    if not ccp4:
        raise unittest.SkipTest('CCP4 environment variable is not set')
    tables_dir = os.path.join(ccp4, 'share', 'acedrg', 'tables')
    if not os.path.isdir(tables_dir):
        raise unittest.SkipTest(f'AceDRG tables directory is not available: {tables_dir}')
    return tables_dir


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


@contextlib.contextmanager
def silence_c_stderr():
    stderr_fd = 2
    saved = os.dup(stderr_fd)
    try:
        with open(os.devnull, 'w', encoding='utf-8') as devnull:
            os.dup2(devnull.fileno(), stderr_fd)
            yield
    finally:
        os.dup2(saved, stderr_fd)
        os.close(saved)


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
            with silence_c_stderr():
                self.prepare(cc)
        self.assertEqual(len(cc.rt.bonds), 1)

    def test_prepare_chemcomp_reuses_deleted_h_name(self):
        cc = gemmi.ChemComp()
        cc.name = 'THREUSE'
        cc.group = gemmi.ChemComp.Group.NonPolymer

        for atom_id, el, chem_type in [
            ('N', 'N', 'N'),
            ('CA', 'C', 'C'),
            ('C', 'C', 'C'),
            ('O', 'O', 'O'),
            ('OXT', 'O', 'O'),
            ('H', 'H', 'H'),
            ('H2', 'H', 'H'),
            ('HXT', 'H', 'H'),
            ('P1', 'P', 'P'),
            ('OP1', 'O', 'O'),
            ('OP2', 'O', 'O'),
            ('OP3', 'O', 'O'),
            ('OP4', 'O', 'O'),
            ('H3', 'H', 'H'),
        ]:
            cc.atoms.append(self.make_atom(atom_id, el, chem_type))

        for a1, a2, btype, value in [
            ('N', 'CA', gemmi.BondType.Single, 1.47),
            ('N', 'H', gemmi.BondType.Single, 1.01),
            ('N', 'H2', gemmi.BondType.Single, 1.01),
            ('CA', 'C', gemmi.BondType.Single, 1.53),
            ('C', 'O', gemmi.BondType.Double, 1.24),
            ('C', 'OXT', gemmi.BondType.Single, 1.31),
            ('OXT', 'HXT', gemmi.BondType.Single, 0.98),
            ('P1', 'OP1', gemmi.BondType.Single, 1.52),
            ('P1', 'OP2', gemmi.BondType.Single, 1.52),
            ('P1', 'OP3', gemmi.BondType.Single, 1.52),
            ('P1', 'OP4', gemmi.BondType.Single, 1.52),
            ('OP1', 'H3', gemmi.BondType.Single, 0.98),
        ]:
            cc.rt.bonds.append(self.make_bond(a1, a2, btype, value, 0.02))

        self.prepare(cc)

        atom_ids = {a.id for a in cc.atoms}
        self.assertNotIn('HXT', atom_ids)
        self.assertIn('H3', atom_ids)
        self.assertNotIn('H4', atom_ids)

        n_h3_bond = any(
            (b.id1.atom == 'N' and b.id2.atom == 'H3') or
            (b.id2.atom == 'N' and b.id1.atom == 'H3')
            for b in cc.rt.bonds
        )
        self.assertTrue(n_h3_bond)


class TestAceDrgBatch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tables = gemmi.AcedrgTables()
        tables_dir = resolve_acedrg_tables_from_ccp4()
        cls.tables.load_tables(tables_dir)

    def test_drg_batch_strict_regression_pack(self):
        rel_inputs = [
            'tests/ccd/ALA.cif',
            'tests/ccd/ATP.cif',
            'tests/ccd/CYS.cif',
            'tests/ccd/HEM.cif',
            'tests/ccd/HIS.cif',
            'tests/ccd/SEC.cif',
            'tests/ccd/TRP.cif',
            'tests/ccd/TYR.cif',
        ]
        inputs = [os.path.join(TOP_DIR, p) for p in rel_inputs]
        for path in inputs:
            if not os.path.isfile(path):
                self.skipTest(f'missing test input: {path}')

        with temp_env('GEMMI_ACE_STRICT', '1'):
            for input_path in inputs:
                doc = gemmi.cif.read(input_path)
                cc = gemmi.make_chemcomp_from_block(doc.sole_block())
                gemmi.prepare_chemcomp(cc, self.tables)
                acedrg_types = self.tables.compute_acedrg_types(cc)

                out_doc = gemmi.cif.Document()
                out_block = out_doc.add_new_block('comp_' + cc.name)
                gemmi.add_chemcomp_to_block(cc, out_block, acedrg_types, False)
                text = out_doc.as_string()

                self.assertIn('_chem_comp_bond.atom_id_1', text, input_path)
                self.assertIn('_chem_comp_angle.atom_id_1', text, input_path)


class TestAcePreparedChemCompInvariants(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tables = gemmi.AcedrgTables()
        cls.tables_dir = resolve_acedrg_tables_from_ccp4()
        cls.tables.load_tables(cls.tables_dir)

    @staticmethod
    def _prepare_from_ccd(comp_id):
        path = os.path.join(CCD_TEST_DIR, f'{comp_id}.cif')
        if not os.path.isfile(path):
            raise unittest.SkipTest(f'missing test input: {path}')
        doc = gemmi.cif.read(path)
        cc = gemmi.make_chemcomp_from_block(doc.sole_block())
        return path, cc

    @staticmethod
    def _make_atom(atom_id, elem, chem_type, charge=0.0):
        atom = gemmi.ChemComp.Atom()
        atom.id = atom_id
        atom.el = gemmi.Element(elem)
        atom.chem_type = chem_type
        atom.charge = charge
        return atom

    @staticmethod
    def _make_nan_bond(a1, a2):
        b = gemmi.Restraints.Bond()
        b.id1 = gemmi.Restraints.AtomId(a1)
        b.id2 = gemmi.Restraints.AtomId(a2)
        b.type = gemmi.BondType.Single
        b.aromatic = False
        b.value = float('nan')
        b.esd = float('nan')
        b.value_nucleus = float('nan')
        b.esd_nucleus = float('nan')
        return b

    @staticmethod
    def _assert_restraint_integrity(cc, source):
        atom_ids = [a.id for a in cc.atoms]
        atom_set = set(atom_ids)
        assert len(atom_ids) == len(atom_set), f'duplicate atom ids in {source}'

        for b in cc.rt.bonds:
            assert b.id1.atom in atom_set and b.id2.atom in atom_set, source
            assert math.isfinite(b.value), source
            assert math.isfinite(b.esd), source
            assert math.isfinite(b.value_nucleus), source
            assert math.isfinite(b.esd_nucleus), source
        for ang in cc.rt.angles:
            assert ang.id1.atom in atom_set and ang.id2.atom in atom_set and ang.id3.atom in atom_set, source
            assert math.isfinite(ang.value), source
            assert math.isfinite(ang.esd), source
        for tor in cc.rt.torsions:
            assert tor.id1.atom in atom_set and tor.id2.atom in atom_set, source
            assert tor.id3.atom in atom_set and tor.id4.atom in atom_set, source
            assert math.isfinite(tor.esd), source
        for chir in cc.rt.chirs:
            assert chir.id_ctr.atom in atom_set and chir.id1.atom in atom_set, source
            assert chir.id2.atom in atom_set and chir.id3.atom in atom_set, source
        for plane in cc.rt.planes:
            plane_atoms = [aid.atom for aid in plane.ids]
            assert len(plane_atoms) == len(set(plane_atoms)), source
            assert len(plane_atoms) >= 3, source
            assert math.isfinite(plane.esd), source
            for atom in plane_atoms:
                assert atom in atom_set, source

    def test_prepare_chemcomp_invariants_on_representative_pack(self):
        for comp_id in ['ALA', 'ATP', 'CYS', 'HEM', 'HIS', 'SEC', 'TRP', 'TYR']:
            path, cc = self._prepare_from_ccd(comp_id)
            gemmi.prepare_chemcomp(cc, self.tables)
            self._assert_restraint_integrity(cc, path)

    def test_prepare_chemcomp_common_motif_charge_sanity(self):
        # Amino-acid zwitterion pattern in free CCD monomer.
        _, ala = self._prepare_from_ccd('ALA')
        gemmi.prepare_chemcomp(ala, self.tables)
        ala_atoms = {a.id: a for a in ala.atoms}
        self.assertLessEqual(ala_atoms['OXT'].charge, -0.5)
        self.assertGreaterEqual(ala_atoms['N'].charge, 0.5)

        # Phosphate-rich ATP should carry multiple negatively charged oxygens.
        _, atp = self._prepare_from_ccd('ATP')
        gemmi.prepare_chemcomp(atp, self.tables)
        atp_neg_o = [
            a.id for a in atp.atoms
            if a.el == gemmi.Element('O') and a.charge <= -0.5
        ]
        self.assertGreaterEqual(len(atp_neg_o), 3)

        # Aspartate sidechain carboxylate oxygen should be deprotonated.
        _, asp = self._prepare_from_ccd('ASP')
        gemmi.prepare_chemcomp(asp, self.tables)
        asp_atoms = {a.id: a for a in asp.atoms}
        self.assertLessEqual(asp_atoms['OD2'].charge, -0.5)

    def test_prepare_chemcomp_proton_hydrogen_element_fallback_sets_nucleus(self):
        cc = gemmi.ChemComp()
        cc.name = 'TPF6H'
        cc.group = gemmi.ChemComp.Group.NonPolymer

        cc.atoms.append(self._make_atom('P1', 'P', 'P'))
        for i in range(1, 7):
            cc.atoms.append(self._make_atom(f'F{i}', 'F', 'F'))
            cc.rt.bonds.append(self._make_nan_bond('P1', f'F{i}'))
        cc.atoms.append(self._make_atom('H1', 'H', 'H'))
        cc.rt.bonds.append(self._make_nan_bond('P1', 'H1'))

        gemmi.prepare_chemcomp(cc, self.tables, {}, True)

        ph = None
        for b in cc.rt.bonds:
            if {b.id1.atom, b.id2.atom} == {'P1', 'H1'}:
                ph = b
                break
        self.assertIsNotNone(ph)
        self.assertAlmostEqual(ph.value, 1.433, places=3)
        self.assertAlmostEqual(ph.value_nucleus, 1.284, places=3)
        self.assertAlmostEqual(ph.esd_nucleus, 0.02, places=3)

    def test_prepare_chemcomp_non_h_bonds_mirror_nucleus_terms(self):
        cc = gemmi.ChemComp()
        cc.name = 'TMIRROR'
        cc.group = gemmi.ChemComp.Group.NonPolymer

        cc.atoms.append(self._make_atom('C1', 'C', 'C'))
        cc.atoms.append(self._make_atom('C2', 'C', 'C'))
        cc.rt.bonds.append(self._make_nan_bond('C1', 'C2'))

        gemmi.prepare_chemcomp(cc, self.tables, {}, True)

        self.assertEqual(len(cc.rt.bonds), 1)
        bb = cc.rt.bonds[0]
        self.assertTrue(math.isfinite(bb.value))
        self.assertTrue(math.isfinite(bb.esd))
        self.assertTrue(math.isfinite(bb.value_nucleus))
        self.assertTrue(math.isfinite(bb.esd_nucleus))
        self.assertAlmostEqual(bb.value_nucleus, bb.value, places=6)
        self.assertAlmostEqual(bb.esd_nucleus, bb.esd, places=6)


class TestChemicalAdjustmentRules(unittest.TestCase):
    """Tests for apply_chemical_adjustments() using CCD examples from chemistry.rst."""

    @staticmethod
    def _load_ccd(comp_id):
        path = os.path.join(CCD_TEST_DIR, f'{comp_id}.cif')
        if not os.path.isfile(path):
            raise unittest.SkipTest(f'missing test input: {path}')
        doc = gemmi.cif.read(path)
        return gemmi.make_chemcomp_from_block(doc.sole_block())

    def _atom_ids(self, cc):
        return {a.id for a in cc.atoms}

    def _charge_of(self, cc, atom_id):
        for a in cc.atoms:
            if a.id == atom_id:
                return a.charge
        self.fail(f'atom {atom_id} not found')

    def _bond_type(self, cc, a1, a2):
        for b in cc.rt.bonds:
            if {b.id1.atom, b.id2.atom} == {a1, a2}:
                return b.type
        return None

    def test_adj_oxoacid_phosphate_atp(self):
        """ATP: phosphate O-H deprotonation removes 4 H, charges 4 O.
        Matches CCP4 monomer library reference (43 atoms, same charges)."""
        cc = self._load_ccd('ATP')
        cc.apply_chemical_adjustments()
        atom_ids = self._atom_ids(cc)
        self.assertEqual(len(cc.atoms), 43)
        for h_id in ['HOA2', 'HOB2', 'HOG2', 'HOG3']:
            self.assertNotIn(h_id, atom_ids)
        for o_id in ['O2A', 'O2B', 'O2G', 'O3G']:
            self.assertAlmostEqual(self._charge_of(cc, o_id), -1.0, places=1)

    def test_adj_oxoacid_sulfate_0sg(self):
        """0SG: sulfate O-H deprotonation removes 2 H, charges 2 O.
        Matches CCP4 monomer library reference (120 atoms, same charges)."""
        cc = self._load_ccd('0SG')
        cc.apply_chemical_adjustments()
        atom_ids = self._atom_ids(cc)
        self.assertEqual(len(cc.atoms), 120)
        for h_id in ['H29', 'H31']:
            self.assertNotIn(h_id, atom_ids)
        for o_id in ['O10', 'O14']:
            self.assertAlmostEqual(self._charge_of(cc, o_id), -1.0, places=1)

    def test_adj_nitro_group_ne5(self):
        """NE5: nitro group N26(=O27)(=O28) → N26(+1)(-O27)(-1)(=O28)."""
        cc = self._load_ccd('NE5')
        cc.apply_chemical_adjustments()
        self.assertAlmostEqual(self._charge_of(cc, 'N26'), 1.0, places=1)
        self.assertAlmostEqual(self._charge_of(cc, 'O27'), -1.0, places=1)
        self.assertEqual(self._bond_type(cc, 'O27', 'N26'), gemmi.BondType.Single)
        self.assertEqual(self._bond_type(cc, 'N26', 'O28'), gemmi.BondType.Double)

    def test_adj_hexafluorophosphate_a9j(self):
        """A9J: PF6 gets an H added (P1-H bond)."""
        cc = self._load_ccd('A9J')
        atoms_before = self._atom_ids(cc)
        cc.apply_chemical_adjustments()
        atoms_after = self._atom_ids(cc)
        added = atoms_after - atoms_before
        self.assertEqual(len(added), 1)
        h_id = added.pop()
        for a in cc.atoms:
            if a.id == h_id:
                self.assertEqual(a.el, gemmi.Element('H'))
        self.assertIsNotNone(self._bond_type(cc, 'P1', h_id))

    def test_adj_carboxy_asp(self):
        """ASP: sidechain carboxylate HD2 removed, OD2 charged; also HXT/OXT.
        Matches CCP4 monomer library reference (15 atoms, same charges)."""
        cc = self._load_ccd('ASP')
        cc.apply_chemical_adjustments()
        atom_ids = self._atom_ids(cc)
        self.assertEqual(len(cc.atoms), 15)
        self.assertNotIn('HD2', atom_ids)
        self.assertNotIn('HXT', atom_ids)
        self.assertAlmostEqual(self._charge_of(cc, 'OD2'), -1.0, places=1)
        self.assertAlmostEqual(self._charge_of(cc, 'OXT'), -1.0, places=1)
        self.assertAlmostEqual(self._charge_of(cc, 'N'), 1.0, places=1)

    def test_adj_terminal_carboxylate_a0g(self):
        """A0G: terminal carboxylate HXT removed, OXT charged.
        Matches CCP4 monomer library reference (11 atoms, same charges)."""
        cc = self._load_ccd('A0G')
        self.assertIn('HXT', self._atom_ids(cc))
        cc.apply_chemical_adjustments()
        self.assertEqual(len(cc.atoms), 11)
        self.assertNotIn('HXT', self._atom_ids(cc))
        self.assertAlmostEqual(self._charge_of(cc, 'OXT'), -1.0, places=1)

    def test_adj_guanidinium_00l(self):
        """00L: guanidinium C=N gets H added, both N charged +1.
        Matches CCP4 monomer library reference (86 atoms, same charges)."""
        cc = self._load_ccd('00L')
        atoms_before = self._atom_ids(cc)
        cc.apply_chemical_adjustments()
        atoms_after = self._atom_ids(cc)
        self.assertEqual(len(cc.atoms), 86)
        added = atoms_after - atoms_before
        self.assertEqual(len(added), 2)
        self.assertAlmostEqual(self._charge_of(cc, 'NH1'), 1.0, places=1)
        self.assertAlmostEqual(self._charge_of(cc, 'N3'), 1.0, places=1)

    def test_adj_amino_ter_amine_00k(self):
        """00K: amino_ter_amine protonates N3, adds H.
        Note: CCP4 reference also protonates NZ (lysine-like terminal amine),
        but our terminal_amine rule does not reach it because the carboxylic
        acid is too far from NZ's alpha carbon."""
        cc = self._load_ccd('00K')
        atoms_before = self._atom_ids(cc)
        cc.apply_chemical_adjustments()
        atoms_after = self._atom_ids(cc)
        added = atoms_after - atoms_before
        self.assertEqual(len(added), 1)
        self.assertAlmostEqual(self._charge_of(cc, 'N3'), 1.0, places=1)

    def test_adj_single_bond_oxide_h1t(self):
        """H1T: single_bond_oxide charges terminal single-bonded oxygens -1.
        H1T is a vanadate compound — 14 terminal oxygens get charge -1.
        CCP4 reference additionally charges 6 bridging oxygens to -2."""
        cc = self._load_ccd('H1T')
        cc.apply_chemical_adjustments()
        self.assertEqual(len(cc.atoms), 27)
        charged_oxygens = ['O01', 'O02', 'O03', 'O04', 'O05', 'O06',
                           'O08', 'O09', 'O10', 'O11', 'O12', 'O13',
                           'O14', 'O19']
        for o_id in charged_oxygens:
            self.assertAlmostEqual(self._charge_of(cc, o_id), -1.0, places=1)

    def test_adj_terminal_amine_lys(self):
        """LYS: terminal_amine protonates N and NZ, terminal_carboxylate
        deprotonates OXT. Matches CCP4 monomer library reference."""
        cc = self._load_ccd('LYS')
        atoms_before = self._atom_ids(cc)
        cc.apply_chemical_adjustments()
        atoms_after = self._atom_ids(cc)
        self.assertEqual(len(cc.atoms), 25)
        self.assertNotIn('HXT', atoms_after)
        self.assertIn('H3', atoms_after - atoms_before)
        self.assertAlmostEqual(self._charge_of(cc, 'N'), 1.0, places=1)
        self.assertAlmostEqual(self._charge_of(cc, 'NZ'), 1.0, places=1)
        self.assertAlmostEqual(self._charge_of(cc, 'OXT'), -1.0, places=1)

    def test_adj_protonated_amide_n_bjs(self):
        """BJS: protonated_amide_n adds H to amide N28, charges N28 +1.
        Matches CCP4 monomer library reference (90 atoms, same charges)."""
        cc = self._load_ccd('BJS')
        atoms_before = self._atom_ids(cc)
        cc.apply_chemical_adjustments()
        atoms_after = self._atom_ids(cc)
        self.assertEqual(len(cc.atoms), 90)
        added = atoms_after - atoms_before
        self.assertEqual(len(added), 1)
        h_id = added.pop()
        self.assertAlmostEqual(self._charge_of(cc, 'N28'), 1.0, places=1)
        self.assertIsNotNone(self._bond_type(cc, 'N28', h_id))


if __name__ == '__main__':
    unittest.main()
