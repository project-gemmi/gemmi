#!/usr/bin/env python

import gzip
import os
import sys
import tempfile
import unittest
import gemmi

def is_written_to_pdb(line):
    if line[:6] in ['COMPND', 'SOURCE', 'MDLTYP', 'AUTHOR', 'REVDAT', 'JRNL  ',
                    'DBREF ', 'SEQADV', 'HET   ', 'HETNAM', 'FORMUL', 'HELIX ',
                    'SHEET ', 'SITE  ', 'MASTER']:
        return False
    if line[:6] == 'REMARK' and line[6:10] != '   2':
        return False
    return True

def full_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)

class TestMol(unittest.TestCase):
    def test_residue(self):
        res = gemmi.Residue()
        self.assertEqual(res.label_seq, None)
        self.assertEqual(res.seq_num, None)
        res.label_seq = 1
        self.assertEqual(res.label_seq, 1)
        res.seq_num = 2
        self.assertEqual(res.seq_num, 2)
        res.label_seq = None
        self.assertEqual(res.label_seq, None)
        res.seq_num = None
        self.assertEqual(res.seq_num, None)

    def test_read_5i55(self):
        cell = gemmi.read_structure(full_path('5i55.cif')).cell
        self.assertAlmostEqual(cell.a, 29.46)
        self.assertAlmostEqual(cell.b, 10.51)
        self.assertAlmostEqual(cell.c, 29.71)
        self.assertEqual(cell.alpha, 90)
        self.assertAlmostEqual(cell.beta, 111.98)
        self.assertEqual(cell.gamma, 90)

    def test_read_5i55_again(self):
        st = gemmi.read_structure(full_path('5i55.cif'))
        a, b, c, d = st[0]
        ent_a = st.find_entity(a.entity_id)
        self.assertEqual(ent_a.entity_type, gemmi.EntityType.Polymer)
        self.assertEqual(ent_a.polymer_type, gemmi.PolymerType.PeptideL)
        ent_b = st.find_entity(b.entity_id)
        self.assertEqual(ent_b.entity_type, gemmi.EntityType.NonPolymer)
        self.assertEqual(ent_b.polymer_type, gemmi.PolymerType.NA)
        ent_d = st.find_entity(d.entity_id)
        self.assertEqual(ent_d.entity_type, gemmi.EntityType.Water)
        self.assertEqual(ent_d.polymer_type, gemmi.PolymerType.NA)

    def read_1pfe(self, filename):
        st = gemmi.read_structure(full_path(filename))
        self.assertAlmostEqual(st.cell.a, 39.374)
        self.assertEqual(st.cell.gamma, 120)
        self.assertEqual(st.name, '1PFE')
        self.assertEqual(st.sg_hm, 'P 63 2 2')
        self.assertEqual(len(st), 1)
        self.assertEqual(len(st[0]), 7)
        label_name_to_auth_name = {ch.name: ch.auth_name for ch in st[0]}
        self.assertEqual(label_name_to_auth_name,
                         dict(A='A', B='B', C='A', D='B', E='B', F='A', G='B'))
        self.assertEqual(len(st[0]['A'][1]), 1)
        chain_a = st[0]['A']
        self.assertEqual(chain_a[1][0].label_seq, 1)
        self.assertEqual(chain_a['1'][0].seq_num, 1)
        b3 = st[0]['B']['3']
        self.assertEqual(repr(b3), repr(st[0]['B'][3]))
        self.assertEqual(len(b3), 2)
        self.assertEqual(b3[0].name, 'N2C')
        self.assertEqual(b3[-1].name, 'NCY')
        chain_c = st[0]['C']
        self.assertEqual(len(chain_c), 1)
        res_cl = list(chain_c)[0]
        self.assertEqual(res_cl.name, 'CL')
        self.assertEqual(len(res_cl), 1)
        atom_cl = res_cl['CL']
        self.assertAlmostEqual(atom_cl.occ, 0.17)
        self.assertEqual(atom_cl.element.name, 'Cl')

    def test_read_1pfe_cif(self):
        self.read_1pfe('1pfe.cif.gz')

    def test_read_1pfe_json(self):
        self.read_1pfe('1pfe.json')

    def test_read_1orc(self):
        st = gemmi.read_structure(full_path('1orc.pdb'))
        self.assertAlmostEqual(st.cell.a, 34.77)
        self.assertEqual(st.cell.alpha, 90)
        self.assertEqual(len(st.ncs), 0)
        model = st[0]
        self.assertEqual(len(model), 2)
        self.assertTrue(all(res.name == 'HOH' for res in model['A_H']))
        A = model['A']
        self.assertTrue(A['3'])
        self.assertFalse(A[3])
        self.assertEqual([res.seq_num for res in A if res.icode], [56] * 5)
        self.assertEqual(len(A['55']), 1)
        self.assertEqual(len(A['55B']), 0)
        self.assertEqual(len(A['56B']), 1)
        self.assertEqual(A['56'][0].icode, '')
        self.assertEqual(A['56c'][0].icode, 'C')

    def write_back_and_compare(self, path):
        st = gemmi.read_structure(path)
        handle, out_name = tempfile.mkstemp()
        os.close(handle)
        st.write_pdb(out_name)
        with open(out_name) as f:
            out_lines = f.readlines()
        os.remove(out_name)
        return out_lines

    def test_read_write_1orc(self):
        path = full_path('1orc.pdb')
        with open(path) as f:
            expected_lines = [line for line in f if is_written_to_pdb(line)]
        out_lines = self.write_back_and_compare(path)
        self.assertEqual(expected_lines, out_lines)

    def test_read_write_1lzh(self):
        path = full_path('1lzh.pdb.gz')
        mode = 'rt' if sys.version_info >= (3,) else 'r'
        with gzip.open(path, mode=mode) as f:
            expected_lines = [line for line in f if is_written_to_pdb(line)]
        out_lines = self.write_back_and_compare(path)
        self.assertEqual(expected_lines[0], out_lines[0])
        # TITLE lines differ because the text is broken at different word
        self.assertEqual(expected_lines[3:], out_lines[3:])

    def test_ncs_in_1lzh(self):
        st = gemmi.read_structure(full_path('1lzh.pdb.gz'))
        self.assertEqual(len(st.ncs), 1)
        A, B = st[0]
        for ra, rb in zip(A, B):
            pa = ra['CA'].pos
            pb = rb['CA'].pos
            image_of_pb = st.ncs[0].apply(pb)
            self.assertTrue(pa.dist(image_of_pb) < 0.01)

    def test_pdb_fragment(self):
        pdb_line = "HETATM 4154 MG    MG A 341       1.384  19.340  11.968" \
                   "  1.00 67.64          MG"
        for line in [pdb_line, pdb_line.strip(' MG')]:
            st = gemmi.read_pdb_string(pdb_line)
            mg_atom = st[0]['A']['341'][0]['MG']
            self.assertEqual(mg_atom.element.name, 'Mg')
            self.assertAlmostEqual(mg_atom.b_iso, 67.64, delta=1e-6)

    def test_ncs(self):
        st = gemmi.read_structure(full_path('5cvz_final.pdb'))
        chain = st[0]['A']
        first_atom = chain['17'][0]['N']
        ne2 = chain['63'][0]['NE2']
        direct_dist = first_atom.pos.dist(ne2.pos)
        self.assertAlmostEqual(direct_dist, 34.89, delta=1e-2)
        nearest_image = st.cell.find_nearest_image(first_atom.pos, ne2.pos)
        nearest_dist = nearest_image.dist()
        self.assertAlmostEqual(nearest_dist, 8.02, delta=1e-2)

if __name__ == '__main__':
    unittest.main()
