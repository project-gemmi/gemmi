#!/usr/bin/env python

import os
import tempfile
import unittest
import gemmi

def is_written_to_pdb(line):
    if line[:6] in [b'COMPND', b'SOURCE', b'AUTHOR', b'REVDAT', b'JRNL  ',
                    b'DBREF ', b'SEQADV', b'FORMUL', b'HELIX ', b'SHEET ',
                    b'MASTER']:
        return False
    if line[:6] == b'REMARK' and line[6:10] != b'   2':
        return False
    return True

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
        path = os.path.join(os.path.dirname(__file__), '5i55.cif')
        cell = gemmi.read_structure(path).cell
        self.assertAlmostEqual(cell.a, 29.46)
        self.assertAlmostEqual(cell.b, 10.51)
        self.assertAlmostEqual(cell.c, 29.71)
        self.assertEqual(cell.alpha, 90)
        self.assertAlmostEqual(cell.beta, 111.98)
        self.assertEqual(cell.gamma, 90)

    def test_read_5i55_again(self):
        path = os.path.join(os.path.dirname(__file__), '5i55.cif')
        st = gemmi.read_structure(path)
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

    def test_read_1pfe(self):
        path = os.path.join(os.path.dirname(__file__), '1pfe.cif.gz')
        st = gemmi.read_structure(path)
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

    def test_read_1orc(self):
        path = os.path.join(os.path.dirname(__file__), '1orc.pdb')
        st = gemmi.read_structure(path)
        self.assertAlmostEqual(st.cell.a, 34.77)
        self.assertEqual(st.cell.alpha, 90)
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

    def test_read_write(self):
        path = os.path.join(os.path.dirname(__file__), '1orc.pdb')
        st = gemmi.read_structure(path)
        with open(path, 'rb') as f:
            expected_lines = [line for line in f if is_written_to_pdb(line)]
        with tempfile.NamedTemporaryFile() as f:
            out = st.write_pdb(f.name)
            out_lines = f.readlines()
        self.assertEqual(expected_lines, out_lines)

if __name__ == '__main__':
    unittest.main()
