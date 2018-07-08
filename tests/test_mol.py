#!/usr/bin/env python

import gzip
import os
import sys
import tempfile
import unittest
import gemmi

def is_written_to_pdb(line, via_cif):
    if line[:6] in ['COMPND', 'SOURCE', 'MDLTYP', 'AUTHOR', 'REVDAT', 'JRNL  ',
                    'DBREF ', 'SEQADV', 'HET   ', 'HETNAM', 'FORMUL', 'HELIX ',
                    'SHEET ', 'SITE  ', 'MASTER']:
        return False
    # currently gemmi reads resolution from mmCIF but does not write it
    # to mmCIF. Otherwise it would be (via_cif and line[6:10] != '   2')
    if line[:6] == 'REMARK' and via_cif:
        return False
    return True

# $ zgrep -P '(^HEADER|CRYST1|SSBOND|ATOM.*SG)' pdb5cfg.ent.gz
SSBOND_FRAGMENT = """\
HEADER    LYASE                                   08-JUL-15   5CFG              
SSBOND   1 CYS A  138    CYS A  138                          1555   2555  2.48  
CRYST1   86.530   45.120   77.980  90.00 105.15  90.00 C 1 2 1       4          
ATOM    166  SG  CYS A  65      17.771 -16.223 -10.059  1.00 39.15           S  
ATOM    396  SG  CYS A  93      10.163 -20.624 -11.577  1.00 25.81           S  
ATOM    444  SG  CYS A  99      12.757 -34.484  -2.237  1.00 44.66           S  
ATOM    744  SG  CYS A 138       0.524 -16.872  -1.123  1.00 38.43           S  
ATOM   1309  SG  CYS A 208      10.668 -20.500 -14.891  1.00 25.60           S  
ATOM   2024  SG  CYS A 296      14.668  -4.407 -16.359  1.00 25.70           S  
ATOM   2128  SG  CYS A 310      24.141 -22.158 -17.213  1.00 20.07           S  
"""  # noqa: W291 - trailing whitespace

def full_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)

class TestMol(unittest.TestCase):
    def test_residue(self):
        res = gemmi.Residue()
        self.assertEqual(res.label_seq, None)
        self.assertEqual(res.seqid.num, None)
        res.label_seq = 1
        self.assertEqual(res.label_seq, 1)
        res.seqid.num = 2
        self.assertEqual(res.seqid.num, 2)
        res.label_seq = None
        self.assertEqual(res.label_seq, None)
        res.seqid.num = None
        self.assertEqual(res.seqid.num, None)
        for name in ['HOH', 'hoh', 'DOD', 'h2o', 'H2O', 'WAT']:
            res.name = name
            self.assertTrue(res.is_water())
        for name in ['SO4', '', 'HO', 'hoho', 'oho']:
            res.name = name
            self.assertFalse(res.is_water())

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
        self.assertEqual(st.info['_entry.id'], '5I55')
        chain, = st[0]
        a, b, c, d = chain.subchains()
        ent_a = st.get_entity_of(a)
        self.assertEqual(ent_a.name, '1')
        self.assertEqual(ent_a.entity_type, gemmi.EntityType.Polymer)
        self.assertEqual(ent_a.polymer_type, gemmi.PolymerType.PeptideL)
        ent_b = st.get_entity_of(b)
        self.assertEqual(ent_b.entity_type, gemmi.EntityType.NonPolymer)
        self.assertEqual(ent_b.polymer_type, gemmi.PolymerType.Unknown)
        ent_d = st.get_entity('4')
        self.assertEqual(ent_d.subchains, ['D'])
        self.assertEqual(ent_d.entity_type, gemmi.EntityType.Water)
        self.assertEqual(ent_d.polymer_type, gemmi.PolymerType.Unknown)

    def test_5i55_predefined_removals(self, clear_entities=False):
        st = gemmi.read_structure(full_path('5i55.cif'))
        if clear_entities:
            self.assertEqual(len(st.entities), 4)
            st.entities = gemmi.VectorEntity()
            self.assertEqual(len(st.entities), 0)
        lys12 = st[0]['A']['12']['LYS']
        count_b = sum(a.altloc == 'B' for a in lys12)
        model = st[0]
        # one author-chain and 4 label-chains: AA, 2 x ligand, waters
        self.assertEqual(len(model), 1)
        self.assertEqual(len(model.subchains()), 4)
        n_waters = len(model.get_subchain('D'))
        n_sites = model.count_atom_sites()
        occ_sum = model.count_occupancies()
        self.assertEqual(n_sites, occ_sum + count_b)
        st.remove_waters()
        self.assertEqual(len(model.subchains()), 3)
        self.assertEqual(model.count_atom_sites(), n_sites - n_waters)
        self.assertEqual(model.count_occupancies(), occ_sum - n_waters)
        st.remove_empty_chains()
        self.assertEqual(len(model.subchains()), 3)
        n_res = len(model['A'])
        st.remove_ligands_and_waters()
        self.assertEqual(len(model['A']), n_res - 2)
        model['A'].trim_to_alanine()
        # ALA has 5 atoms, except the last one which has OXT (hence +1)
        expected_count = sum(4 + (r.name != 'GLY') for r in model['A']) + 1
        self.assertEqual(model.count_occupancies(), expected_count)

    def test_5i55_predefined_removals2(self):
        self.test_5i55_predefined_removals(clear_entities=True)

    def test_rnase_predefined_removals(self, add_entities=False):
        st = gemmi.read_structure(full_path('rnase_frag.pdb'))
        if add_entities:
            self.assertEqual(len(st.entities), 0)
            st.assign_subchains()
            st.ensure_entities()
            self.assertEqual(len(st.entities), 4)
        model = st[0]
        nres_a = len(model['A'])
        nres_b = len(model['B'])
        st.remove_ligands_and_waters()  # removes SO4 from each chain
        self.assertEqual(len(model['A']), nres_a - 1)
        self.assertEqual(len(model['B']), nres_b - 1)
        self.assertEqual(len(model['W']), 0)
        self.assertEqual(len(model), 3)
        st.remove_empty_chains()
        self.assertEqual(len(model), 2)

    def test_rnase_predefined_removals2(self):
        self.test_rnase_predefined_removals(add_entities=True)

    def read_1pfe(self, filename):
        st = gemmi.read_structure(full_path(filename))
        self.assertAlmostEqual(st.cell.a, 39.374)
        self.assertEqual(st.cell.gamma, 120)
        self.assertEqual(st.name, '1PFE')
        self.assertEqual(st.spacegroup_hm, 'P 63 2 2')
        self.assertEqual(st.info['_entry.id'], '1PFE')
        self.assertEqual(st.info['_exptl.method'], 'X-RAY DIFFRACTION')
        self.assertEqual(len(st), 1)
        self.assertEqual(len(st[0]), 2)
        label_to_auth_name = {sub.name(): ch.name for ch in st[0]
                              for sub in ch.subchains()}
        self.assertEqual(label_to_auth_name,
                         dict(A='A', B='B', C='A', D='B', E='B', F='A', G='B'))
        self.assertEqual(len(st[0]['A']['1']), 1)
        chain_a = st[0]['A']
        self.assertEqual(chain_a[0].label_seq, 1)
        self.assertEqual(chain_a['1'][0].seqid.num, 1)
        b3 = st[0]['B']['3']
        self.assertEqual(len(b3), 2)
        self.assertEqual(repr(b3[0]), repr(st[0]['B'][2]))
        self.assertEqual(b3[0].name, 'N2C')
        self.assertEqual(b3[-1].name, 'NCY')
        chain_c = st[0].get_subchain('C')
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
        self.assertEqual(len(model), 1)
        self.assertEqual(len(model.subchains()), 2)
        A = model['A']
        waters = A.get_waters()
        self.assertEqual(len(waters), 57)  # FORMUL   2  HOH   *57(H2 O)
        self.assertTrue(all(res.name == 'HOH' for res in waters))
        self.assertTrue(A['3'])
        self.assertFalse(A['0'])
        self.assertEqual([res.seqid.num for res in A if res.seqid.icode != ' '],
                         [56] * 5)
        self.assertEqual(len(A['55']), 1)
        self.assertEqual(len(A['55B']), 0)
        self.assertEqual(len(A['56B']), 1)
        self.assertEqual(A['56'][0].seqid.icode, ' ')
        self.assertEqual(A['56c'][0].seqid.icode, 'C')

    def write_back_and_compare(self, path, via_cif):
        st = gemmi.read_structure(path)
        if via_cif:
            doc = st.make_mmcif_document()
            st = gemmi.make_structure_from_block(doc[0])
        handle, out_name = tempfile.mkstemp()
        os.close(handle)
        st.write_pdb(out_name)
        with open(out_name) as f:
            out_lines = f.readlines()
        os.remove(out_name)
        return out_lines

    def test_read_write_1orc(self, via_cif=False):
        path = full_path('1orc.pdb')
        with open(path) as f:
            expected = [line for line in f if is_written_to_pdb(line, via_cif)]
        out_lines = self.write_back_and_compare(path, via_cif)
        self.assertEqual(expected, out_lines)

    def test_read_write_1orc_via_cif(self):
        self.test_read_write_1orc(via_cif=True)

    def test_read_write_1lzh(self, via_cif=False):
        path = full_path('1lzh.pdb.gz')
        mode = 'rt' if sys.version_info >= (3,) else 'r'
        with gzip.open(path, mode=mode) as f:
            expected = [line for line in f if is_written_to_pdb(line, via_cif)]
        out_lines = self.write_back_and_compare(path, via_cif)
        self.assertEqual(expected[0], out_lines[0])
        # TITLE lines differ because the text is broken at different word
        self.assertEqual(expected[3:], out_lines[3:])

    def test_read_write_1lzh_via_cif(self):
        self.test_read_write_1lzh(via_cif=True)

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
        for line in [pdb_line, pdb_line.strip(' MG'), pdb_line[:-2] + '  ']:
            st = gemmi.read_pdb_string(line)
            mg_atom = st[0].sole_residue('A', 341, ' ')['MG']
            self.assertEqual(mg_atom.element.name, 'Mg')
            self.assertAlmostEqual(mg_atom.b_iso, 67.64, delta=1e-6)

    def test_4hhh_frag(self):
        path = full_path('4hhh_frag.pdb')
        with open(path) as f:
            frag = f.read()
        st = gemmi.read_pdb_string(frag)
        in_headers = frag.splitlines()
        out_headers = st.make_pdb_headers().splitlines()
        self.assertEqual(in_headers[0], out_headers[0])
        # the difference 4555 vs 2555 doesn't matter for us
        self.assertEqual(in_headers[1], out_headers[1].replace(' 4555 ',
                                                               ' 2555 '))
        self.assertEqual(in_headers[2], out_headers[2])

    def test_ncs(self):
        st = gemmi.read_structure(full_path('5cvz_final.pdb'))
        first_atom = st[0].sole_residue('A', 17, ' ')['N']
        ne2 = st[0].sole_residue('A', 63, ' ')['NE2']
        direct_dist = first_atom.pos.dist(ne2.pos)
        self.assertAlmostEqual(direct_dist, 34.89, delta=1e-2)
        nearest_image = st.cell.find_nearest_image(first_atom.pos, ne2.pos)
        nearest_dist = nearest_image.dist()
        self.assertAlmostEqual(nearest_dist, 8.02, delta=1e-2)

    def test_ssbond(self):
        st = gemmi.read_pdb_string(SSBOND_FRAGMENT)
        out = st.make_pdb_headers()
        self.assertEqual(out.splitlines(), SSBOND_FRAGMENT.splitlines()[:3])

    def test_ssbond_again(self):
        st = gemmi.read_pdb_string(SSBOND_FRAGMENT)
        doc = st.make_mmcif_document()
        st2 = gemmi.make_structure_from_block(doc[0])
        out = st2.make_pdb_headers()
        self.assertEqual(out.splitlines(), SSBOND_FRAGMENT.splitlines()[:3])

    def test_remove_atom(self):
        st = gemmi.read_pdb_string(SSBOND_FRAGMENT)
        res = st[0].sole_residue('A', 310, ' ')
        self.assertEqual(len(res), 1)
        del res['SG']
        self.assertEqual(len(res), 0)
        self.assertEqual(len(st), 1)
        self.assertEqual(st[0].name, '1')
        del st['1']
        self.assertEqual(len(st), 0)

    def test_extract_sequence_info(self):
        st = gemmi.read_structure(full_path('5cvz_final.pdb'))
        st.add_entity_types()
        polymer = st[0][0].get_polymer()
        self.assertEqual(polymer.check_polymer_type(),
                         gemmi.PolymerType.PeptideL)
        expected = ('AAATSLVYDTCYVTLTERATTSFQRQSFPTLKGMGDRAFQVVAFTIQGVS'
                    'AAPLMYNARLYNPGDTDSVHATGVQLMGTVPRTVRLTPRVGQNNWFFGNT'
                    'EEAETILAIDGLVSTKGANAPSNTVIVTGCFRLAPSELQSS')
        self.assertEqual(polymer.make_one_letter_sequence(), expected)


if __name__ == '__main__':
    unittest.main()
