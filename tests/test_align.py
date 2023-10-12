#!/usr/bin/env python

import unittest
import gemmi
from common import full_path

class TestAlignment(unittest.TestCase):
    def test_string_align(self):
        result = gemmi.align_string_sequences(list('AABCC'),
                                              list('ABC'), [0])
        self.assertEqual(result.score, 0)
        self.assertEqual(result.cigar_str(), '1I3M1I')
        self.assertEqual(result.add_gaps('AABCC', 1), 'AABCC')
        self.assertEqual(result.add_gaps('ABC', 2), '-ABC-')
        self.assertEqual(result.calculate_identity(), 100.)
        self.assertEqual(result.calculate_identity(1), 60.)
        self.assertEqual(result.calculate_identity(2), 100.)
        self.assertEqual(result.match_string, ' ||| ')
        result = gemmi.align_string_sequences(list('SIMILARITY'),
                                              list('PILLAR'), [])
        self.assertEqual(result.match_count, 4)
        self.assertEqual(result.cigar_str(), '3M1I3M3I')

    def test_hemoglobin_alignment(self):
        # based on example from
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html
        hba_human = ("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQ"
                     "VKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTL"
                     "AAHLPAEFTPAVHASLDKFLASVSTVLTSKYR")
        hbb_human = ("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAV"
                     "MGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNV"
                     "LVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH")
        AA = gemmi.ResidueKind.AA
        hba_seq = gemmi.expand_one_letter_sequence(hba_human, AA)
        hbb_seq = gemmi.expand_one_letter_sequence(hbb_human, AA)
        id_score = gemmi.AlignmentScoring()
        id_score.match = 1
        id_score.mismatch = 0
        id_score.gapo = 0
        id_score.gape = 0
        result = gemmi.align_string_sequences(hba_seq, hbb_seq, [], id_score)
        # "80 different alignments with the score 72"
        self.assertEqual(result.score, 72)
        blosum62 = gemmi.AlignmentScoring('b')
        blosum62.gapo = -9
        result = gemmi.align_string_sequences(hba_seq, hbb_seq, [], blosum62)
        # BioPython equivalent is:
        # pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -10, -1)
        self.assertEqual(result.score, 290)

    def test_assign_best_sequences(self):
        st = gemmi.read_structure(full_path('1lzh.pdb.gz'))
        st.setup_entities()
        seq3 = st.entities[0].full_sequence
        self.assertEqual(len(seq3), 129)
        seq1 = gemmi.one_letter_code(seq3)
        self.assertEqual(len(seq1), 129)
        st.clear_sequences()
        self.assertEqual(len(st.entities[0].full_sequence), 0)
        st.assign_best_sequences([seq1+'A', seq1])
        self.assertEqual(len(st.entities[0].full_sequence), 129)
        self.assertEqual(st.entities[0].subchains, ['Axp', 'Bxp'])

    def test_superposition(self):
        model = gemmi.read_structure(full_path('4oz7.pdb'))[0]
        poly1 = model['A'].get_polymer()
        poly2 = model['B'].get_polymer()
        ptype = poly1.check_polymer_type()
        S = gemmi.SupSelect
        s1 = gemmi.calculate_superposition(poly1, poly2, ptype, S.CaP)
        s2 = gemmi.calculate_superposition(poly1, poly2, ptype, S.MainChain)
        s3 = gemmi.calculate_superposition(poly1, poly2, ptype, S.All)
        self.assertEqual(s1.count, 10)
        self.assertEqual(s2.count, 39)
        self.assertEqual(s3.count, 77)
        self.assertAlmostEqual(s1.rmsd, 0.146, places=3)
        self.assertAlmostEqual(s2.rmsd, 0.174, places=3)
        self.assertAlmostEqual(s3.rmsd, 0.400, places=3)
        for s in [s1, s2, s3]:
            self.assertAlmostEqual(s.transform.vec.y, 17.0, places=1)

if __name__ == '__main__':
    unittest.main()
