#!/usr/bin/env python

import unittest
import gemmi
from common import full_path


def _make_peptide_chain(model, chain_name, subchain_id,
                            residue_names, spacing=3.8):
    """Helper: build a polypeptide chain with backbone atoms (N, CA, C)."""
    chain = model.add_chain(chain_name)
    for i, name in enumerate(residue_names):
        res = gemmi.Residue()
        res.name = name
        res.seqid = gemmi.SeqId(str(i + 1))
        res.entity_type = gemmi.EntityType.Polymer
        res.subchain = subchain_id
        x = (i + 1) * spacing
        for aname, el_str, dx in [('N', 'N', -1.3),
                                    ('CA', 'C', 0.0),
                                    ('C', 'C', 1.3)]:
            atom = gemmi.Atom()
            atom.name = aname
            atom.element = gemmi.Element(el_str)
            atom.pos = gemmi.Position(x + dx, 0, 0)
            res.add_atom(atom)
        chain.add_residue(res)
    return chain

def _make_cation_residue(chain, name, seqid, subchain_id):
    """Helper: add a cation residue to a chain as part of the polymer.

    In many real PDB files (e.g. from refinement software), cation
    residues like Zn and Ca end up inside the polymer subchain.  The
    alignment code must filter them out."""
    res = gemmi.Residue()
    res.name = name
    res.seqid = gemmi.SeqId(str(seqid))
    res.entity_type = gemmi.EntityType.Polymer
    res.subchain = subchain_id
    atom = gemmi.Atom()
    atom.name = name
    atom.element = gemmi.Element(name[:2].strip())
    atom.pos = gemmi.Position(100, 0, 0)
    res.add_atom(atom)
    chain.add_residue(res)

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


    def test_assign_sequences_ignores_cations(self):
        """Cations (Ca, Zn, etc.) that are part of the polymer subchain
        (as happens in real PDB files from refinement) should be filtered
        out during alignment and not prevent FASTA assignment."""
        st = gemmi.Structure()
        model = gemmi.Model('1')
        protein_residues = ['ALA'] * 10  # 10-residue polyalanine
        chain = _make_peptide_chain(model, 'A', 'Axp', protein_residues)
        # Cations inside the polymer subchain
        _make_cation_residue(chain, 'Zn', 11, 'Axp')
        _make_cation_residue(chain, 'Ca', 12, 'Axp')  # calcium
        st.add_model(model)

        # Set up entity manually to have known polymer_type
        ent = gemmi.Entity('1')
        ent.entity_type = gemmi.EntityType.Polymer
        ent.polymer_type = gemmi.PolymerType.PeptideL
        ent.subchains = ['Axp']
        st.entities.append(ent)

        fasta = 'AAAAAAAAAA'  # exactly 10 Ala
        st.assign_best_sequences([fasta])

        assigned = st.entities[0].full_sequence
        self.assertEqual(len(assigned), 10,
                         f'Expected 10 residues but got {len(assigned)}: {assigned}')
        self.assertTrue(all(r == 'ALA' for r in assigned))

    def test_assign_sequences_longer_fasta_no_spillover(self):
        """When FASTA is longer than the chain (unmodeled N/C-termini),
        the entire FASTA should be assigned, even when cations are inside the polymer 
        subchain."""
        st = gemmi.Structure()
        model = gemmi.Model('1')
        # Structure has residues 3-8 of a 10-residue protein
        protein_residues = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET']
        chain = _make_peptide_chain(model, 'A', 'Axp', protein_residues)
        # Cation inside the polymer subchain (realistic PDB scenario)
        _make_cation_residue(chain, 'Zn', 7, 'Axp')
        st.add_model(model)

        ent = gemmi.Entity('1')
        ent.entity_type = gemmi.EntityType.Polymer
        ent.polymer_type = gemmi.PolymerType.PeptideL
        ent.subchains = ['Axp']
        st.entities.append(ent)

        # FASTA covers the full biological sequence including unmodeled termini
        fasta = 'AAGAVLIMAA'  # 10-mer: AA + GAVLIM + AA
        st.assign_best_sequences([fasta])

        assigned = st.entities[0].full_sequence
        # The full FASTA should be assigned as the SEQRES
        self.assertEqual(len(assigned), 10)


    def test_assign_sequences_resolves_unknown_polymer_type(self):
        """An entity with PolymerType.Unknown should have its type resolved
        from chain content before alignment, so that it can be matched to
        the correct FASTA sequence.  On main (before the fix), Unknown
        entities are never processed because none of the three polymer-type
        passes (PeptideL, Rna, Dna) match Unknown."""
        st = gemmi.Structure()
        model = gemmi.Model('1')

        # Protein chain with standard amino acids
        protein_residues = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE',
                            'MET', 'PHE', 'TRP']
        _make_peptide_chain(model, 'A', 'Axp', protein_residues)
        st.add_model(model)

        # Entity deliberately set to Unknown (as happens when
        # setup_entities() wasn't called or detection was incomplete)
        ent = gemmi.Entity('1')
        ent.entity_type = gemmi.EntityType.Polymer
        ent.polymer_type = gemmi.PolymerType.Unknown
        ent.subchains = ['Axp']
        st.entities.append(ent)

        fasta = 'GAVLIMFW'  # matches the 8 residues exactly
        st.assign_best_sequences([fasta])

        assigned = st.entities[0].full_sequence
        self.assertTrue(len(assigned) > 0,
                        'Entity with Unknown polymer_type got no sequence '
                        'assigned — type resolution is not working')
        self.assertEqual(len(assigned), 8)
        expected = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP']
        self.assertEqual(assigned, expected)

    def test_assign_sequences_connectivity_aware_gap_penalties(self):
        """Gap penalties should be based on backbone connectivity.

        When there is a break in the backbone (large CA-CA distance,
        simulating missing/unmodeled residues), gap opening at that
        position should be free (good_gapo=0) so that extra FASTA
        residues align into the break.

        We build a 6-residue chain in two segments separated by a
        large gap (50 A between segment CAs), then provide a 10-residue
        FASTA with 4 extra residues that belong in the break."""
        st = gemmi.Structure()
        model = gemmi.Model('1')
        chain = model.add_chain('A')

        # Segment 1: residues 1-3, closely spaced (3.8 A apart -> connected)
        seg1 = ['GLY', 'ALA', 'VAL']
        # Segment 2: residues 4-6, closely spaced but far from segment 1
        seg2 = ['LEU', 'ILE', 'MET']
        seg1_x_start = 0.0
        seg2_x_start = 50.0  # 50 A away -> CA-CA >> 5 A -> disconnected
        spacing = 3.8
        seqid_counter = 0
        for seg_names, x_start in [(seg1, seg1_x_start),
                                    (seg2, seg2_x_start)]:
            for j, name in enumerate(seg_names):
                seqid_counter += 1
                res = gemmi.Residue()
                res.name = name
                res.seqid = gemmi.SeqId(str(seqid_counter))
                res.entity_type = gemmi.EntityType.Polymer
                res.subchain = 'Axp'
                x = x_start + j * spacing
                for aname, el_str, dx in [('N', 'N', -1.3),
                                           ('CA', 'C', 0.0),
                                           ('C', 'C', 1.3)]:
                    atom = gemmi.Atom()
                    atom.name = aname
                    atom.element = gemmi.Element(el_str)
                    atom.pos = gemmi.Position(x + dx, 0, 0)
                    res.add_atom(atom)
                chain.add_residue(res)
        st.add_model(model)
        ent = gemmi.Entity('1')
        ent.entity_type = gemmi.EntityType.Polymer
        ent.polymer_type = gemmi.PolymerType.PeptideL
        ent.subchains = ['Axp']
        st.entities.append(ent)
        # FASTA: GAV + PFYW (4 unmodeled residues in the break) + LIM = 10
        fasta = 'GAVPFYWLIM'
        st.assign_best_sequences([fasta])

        assigned = st.entities[0].full_sequence
        self.assertEqual(len(assigned), 10,
                         f'Expected 10-residue FASTA assigned across the '
                         f'backbone break, got {len(assigned)}: {assigned}')
        expected = ['GLY', 'ALA', 'VAL', 'PRO', 'PHE', 'TYR', 'TRP',
                    'LEU', 'ILE', 'MET']
        self.assertEqual(assigned, expected)

    def test_assign_sequences_many_trailing_cations(self):
        """Ensure a long protein chain with several cations (Zn, Ca) appended to the 
        end of a polymer subchain still has proper sequence assignment."""
        st = gemmi.Structure()
        model = gemmi.Model('1')
        protein_residues = ['ALA'] * 50
        chain = _make_peptide_chain(model, 'A', 'Axp', protein_residues)
        # Append multiple cations as part of the polymer subchain
        for j, ion in enumerate(['Zn', 'Zn', 'Ca', 'Ca', 'Ca']):
            _make_cation_residue(chain, ion, 51 + j, 'Axp')
        st.add_model(model)

        ent = gemmi.Entity('1')
        ent.entity_type = gemmi.EntityType.Polymer
        ent.polymer_type = gemmi.PolymerType.PeptideL
        ent.subchains = ['Axp']
        st.entities.append(ent)

        fasta = 'A' * 50
        st.assign_best_sequences([fasta])

        assigned = st.entities[0].full_sequence
        self.assertEqual(len(assigned), 50,
                         f'Expected 50 residues but got {len(assigned)}')

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
