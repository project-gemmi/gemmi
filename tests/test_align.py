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
    residues like ZN and CA end up inside the polymer subchain.  The
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
        """Cations (CA, ZN, etc.) that are part of the polymer subchain
        (as happens in real PDB files from refinement) should be filtered
        out during alignment and not prevent FASTA assignment."""
        st = gemmi.Structure()
        model = gemmi.Model('1')
        protein_residues = ['ALA'] * 10  # 10-residue polyalanine
        chain = _make_peptide_chain(model, 'A', 'Axp', protein_residues)
        # Cations inside the polymer subchain (the real-world problematic case)
        _make_cation_residue(chain, 'ZN', 11, 'Axp')
        _make_cation_residue(chain, 'CA', 12, 'Axp')  # calcium
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
        _make_cation_residue(chain, 'ZN', 7, 'Axp')
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

    def test_assign_sequences_multi_chain_truncations(self):
        """Entity with multiple subchains (asymmetric unit copies) having
        different truncation levels should still pick the best FASTA."""
        st = gemmi.Structure()
        model = gemmi.Model('1')
        # Chain A: shorter truncation (residues 3-8 of a 10-residue protein)
        short_residues = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET']
        _make_peptide_chain(model, 'A', 'Axp', short_residues)
        # Chain B: longer truncation (residues 2-9)
        long_residues = ['ALA', 'GLY', 'ALA', 'VAL', 'LEU', 'ILE',
                         'MET', 'ALA']
        _make_peptide_chain(model, 'B', 'Bxp', long_residues)
        st.add_model(model)

        # Both subchains belong to the same entity
        ent = gemmi.Entity('1')
        ent.entity_type = gemmi.EntityType.Polymer
        ent.polymer_type = gemmi.PolymerType.PeptideL
        ent.subchains = ['Axp', 'Bxp']
        st.entities.append(ent)

        correct_fasta = 'AAGAVLIMAA'
        wrong_fasta = 'MMMMMMMMMM'
        st.assign_best_sequences([wrong_fasta, correct_fasta])

        assigned = st.entities[0].full_sequence
        self.assertEqual(len(assigned), 10)
        # Verify the correct FASTA was picked, not the wrong one
        assigned_1letter = gemmi.one_letter_code(assigned)
        self.assertEqual(assigned_1letter, correct_fasta)

    def test_assign_sequences_protein_not_assigned_to_dna(self):
        """A protein FASTA should not be assigned to a DNA entity."""
        st = gemmi.Structure()
        model = gemmi.Model('1')

        # Protein chain
        protein_residues = ['ALA', 'GLY', 'VAL', 'LEU', 'ILE']
        _make_peptide_chain(model, 'A', 'Axp', protein_residues)

        # DNA chain
        chain_b = model.add_chain('B')
        for i, name in enumerate(['DA', 'DT', 'DG', 'DC', 'DA'], start=1):
            res = gemmi.Residue()
            res.name = name
            res.seqid = gemmi.SeqId(str(i))
            res.entity_type = gemmi.EntityType.Polymer
            res.subchain = 'Bxp'
            atom = gemmi.Atom()
            atom.name = 'P'
            atom.element = gemmi.Element('P')
            atom.pos = gemmi.Position(i * 6.0, 0, 0)
            res.add_atom(atom)
            chain_b.add_residue(res)

        st.add_model(model)

        # Protein entity
        ent_prot = gemmi.Entity('1')
        ent_prot.entity_type = gemmi.EntityType.Polymer
        ent_prot.polymer_type = gemmi.PolymerType.PeptideL
        ent_prot.subchains = ['Axp']
        st.entities.append(ent_prot)

        # DNA entity
        ent_dna = gemmi.Entity('2')
        ent_dna.entity_type = gemmi.EntityType.Polymer
        ent_dna.polymer_type = gemmi.PolymerType.Dna
        ent_dna.subchains = ['Bxp']
        st.entities.append(ent_dna)

        # Only provide a protein FASTA (no DNA FASTA)
        protein_fasta = 'AGVLI'
        st.assign_best_sequences([protein_fasta])

        # Protein entity should get the sequence
        self.assertEqual(len(st.entities[0].full_sequence), 5)
        # DNA entity should NOT get assigned the protein FASTA
        self.assertEqual(len(st.entities[1].full_sequence), 0,
                         'Protein FASTA was incorrectly assigned to DNA entity')

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

    def test_assign_sequences_many_trailing_cations(self):
        """Realistic scenario: a long protein chain with several cations
        (ZN, CA) appended inside the polymer subchain, as seen in PDB files
        from refinement software.  On master, align_sequence_to_polymer
        includes these cations, causing the alignment score to be too low
        for assignment."""
        st = gemmi.Structure()
        model = gemmi.Model('1')
        protein_residues = ['ALA'] * 50
        chain = _make_peptide_chain(model, 'A', 'Axp', protein_residues)
        # Append multiple cations as part of the polymer subchain
        for j, ion in enumerate(['ZN', 'ZN', 'CA', 'CA', 'CA']):
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

    def test_assign_sequences_unknown_polymer_type(self):
        """When an entity has polymer_type=Unknown (e.g. setup_entities was
        not called), assign_best_sequences should auto-detect the polymer
        type from the chain content and still assign the FASTA."""
        st = gemmi.Structure()
        model = gemmi.Model('1')
        protein_residues = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET']
        _make_peptide_chain(model, 'A', 'Axp', protein_residues)
        st.add_model(model)

        ent = gemmi.Entity('1')
        ent.entity_type = gemmi.EntityType.Polymer
        ent.polymer_type = gemmi.PolymerType.Unknown  # not set
        ent.subchains = ['Axp']
        st.entities.append(ent)

        fasta = 'AAGAVLIMAA'  # longer FASTA with extra N/C-terminal residues
        st.assign_best_sequences([fasta])

        self.assertEqual(st.entities[0].polymer_type, gemmi.PolymerType.PeptideL,
                         'polymer_type should be auto-detected as PeptideL')
        assigned = st.entities[0].full_sequence
        self.assertEqual(len(assigned), 10,
                         f'Expected 10 residues but got {len(assigned)}: {assigned}')
        assigned_1letter = gemmi.one_letter_code(assigned)
        self.assertEqual(assigned_1letter, fasta)

    def test_assign_sequences_fasta_with_insertion(self):
        """FASTA with extra residues (prefix + mid-chain insertion) compared
        to the model should be fully assigned as SEQRES."""
        st = gemmi.Structure()
        model = gemmi.Model('1')
        # Model has GAVLIM (6 residues)
        protein_residues = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET']
        _make_peptide_chain(model, 'A', 'Axp', protein_residues)
        st.add_model(model)

        ent = gemmi.Entity('1')
        ent.entity_type = gemmi.EntityType.Polymer
        ent.polymer_type = gemmi.PolymerType.PeptideL
        ent.subchains = ['Axp']
        st.entities.append(ent)

        # FASTA: AA prefix + GAV + extra S + LIM + AA suffix = 13 residues
        fasta = 'AAGAVSLIMAA'
        st.assign_best_sequences([fasta])

        assigned = st.entities[0].full_sequence
        assigned_1letter = gemmi.one_letter_code(assigned)
        self.assertEqual(len(assigned), len(fasta),
                         f'Expected {len(fasta)} residues but got {len(assigned)}')
        self.assertEqual(assigned_1letter, fasta,
                         'SEQRES should match the full FASTA including extra residues')


if __name__ == '__main__':
    unittest.main()
