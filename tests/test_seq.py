#!/usr/bin/env python

import io
import unittest
import gemmi
try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None

# https://rest.uniprot.org/uniprotkb/P00698.fasta
FASTA1 = """\
>sp|P00698|LYSC_CHICK Lysozyme C OS=Gallus gallus OX=9031 GN=LYZ PE=1 SV=1
MRSLLILVLCFLPLAALGKVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQA
TNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDG
NGMNAWVAWRNRCKGTDVQAWIRGCRL
"""  # noqa: W291 - trailing whitespace

# https://en.wikipedia.org/wiki/FASTA_format
FASTA2 = """\
>MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken
MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTID
FPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREA
DIDGDGQVNYEEFVQMMTAK*

>gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV
EWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG
LLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVIL
GLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX
IENY
"""  # noqa: W291 - trailing whitespace

# https://www.metagenomics.wiki/tools/fastq/multi-fasta-format
FASTA3 = """\
>sequenceID-001 description
AAGTAGGAATAATATCTTATCATTATAGATAAAAACCTTCTGAATTTGCTTAGTGTGTAT
ACGACTAGACATATATCAGCTCGCCGATTATTTGGATTATTCCCTG
>sequenceID-002 description
CAGTAAAGAGTGGATGTAAGAACCGTCCGATCTACCAGATGTGATAGAGGTTGCCAGTAC
AAAAATTGCATAATAATTGATTAATCCTTTAATATTGTTTAGAATATATCCGTCAGATAA
TCCTAAAAATAACGATATGATGGCGGAAATCGTC
>sequenceID-003 description
CTTCAATTACCCTGCTGACGCGAGATACCTTATGCATCGAAGGTAAAGCGATGAATTTAT
CCAAGGTTTTAATTTG
"""  # noqa: W291 - trailing whitespace

# from _entity_poly.pdbx_seq_one_letter_code from 5I55 and 1PFE
FASTA4 = """\
>
(MSE)EFVAKLFKFFKDLLGKFLGNN
>
(DSN)A(N2C)(MVA)(DSN)A(NCY)(MVA)
"""

# https://biopython.org/docs/1.81/api/Bio.SeqIO.PirIO.html
PIR1 = """\
>P1;S27231
rhodopsin - northern leopard frog
MNGTEGPNFY IPMSNKTGVV RSPFDYPQYY LAEPWKYSVL AAYMFLLILL GLPINFMTLY
VTIQHKKLRT PLNYILLNLG VCNHFMVLCG FTITMYTSLH GYFVFGQTGC YFEGFFATLG
GEIALWSLVV LAIERYIVVC KPMSNFRFGE NHAMMGVAFT WIMALACAVP PLFGWSRYIP
EGMQCSCGVD YYTLKPEVNN ESFVIYMFVV HFLIPLIIIS FCYGRLVCTV KEAAAQQQES
ATTQKAEKEV TRMVIIMVIF FLICWVPYAY VAFYIFTHQG SEFGPIFMTV PAFFAKSSAI
YNPVIYIMLN KQFRNCMITT LCCGKNPFGD DDASSAATSK TEATSVSTSQ VSPA*

>P1;I51200
rhodopsin - African clawed frog
MNGTEGPNFY VPMSNKTGVV RSPFDYPQYY LAEPWQYSAL AAYMFLLILL GLPINFMTLF
VTIQHKKLRT PLNYILLNLV FANHFMVLCG FTVTMYTSMH GYFIFGPTGC YIEGFFATLG
GEVALWSLVV LAVERYIVVC KPMANFRFGE NHAIMGVAFT WIMALSCAAP PLFGWSRYIP
EGMQCSCGVD YYTLKPEVNN ESFVIYMFIV HFTIPLIVIF FCYGRLLCTV KEAAAQQQES
LTTQKAEKEV TRMVVIMVVF FLICWVPYAY VAFYIFTHQG SNFGPVFMTV PAFFAKSSAI
YNPVIYIVLN KQFRNCLITT LCCGKNPFGD EDGSSAATSK TEASSVSSSQ VSPA*

>P1;JN0120
rhodopsin - Japanese lamprey
MNGTEGDNFY VPFSNKTGLA RSPYEYPQYY LAEPWKYSAL AAYMFFLILV GFPVNFLTLF
VTVQHKKLRT PLNYILLNLA MANLFMVLFG FTVTMYTSMN GYFVFGPTMC SIEGFFATLG
GEVALWSLVV LAIERYIVIC KPMGNFRFGN THAIMGVAFT WIMALACAAP PLVGWSRYIP
EGMQCSCGPD YYTLNPNFNN ESYVVYMFVV HFLVPFVIIF FCYGRLLCTV KEAAAAQQES
ASTQKAEKEV TRMVVLMVIG FLVCWVPYAS VAFYIFTHQG SDFGATFMTL PAFFAKSSAL
YNPVIYILMN KQFRNCMITT LCCGKNPLGD DE-SGASTSKT EVSSVSTSPV SPA*
"""  # noqa: W291 - trailing whitespace

# https://www.bioinformatics.nl/tools/crab_pir.html
PIR2 = """\
>P1;CRAB_ANAPL
ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN).
  MDITIHNPLI RRPLFSWLAP SRIFDQIFGE HLQESELLPA SPSLSPFLMR 
  SPIFRMPSWL ETGLSEMRLE KDKFSVNLDV KHFSPEELKV KVLGDMVEIH 
  GKHEERQDEH GFIAREFNRK YRIPADVDPL TITSSLSLDG VLTVSAPRKQ 
  SDVPERSIPI TREEKPAIAG AQRK*

>P1;CRAB_BOVIN
ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN).
  MDIAIHHPWI RRPFFPFHSP SRLFDQFFGE HLLESDLFPA STSLSPFYLR 
  PPSFLRAPSW IDTGLSEMRL EKDRFSVNLD VKHFSPEELK VKVLGDVIEV 
  HGKHEERQDE HGFISREFHR KYRIPADVDP LAITSSLSSD GVLTVNGPRK 
  QASGPERTIP ITREEKPAVT AAPKK*
"""  # noqa: W291 - trailing whitespace

@unittest.skipIf(SeqIO is None, "BioPython not installed.")
class TestReadingSeq(unittest.TestCase):
    def test_fasta(self):
        for string in (FASTA1, FASTA2, FASTA3, FASTA4):
            gemmi_seqs = gemmi.read_pir_or_fasta(string)
            biopy_seqs = list(SeqIO.parse(io.StringIO(string), 'fasta'))
            self.assertEqual(len(gemmi_seqs), len(biopy_seqs))
            for (gseq, bseq) in zip(gemmi_seqs, biopy_seqs):
                self.assertEqual(gseq.header, bseq.description)
                self.assertEqual(gseq.seq, bseq.seq.rstrip('*'))

    def test_pir(self):
        for string in (PIR1, PIR2):
            gemmi_seqs = gemmi.read_pir_or_fasta(string)
            biopy_seqs = list(SeqIO.parse(io.StringIO(string), 'pir'))
            self.assertEqual(len(gemmi_seqs), len(biopy_seqs))
            for (gseq, bseq) in zip(gemmi_seqs, biopy_seqs):
                g_first, g_desc = gseq.header.splitlines()
                g_type, _, g_id = g_first.partition(';')
                self.assertEqual(g_type, bseq.annotations["PIR-type"])
                self.assertEqual(g_id, bseq.id)
                self.assertEqual(g_desc, bseq.description)
                self.assertEqual(gseq.seq, bseq.seq)

    def test_code_conversion_aa(self):
        seq1 = gemmi.read_pir_or_fasta(FASTA1)[0].seq
        seq3 = gemmi.expand_one_letter_sequence(seq1, gemmi.ResidueKind.AA)
        self.assertEqual(seq1, gemmi.one_letter_code(seq3))

    def test_code_conversion_dna(self):
        seq1 = gemmi.read_pir_or_fasta(FASTA3)[0].seq
        kind = gemmi.ResidueKind.DNA
        seq3 = gemmi.expand_one_letter_sequence(seq1, kind)
        self.assertEqual(seq1, gemmi.one_letter_code(seq3))
        self.assertEqual(seq1, gemmi.pdbx_one_letter_code(seq3, kind))

    def test_code_conversion_rna(self):
        seq1 = 'GGCGAUACCAGCCGAAAGGCCCUUGGCAGCGCC'  # from 8d2b
        kind = gemmi.ResidueKind.RNA
        seq3 = gemmi.expand_one_letter_sequence(seq1, kind)
        self.assertEqual(seq1, gemmi.one_letter_code(seq3))
        self.assertEqual(seq1, gemmi.pdbx_one_letter_code(seq3, kind))

    def test_code_with_brackets(self):
        # test 1PFE _entity_poly.pdbx_seq_one_letter_code[_can]
        seq1 = gemmi.read_pir_or_fasta(FASTA4)[1].seq
        kind = gemmi.ResidueKind.AA
        seq3 = gemmi.expand_one_letter_sequence(seq1, kind)
        self.assertEqual(seq3, ['DSN', 'ALA', 'N2C', 'MVA',
                                'DSN', 'ALA', 'NCY', 'MVA'])
        self.assertEqual(seq1, gemmi.pdbx_one_letter_code(seq3, kind))

if __name__ == '__main__':
    unittest.main()
