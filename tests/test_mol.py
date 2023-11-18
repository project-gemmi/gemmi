#!/usr/bin/env python

import gzip
from io import StringIO
import os
import sys
import unittest
import gemmi
from common import full_path, get_path_for_tempfile
try:
    from Bio import PDB
except ImportError:
    PDB = None

def is_written_to_pdb(line, via_cif):
    if line[:6] in ['COMPND', 'SOURCE', 'MDLTYP', 'AUTHOR', 'REVDAT', 'JRNL  ',
                    'SEQADV', 'HET   ', 'HETNAM', 'FORMUL',
                    'SITE  ', 'MASTER', 'CONECT']:
        return False
    # ORIGX is written only if it is a non-identity matrix
    # SCALE is written only if it is non-default
    if line[:5] in ['ORIGX', 'SCALE']:
        return False
    if line[:6] == 'REMARK' and via_cif and line[7:10] not in ['  2', '350']:
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

SHORT_SSBOND = """\
SSBOND   1 CYS A    6    CYS A   11
ATOM     37  SG  CYS A   6      38.416  25.985  18.085  1.00 14.07  16       S
ATOM     69  SG  CYS A  11      36.989  25.994  19.570  1.00 14.23  16       S
"""

# fragment of $CCP4/examples/toxd/toxd_mod_p1.pdb without chain names
BLANK_CHAIN_FRAGMENT = """\
CRYST1   50.000   50.00    50.00   90.00  90.00   90.00
ATOM      1  N   GLN     1      16.207  -7.425   8.244  1.00 35.23
ATOM      2  CA  GLN     1      16.175  -6.019   8.805  1.00 35.13
ATOM      3  C   GLN     1      14.890  -5.752   9.585  1.00 33.49
ATOM      4  O   GLN     1      14.563  -6.517  10.541  1.00 35.33
ATOM      5  N   PRO     2      14.235  -4.626   9.270  1.00 30.05
ATOM      6  CA  PRO     2      12.967  -4.255   9.939  1.00 26.63
ATOM      7  C   PRO     2      11.955  -5.328   9.522  1.00 24.15
ATOM      8  O   PRO     2      12.217  -6.028   8.540  1.00 22.83
ATOM      9  CB  PRO     2      12.617  -2.895   9.427  1.00 26.18
ATOM     10  CG  PRO     2      13.342  -2.761   8.127  1.00 28.07
ATOM     11  CD  PRO     2      14.574  -3.663   8.212  1.00 29.12
ATOM     12  N   ARG     3      10.892  -5.446  10.311  1.00 21.86
ATOM     13  CA  ARG     3       9.875  -6.426  10.119  1.00 19.82
"""

# from https://github.com/project-gemmi/gemmi/issues/37
AMBER_FRAGMENT = """\
ATOM      7  CB  VAL     1     -14.375 -11.856  27.866  1.00  0.00
ATOM      8  HB  VAL     1     -14.217 -11.118  27.080  1.00  0.00
ATOM      9  CG1 VAL     1     -13.033 -12.232  28.471  1.00  0.00
ATOM     10 HG11 VAL     1     -12.398 -12.673  27.702  1.00  0.00
"""

FRAGMENT_WITH_HG = """\
HETATM 4406 HG    HG P 693      28.820  31.751  40.919  0.20 25.99
HETATM 4407 HG1   HG P 694      27.455  32.086  39.686  0.20 35.18
"""

TRJCONV_FRAGMENT = """\
ATOM  12609 5C'N NPH  6378     421.300 400.390 491.570  1.00  0.00
ATOM  12610 1HN5 NPH  6378     422.020 400.650 490.790  1.00  0.00
ATOM  12611 2HN5 NPH  6378     420.490 399.940 490.980  1.00  0.00
ATOM  12612 5O'N NPH  6378     421.940 399.470 492.410  1.00  0.00
ATOM  12613  PN  NPH  6378     422.280 397.980 491.900  1.00  0.00
ATOM  12614 1OPN NPH  6378     422.690 397.240 493.120  1.00  0.00
ATOM  12615 2OPN NPH  6378     423.260 398.170 490.810  1.00  0.00
ATOM  12616  O3P NPH  6378     420.850 397.490 491.350  1.00  0.00
"""


# from $CCP4/examples/data/insulin.pdb
UNORDERED_ALTLOC_FRAGMENT = """\
ATOM     54  CB  THR A   8      21.486  49.557  34.680  1.00 17.33           C  
ANISOU   54  CB  THR A   8     2677   1737   2168   -147    348    -77       C  
ATOM     55  N  ASER A   9      22.340  46.718  36.111  0.50 15.13           N  
ANISOU   55  N  ASER A   9     2397   1624   1726    -74    278    -66       N  
ATOM     56  CA ASER A   9      22.426  45.878  37.287  0.50 14.83           C  
ANISOU   56  CA ASER A   9     2219   1729   1687    -83    202    -92       C  
ATOM     57  C  ASER A   9      23.045  44.544  36.849  0.50 13.99           C  
ANISOU   57  C  ASER A   9     2072   1701   1540   -107    221    -92       C  
ATOM     58  O  ASER A   9      23.264  44.329  35.667  0.50 13.72           O  
ANISOU   58  O  ASER A   9     2090   1725   1397   -111    331    -63       O  
ATOM     59  CB ASER A   9      21.034  45.674  37.848  0.50 15.28           C  
ANISOU   59  CB ASER A   9     2236   1828   1740      0    251   -242       C  
ATOM     60  OG ASER A   9      20.347  46.901  38.012  0.50 19.98           O  
ANISOU   60  OG ASER A   9     2799   2271   2521    150    548     90       O  
ATOM     61  N  BSER A   9      22.313  46.728  36.138  0.50 14.71           N  
ANISOU   61  N  BSER A   9     2322   1562   1703    -78    282    -86       N  
ATOM     62  CA BSER A   9      22.511  45.877  37.303  0.50 14.08           C  
ANISOU   62  CA BSER A   9     2078   1645   1625   -106    203   -118       C  
ATOM     63  C  BSER A   9      22.951  44.499  36.837  0.50 13.55           C  
ANISOU   63  C  BSER A   9     1971   1640   1538   -117    232   -119       C  
ATOM     64  O  BSER A   9      22.899  44.194  35.639  0.50 12.89           O  
ANISOU   64  O  BSER A   9     1807   1651   1438    -89    415    -90       O  
ATOM     65  CB BSER A   9      21.244  45.783  38.133  0.50 14.46           C  
ANISOU   65  CB BSER A   9     1987   1655   1851    -61    227   -231       C  
ATOM     66  OG BSER A   9      20.199  45.110  37.441  0.50 13.08           O  
ANISOU   66  OG BSER A   9     1812   1758   1398    -16    450     56       O  
ATOM     67  N   VAL A  10      23.342  43.662  37.798  1.00 13.44           N  
ANISOU   67  N   VAL A  10     2079   1653   1371   -123    133   -148       N  
"""  # noqa: W291 - trailing whitespace

# from https://cci.lbl.gov/hybrid_36/
HY36_EXAMPLE = """\
ATOM  99998  SD  MET L9999      48.231 -64.383  -9.257  1.00 11.54           S
ATOM  99999  CE  MET L9999      49.398 -63.242 -10.211  1.00 14.60           C
ATOM  A0000  N   VAL LA000      52.228 -67.689 -12.196  1.00  8.76           N
ATOM  A0001  CA  VAL LA000      53.657 -67.774 -12.458  1.00  3.40           C
"""

# from 1keb.pdb
TER_EXAMPLE = """\
ATOM   1643  OXT ALA B 108      -2.885   7.940  32.034  1.00 44.47           O  
TER    1644      ALA B 108                                                      
HETATM 1645 CU    CU A 109       3.422  -7.286  12.794  1.00 20.53          CU  
HETATM 1646 CU    CU B 109       3.446  18.122   8.201  1.00 15.05          CU  
"""  # noqa: W291 - trailing whitespace

def read_lines_and_remove(path):
    with open(path) as f:
        out_lines = f.readlines()
    os.remove(path)
    return out_lines


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
        block = gemmi.cif.read(full_path('5i55.cif'))[0]
        st = gemmi.make_structure_from_block(block)
        self.assertEqual(st.info['_entry.id'], '5I55')

        center = st[0].calculate_center_of_mass()
        # PyMOL>print cmd.centerofmass()
        pymol_ctr = [15.468438991742687, 4.8312495347721045, 20.607400844016833]
        self.assertTrue(center.dist(gemmi.Position(*pymol_ctr)) < 1e-7)

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

        output_block = st.make_mmcif_document().sole_block()
        cnames = block.get_mmcif_category_names()
        cnames_out = [name for name in output_block.get_mmcif_category_names()
                      if len(output_block.find_mmcif_category(name)) > 0]
        common_categories = [name for name in cnames_out if name in cnames]
        common_categories.sort()
        cc = ['_atom_site.', '_atom_type.', '_audit_author.',
              '_cell.', '_chem_comp.',
              '_diffrn.', '_diffrn_detector.', '_diffrn_radiation.',
              '_diffrn_source.',
              '_entity.', '_entity_poly.', '_entity_poly_seq.', '_entry.',
              '_exptl.', '_exptl_crystal.', '_pdbx_database_status.',
              '_pdbx_struct_assembly.', '_pdbx_struct_assembly_gen.',
              '_pdbx_struct_mod_residue.',
              '_pdbx_struct_oper_list.', '_refine.', '_reflns.', '_software.',
              '_struct.', '_struct_asym.', '_struct_conf.',
              '_struct_conf_type.', '_struct_conn.', '_struct_conn_type.',
              '_struct_keywords.', '_struct_ref.', '_struct_ref_seq.',
              '_symmetry.']
        self.assertEqual(common_categories, cc)
        for name in common_categories:
            cat_in = block.get_mmcif_category(name)
            cat_out = output_block.get_mmcif_category(name)
            for tag, values_out in cat_out.items():
                if tag == 'ccp4_link_id':
                    continue
                values_in = cat_in[tag]
                self.assertEqual(len(values_in), len(values_out))
                for (a, b) in zip(values_in, values_out):
                    try:
                        if a == b or abs(float(a) - float(b)) < 2e-4:
                            continue
                    except ValueError:
                        pass
                    self.assertTrue(name+tag in ['_struct_conf.id',
                                                 '_chem_comp.type'])
        for name_out in cnames_out:
            self.assertTrue(name_out in cnames)

    def test_5i55_predefined_removals(self, clear_entities=False):
        st = gemmi.read_structure(full_path('5i55.cif'))
        if clear_entities:
            self.assertEqual(len(st.entities), 4)
            st.entities = gemmi.EntityList()
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
        mse = model['A'][0]
        self.assertEqual(len(mse), 8)
        mse.trim_to_alanine()
        self.assertEqual([a.name for a in mse], ['N', 'CA', 'C', 'O', 'CB'])
        model['A'].trim_to_alanine()
        # ALA has 5 atoms, except the last one which has OXT (hence +1)
        expected_count = sum(4 + (r.name != 'GLY') for r in model['A']) + 1
        self.assertEqual(model.count_occupancies(), expected_count)
        st.remove_alternative_conformations()
        self.assertEqual(model.count_occupancies(), expected_count - 5/2.)
        self.assertTrue(not any(a.has_altloc() for a in lys12))

    def test_5i55_predefined_removals2(self):
        self.test_5i55_predefined_removals(clear_entities=True)

    def test_rnase_predefined_removals(self, add_entities=False):
        st = gemmi.read_structure(full_path('rnase_frag.pdb'))
        if add_entities:
            self.assertEqual(len(st.entities), 0)
            st.add_entity_types()
            st.assign_subchains()
            st.ensure_entities()
            self.assertEqual(len(st.entities), 4)
        model = st[0]
        nres_a = len(model['A'])
        nres_b = len(model['B'])
        st.add_entity_types()
        st.remove_ligands_and_waters()  # removes SO4 from each chain
        self.assertEqual(len(model['A']), nres_a - 1)
        self.assertEqual(len(model['B']), nres_b - 1)
        self.assertEqual(len(model['W']), 0)
        self.assertEqual(len(model), 3)
        st.remove_empty_chains()
        self.assertEqual(len(model), 2)

    def test_rnase_predefined_removals2(self):
        self.test_rnase_predefined_removals(add_entities=True)

    def test_3dg1(self):
        st = gemmi.read_structure(full_path('3dg1_final.cif'))
        self.assertEqual(st.info['_entry.id'], '3DG1')
        self.assertEqual(len(st[0]), 1)
        chain = st[0]['A']
        for res in chain[-2:]:
            self.assertEqual(res.name, 'HOH')
        for res in chain:
            for atom in res:
                n_images = st.cell.is_special_position(atom.pos)
                self.assertEqual(atom.occ * (n_images + 1), 1.0)

    def test_3wup(self):
        st = gemmi.read_structure(full_path('3wup.json.gz'))
        self.assertEqual(st.info['_entry.id'], '3WUP')
        model = st[0]
        self.assertEqual(len(model), 1)
        self.assertEqual(len(st.assemblies), 2)

        special_count = 0
        for res in st[0]['A'].get_ligands():
            for atom in res:
                n_images = st.cell.is_special_position(atom.pos)
                if n_images:
                    special_count += 1
                    self.assertEqual(round(atom.occ, 6), 0.33)
        self.assertEqual(special_count, 2)

        site_count = model.count_atom_sites()
        self.assertEqual(site_count, 279)
        how = gemmi.HowToNameCopiedChain.Short
        a1 = gemmi.make_assembly(st.assemblies[0], model, how)
        self.assertEqual(a1.count_atom_sites(), site_count)
        a2 = gemmi.make_assembly(st.assemblies[1], model, how)
        self.assertEqual(a2.count_atom_sites(), site_count * 3)
        gemmi.merge_atoms_in_expanded_model(a2, gemmi.UnitCell())
        # 3 atoms are on a 3-fold rotation axis
        self.assertEqual(a2.count_atom_sites(), (site_count - 3) * 3 + 3 * 1)

    def test_software_category(self):
        doc = gemmi.cif.read_file(full_path('3dg1_final.cif'))
        input_block = doc.sole_block()
        st = gemmi.make_structure_from_block(input_block)
        output_block = st.make_mmcif_document().sole_block()
        software = output_block.get_mmcif_category('_software')
        del software['date']
        assert software['version'][0] is False
        software['version'][0] = None
        self.assertEqual(input_block.get_mmcif_category('_software'), software)

    def test_5moo_header(self):
        st = gemmi.read_structure(full_path('5moo_header.pdb'))
        block = st.make_mmcif_document().sole_block()
        refine = block.get_mmcif_category('_refine')
        self.assertEqual(refine['ls_d_res_high'], ['1.44', '1.43'])
        self.assertEqual(refine['pdbx_starting_model'], ['4I8H'] * 2)
        self.assertEqual(list(block.find_values('_diffrn.ambient_temp')),
                         ['295', '295'])

    def check_1pfe(self, st):
        self.assertAlmostEqual(st.cell.a, 39.374)
        self.assertEqual(st.cell.gamma, 120)
        self.assertEqual(st.name, '1PFE')
        self.assertEqual(st.spacegroup_hm, 'P 63 2 2')
        self.assertEqual(st.info['_entry.id'], '1PFE')
        self.assertEqual(len(st), 1)
        self.assertEqual(len(st[0]), 2)
        label_to_auth_name = {sub.subchain_id(): ch.name for ch in st[0]
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
        atom_cl = res_cl.find_atom('CL', '*')
        self.assertAlmostEqual(atom_cl.occ, 0.17)
        self.assertEqual(atom_cl.element.name, 'Cl')
        return st

    def test_previous_next_residue(self):
        st = gemmi.read_structure(full_path('1pfe.cif.gz'),
                                  merge_chain_parts=False)
        chain_b = st[0]['B']
        res = chain_b['6']['ALA']
        res = chain_b.next_residue(res)
        self.assertEqual(res.name, 'NCY')
        res = chain_b.next_residue(res)
        self.assertEqual(res.name, 'MVA')
        self.assertIsNone(chain_b.next_residue(res))
        res = chain_b.previous_residue(res)
        self.assertEqual(res.name, 'NCY')
        res = chain_b.previous_residue(res)
        self.assertEqual(res.name, 'ALA')
        res = chain_b.next_residue(chain_b['7']['N2C'])
        self.assertEqual(res.name, 'MVA')
        self.assertEqual(chain_b.previous_residue(res).name, 'NCY')
        self.assertIsNone(chain_b.previous_residue(chain_b[0]))

    def test_read_1pfe_cif(self):
        st = gemmi.read_structure(full_path('1pfe.cif.gz'))
        self.check_1pfe(st)

        # write structure to cif and read it back
        out_name = get_path_for_tempfile(suffix='.cif')
        st.make_mmcif_document().write_file(out_name)
        st2 = gemmi.read_structure(out_name)
        os.remove(out_name)
        self.check_1pfe(st2)

    def test_read_1pfe_json(self):
        st = gemmi.read_structure(full_path('1pfe.json'))
        self.check_1pfe(st)

    def test_read_1orc(self):
        st = gemmi.read_structure(full_path('1orc.pdb'))
        self.assertEqual(st.resolution, 1.54)
        self.assertAlmostEqual(st.cell.a, 34.77)
        self.assertEqual(st.cell.alpha, 90)
        self.assertEqual(len(st.ncs), 0)
        self.assertEqual(st.meta.authors,
                         ['ALBRIGHT, R.A.', 'MOSSING, M.C.', 'MATTHEWS, B.W.'])
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

        result = gemmi.align_sequence_to_polymer(st.entities[0].full_sequence,
                                                 A.get_polymer(),
                                                 gemmi.PolymerType.Unknown)
        self.assertEqual(result.cigar_str(), '2I64M5I')

    def write_and_read(self, st, via_cif):
        if via_cif:
            st.setup_entities()
            st.assign_label_seq_id()
            doc = st.make_mmcif_document()
            st = gemmi.make_structure_from_block(doc[0])
        out_name = get_path_for_tempfile()
        st.write_pdb(out_name)
        return read_lines_and_remove(out_name)

    def test_read_write_1orc(self, via_cif=False):
        path = full_path('1orc.pdb')
        with open(path) as f:
            expected = [line for line in f if is_written_to_pdb(line, via_cif)]
        st = gemmi.read_structure(path)
        out_lines = self.write_and_read(st, via_cif)
        self.assertEqual(expected, out_lines)

    def test_read_write_1orc_via_cif(self):
        self.test_read_write_1orc(via_cif=True)

    def test_read_write_1lzh(self, via_cif=False):
        path = full_path('1lzh.pdb.gz')
        mode = 'rt' if sys.version_info >= (3,) else 'r'
        with gzip.open(path, mode=mode) as f:
            expected = [line for line in f if is_written_to_pdb(line, via_cif)]
        out_lines = self.write_and_read(gemmi.read_structure(path), via_cif)
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
            pa = ra['CA'][0].pos
            pb = rb['CA'][0].pos
            image_of_pb = st.ncs[0].apply(pb)
            self.assertTrue(pa.dist(image_of_pb) < 0.01)

    def test_read_write_5cvz_final(self, via_cif=False):
        path = full_path('5cvz_final.pdb')
        with open(path) as f:
            expected = [line.rstrip() for line in f
                        if is_written_to_pdb(line, via_cif)
                        # SCALE is not written b/c CRYST1 has more precision.
                        and line[:5] != 'SCALE']
        st = gemmi.read_structure(path)
        if via_cif:
            # input file w/o TER record -> subchains not setup automatically
            st.setup_entities()
            doc = st.make_mmcif_document()
            st = gemmi.make_structure_from_block(doc[0])
        out_lines = self.write_and_read(st, via_cif=False)
        if via_cif:
            out_lines = [line for line in out_lines
                         # input file has no REMARK 2, but it gets generated
                         # from REMARK 3 when going pdb->cif->pdb
                         if line[:10] != 'REMARK   2'
                         and line[:5] != 'TER  ']
        self.assertEqual(expected, [line.rstrip() for line in out_lines])

    def test_read_write_5cvz_final_via_cif(self):
        self.test_read_write_5cvz_final(via_cif=True)

    def test_read_write_4oz7(self, via_cif=False):
        path = full_path('4oz7.pdb')
        with open(path) as f:
            expected = [line for line in f if is_written_to_pdb(line, via_cif)]
        st = gemmi.read_structure(path, merge_chain_parts=False)
        out_lines = self.write_and_read(st, via_cif)
        self.assertEqual(expected, out_lines)

    def test_read_write_4oz7_via_cif(self):
        self.test_read_write_4oz7(via_cif=True)

    @unittest.skipIf(PDB is None, "BioPython not installed.")
    def test_reading_output_mmcif_with_biopython(self):
        path = full_path('4oz7.pdb')
        st = gemmi.read_structure(path)
        groups = gemmi.MmcifOutputGroups(True)
        # BioPython parser chokes without _atom_site.group_PDB
        # groups.group_pdb = False
        doc = st.make_mmcif_document(groups)
        parser = PDB.MMCIFParser(QUIET=True)
        parser.get_structure("none", StringIO(doc.as_string()))

    def test_pdb_element_names(self):
        pdb_line = "HETATM 4154 MG    MG A 341       1.384  19.340  11.968" \
                   "  1.00 67.64          MG"
        for line in [pdb_line, pdb_line.strip(' MG'), pdb_line[:-2] + '  ']:
            st = gemmi.read_pdb_string(line)
            residue = st[0].sole_residue('A', gemmi.SeqId(341, ' '))
            mg_atom = residue.sole_atom('MG')
            self.assertEqual(mg_atom.element.name, 'Mg')
            self.assertEqual(mg_atom.padded_name(), 'MG')
            self.assertAlmostEqual(mg_atom.b_iso, 67.64, delta=1e-6)
        mg_atom.element = gemmi.Element('Cu')
        self.assertEqual(mg_atom.element.name, 'Cu')

    def test_pdb_misaligned_element(self):
        pdb_line = "ATOM      7 S    SUB A   7      34.489 -14.293  34.343" \
                   "  0.29 43.77          S"
        for line in [pdb_line, pdb_line + '\n', pdb_line + '\r\n']:
            st = gemmi.read_pdb_string(line)
            atom = st[0].sole_residue('A', gemmi.SeqId('7')).sole_atom('S')
            self.assertEqual(atom.element.name, 'S')

    def test_pdb_element_names_from_amber(self):
        st = gemmi.read_pdb_string(AMBER_FRAGMENT)
        residue = st[0][''][0]
        self.assertEqual(residue.sole_atom('CB').element, gemmi.Element('C'))
        self.assertEqual(residue.sole_atom('HB').element, gemmi.Element('H'))
        self.assertEqual(residue.sole_atom('CG1').element, gemmi.Element('C'))
        self.assertEqual(residue.sole_atom('HG11').element, gemmi.Element('H'))
        lines = AMBER_FRAGMENT.splitlines()
        for n, atom in enumerate(residue):
            self.assertEqual(atom.padded_name(), lines[n][12:16].rstrip())
        chain = gemmi.read_pdb_string(FRAGMENT_WITH_HG)[0]['P']
        self.assertEqual(chain[0].sole_atom('HG').element, gemmi.Element('Hg'))
        self.assertEqual(chain[1].sole_atom('HG1').element, gemmi.Element('Hg'))

    def test_pdb_element_names_from_trjconv(self):
        st = gemmi.read_pdb_string(TRJCONV_FRAGMENT)
        residue = st[0][''][0]
        expected = ['C', 'H', 'H', 'O', 'P', 'O', 'O', 'O']
        lines = TRJCONV_FRAGMENT.splitlines()
        for n, atom in enumerate(residue):
            self.assertEqual(atom.element.name, expected[n])
            self.assertEqual(atom.padded_name(), lines[n][12:16].rstrip())

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

    def test_blank_mmcif(self):
        input_block = gemmi.cif.Block('empty')
        st = gemmi.make_structure_from_block(input_block)
        self.assertEqual(st.name, 'empty')
        output_block = st.make_mmcif_document().sole_block()
        self.assertEqual(output_block.get_mmcif_category_names(), [])

    def test_blank_chain(self):
        st = gemmi.read_pdb_string(BLANK_CHAIN_FRAGMENT)
        out_name = get_path_for_tempfile()
        st.write_minimal_pdb(out_name)
        out = read_lines_and_remove(out_name)
        # CRYST1 differs (50.000 not 50.00 and added P1).
        # ATOM lines have added element names.
        trimmed_out = [line[:66] for line in out[1:]]
        self.assertEqual(trimmed_out, BLANK_CHAIN_FRAGMENT.splitlines()[1:])

    def test_ncs(self):
        st = gemmi.read_structure(full_path('5cvz_final.pdb'))
        self.assertEqual(len(st.cell.images), 20 * 12 - 1)
        self.assertEqual(st.resolution, 3.29)
        first_atom = st[0].sole_residue('A', gemmi.SeqId(17, ' '))[0]
        ne2 = st[0].sole_residue('A', gemmi.SeqId('63')).sole_atom('NE2')
        direct_dist = first_atom.pos.dist(ne2.pos)
        self.assertAlmostEqual(direct_dist, 34.89, delta=1e-2)
        nearest_image = st.cell.find_nearest_image(first_atom.pos, ne2.pos)
        nearest_dist = nearest_image.dist()
        self.assertAlmostEqual(nearest_dist, 8.02, delta=1e-2)

        # test __getitem__(splice) - unrelated to NCS (sometimes we put
        # unrelated tests together to avoid the same file again)
        chain = st[0]['A']
        res = chain.next_residue(chain[:5][0])
        self.assertTrue(res is chain[1])

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

    def test_short_ssbond(self):
        st = gemmi.read_pdb_string(SHORT_SSBOND)
        out = st.make_pdb_headers()
        self.assertEqual(out.splitlines()[0],
                         SHORT_SSBOND.splitlines()[0]
                         + "                          1555   1555  2.06  ")

    def test_add_remove(self):
        st = gemmi.read_pdb_string(SSBOND_FRAGMENT)
        st.add_model(st[0])
        st.renumber_models()
        res = st[0].sole_residue('A', gemmi.SeqId('310'))
        self.assertEqual(len(res), 1)
        res.remove_atom('SG', ' ')
        self.assertEqual(len(res), 0)
        res = st[1].sole_residue('A', gemmi.SeqId('310'))
        self.assertEqual(len(res), 1)
        self.assertEqual(len(st[0]['A']), 7)
        del st[0]['A'][3]
        self.assertEqual(len(st[0]['A']), 6)
        self.assertEqual(len(st), 2)
        self.assertEqual(st[0].name, '1')
        del st['1']
        self.assertEqual(len(st), 1)
        self.assertEqual(st[0].name, '2')
        st.renumber_models()
        self.assertEqual(st[0].name, '1')
        st.add_model(st[0])
        st.add_model(st[0])
        st.renumber_models()
        self.assertEqual(st[0].name, '1')
        self.assertEqual(st[-1].name, '3')
        del st[:-1]
        self.assertEqual(len(st), 1)
        self.assertEqual(st[0].name, '3')
        del st[0]
        self.assertEqual(len(st), 0)

    def test_remove2(self):
        # test also save_doc
        saved_doc = gemmi.cif.Document()
        st = gemmi.read_structure(full_path('1pfe.cif.gz'), save_doc=saved_doc)
        self.assertEqual(saved_doc[0].name, '1PFE')
        self.assertEqual(saved_doc[0].find_value('_entity_name_com.name'),
                         "'QUINOMYCIN A'")
        model = st[0]
        self.assertEqual(len(model), 2)
        b = model['B']
        self.assertEqual(b[0].name, 'DSN')
        del b['1']['DSN']
        self.assertEqual(b[0].name, 'ALA')
        del b[0]
        self.assertEqual(b[0].name, 'N2C')

        # test append_residues()
        self.assertEqual(len(b), 20)
        b.append_residues(b[:5], min_sep=10)
        self.assertEqual(len(b), 25)

        # test append_residues() with empty chain
        new_chain = gemmi.Chain('X')
        new_chain.append_residues(b[:5], min_sep=1)
        self.assertEqual(len(new_chain), 5)

        # test adding and removing chains
        model.add_chain(new_chain, unique_name=False)
        model.add_chain(new_chain, unique_name=True)
        model.add_chain(new_chain)
        self.assertEqual([chain.name for chain in model], list('ABXCX'))
        del model[2:]
        model.add_chain(new_chain, unique_name=True)
        self.assertEqual([chain.name for chain in model], list('ABX'))
        del model[-1]
        del model['A']
        self.assertEqual(len(model), 1)
        self.assertEqual(model[0].name, 'B')
        doc = st.make_mmcif_document()
        ref_seq = doc[0].get_mmcif_category('_struct_ref_seq')
        self.assertEqual(ref_seq['pdbx_strand_id'], ['B'])

    def test_first_conformer(self):
        model = gemmi.read_structure(full_path('1pfe.cif.gz'))[0]
        b = model['B']
        self.assertEqual([res.name for res in b if not res.is_water()],
                         ['DSN', 'ALA', 'N2C', 'NCY', 'MVA', 'DSN',
                          'ALA', 'NCY', 'N2C', 'MVA', 'QUI', 'QUI'])
        self.assertEqual([res.name for res in b.first_conformer()
                          if not res.is_water()],
                         ['DSN', 'ALA', 'N2C', 'MVA', 'DSN',
                          'ALA', 'NCY', 'MVA', 'QUI', 'QUI'])
        polymer = b.get_polymer()
        self.assertEqual([res.name for res in polymer],
                         ['DSN', 'ALA', 'N2C', 'NCY', 'MVA', 'DSN',
                          'ALA', 'NCY', 'N2C', 'MVA'])
        self.assertEqual([res.name for res in polymer.first_conformer()],
                         ['DSN', 'ALA', 'N2C', 'MVA', 'DSN',
                          'ALA', 'NCY', 'MVA'])
        self.assertEqual(len(polymer), 10)
        self.assertEqual(polymer.length(), 8)
        # The bond between 4 MVA(v) and 5 DSN(s) is between C and OG
        # (so it's not a peptide bond as expected). Depending on the heuristic
        # used to determine gaps, this sequence could have a gap in the middle.
        self.assertEqual(polymer.make_one_letter_sequence(), 'sAXvsAXv')
        self.assertEqual([res.name for res in b.get_ligands()], ['QUI', 'QUI'])
        res1 = model.sole_residue('A', gemmi.SeqId('1'))
        self.assertEqual([atom.name for atom in res1.first_conformer()],
                         [atom.name for atom in res1 if atom.altloc != 'B'])

    def test_model_all(self):
        model = gemmi.Model('1')
        for name in 'ABCDEFG':
            model.add_chain(gemmi.Chain(name))
        expected = []
        for cname in 'BCF':
            chain = model[cname]
            for _ in range(7):
                chain.add_residue(gemmi.Residue())
            for (r, name) in [(2, '0'), (2, '1'), (3, '2'), (5, '3')]:
                a = gemmi.Atom()
                a.name = cname + name
                expected.append(a.name)
                chain[r].add_atom(a)
        self.assertEqual([cra.atom.name for cra in model.all()], expected)
        st = gemmi.Structure()
        st.add_model(model)
        st.remove_empty_chains()
        self.assertEqual([cra.atom.name for cra in model.all()], expected)

    def test_different_altloc_order(self):
        st = gemmi.read_pdb_string(UNORDERED_ALTLOC_FRAGMENT)
        chain = st[0]['A']
        cb = chain['9']['SER']['CB']
        cb_numbers = [atom.serial for atom in cb]
        self.assertEqual(cb_numbers, [59, 65])
        self.assertEqual(cb[0].serial, 59)
        self.assertEqual(cb[1].serial, 65)
        self.assertEqual(cb[-1].serial, 65)
        self.assertEqual(cb[-2].serial, 59)
        self.assertEqual(chain.count_atom_sites(), 14)
        self.assertEqual(chain.count_occupancies(), 8)
        st.remove_alternative_conformations()
        self.assertEqual(chain.count_atom_sites(), 8)

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

    def test_assembly(self):
        st = gemmi.read_structure(full_path('1pfe.cif.gz'),
                                  merge_chain_parts=False)
        model = st[0]
        ch_names = [ch.name for ch in model]
        self.assertEqual(['A', 'B', 'A', 'B', 'A', 'B'], ch_names)
        model_mass = model.calculate_mass()

        self.assertEqual(len(st.assemblies), 1)
        asem = st.assemblies[0]
        bio = gemmi.make_assembly(asem, model,
                                  gemmi.HowToNameCopiedChain.Short)
        new_naming = {'A': 'C', 'B': 'D'}
        self.assertEqual([ch.name for ch in bio],
                         ch_names + [new_naming[x] for x in ch_names])
        self.assertAlmostEqual(bio.calculate_mass(), 2 * model_mass)

        bio = gemmi.make_assembly(asem, model,
                                  gemmi.HowToNameCopiedChain.AddNumber)
        self.assertEqual([ch.name for ch in bio],
                         [x+'1' for x in ch_names] + [x+'2' for x in ch_names])

    def test_assembly_naming(self):
        st = gemmi.read_structure(full_path('4oz7.pdb'))
        model = st[0]
        a1 = st.assemblies[1]
        bio = gemmi.make_assembly(a1, model,
                                  gemmi.HowToNameCopiedChain.AddNumber)
        self.assertEqual([ch.name for ch in bio], ['B1'])
        bio = gemmi.make_assembly(a1, model, gemmi.HowToNameCopiedChain.Dup)
        self.assertEqual([ch.name for ch in bio], ['B'])
        bio = gemmi.make_assembly(a1, model, gemmi.HowToNameCopiedChain.Short)
        self.assertEqual([ch.name for ch in bio], ['B'])

    def test_hybrid36(self):
        st = gemmi.read_pdb_string(HY36_EXAMPLE)
        nums = [(cra.atom.serial, cra.residue.seqid.num) for cra in st[0].all()]
        self.assertEqual(nums, [(99998, 9999), (99999, 9999),
                                (100000, 10000), (100001, 10000)])
        minimal_opt = gemmi.PdbWriteOptions(minimal=True)
        out = st.make_pdb_string(minimal_opt).splitlines()
        # original serial numbers are lost when writing pdb, check only seqid
        self.assertEqual(out[1][17:30], 'MET L9999    ')
        self.assertEqual(out[4][17:30], 'VAL LA000    ')

    def test_ter_writing(self):
        st = gemmi.read_pdb_string(TER_EXAMPLE)
        def write_and_read(**kwargs):
            opt = gemmi.PdbWriteOptions(**kwargs)
            kwargs['minimal'] = True
            st_out = gemmi.read_pdb_string(st.make_pdb_string(opt))
            return [(cra.chain.name + ' ' + cra.atom.name, cra.atom.serial)
                    for cra in st_out[0].all()]
        self.assertEqual(write_and_read(),
                         [('B OXT', 1), ('A CU', 3), ('B CU', 4)])
        st.merge_chain_parts()
        self.assertEqual(write_and_read(),
                         [('B OXT', 1), ('B CU', 3), ('A CU', 4)])
        self.assertEqual(write_and_read(numbered_ter=False),
                         [('B OXT', 1), ('B CU', 2), ('A CU', 3)])
        self.assertEqual(write_and_read(preserve_serial=True),
                         [('B OXT', 1643), ('B CU', 1646), ('A CU', 1645)])
        st.assign_serial_numbers(numbered_ter=True)
        self.assertEqual(write_and_read(preserve_serial=True),
                         [('B OXT', 1), ('B CU', 3), ('A CU', 4)])

if __name__ == '__main__':
    unittest.main()
