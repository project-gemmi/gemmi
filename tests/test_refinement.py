import math
import unittest

import gemmi


ATOM_SITE = """
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_asym_id
1 C CA . ALA A 1 0 0 0 1 A
"""


SCALAR_FIELDS = [
    ('pdbx_ls_cross_valid_method', 'cross_validation_method', 'THROUGHOUT'),
    ('pdbx_R_Free_selection_details', 'rfree_selection_method', 'random set'),
    ('B_iso_mean', 'mean_b', 12.5),
    ('pdbx_overall_SU_R_Blow_DPI', 'dpi_blow_r', 0.125),
    ('pdbx_overall_SU_R_free_Blow_DPI', 'dpi_blow_rfree', 0.25),
    ('overall_SU_R_Cruickshank_DPI', 'dpi_cruickshank_r', 0.375),
    ('pdbx_overall_SU_R_free_Cruickshank_DPI', 'dpi_cruickshank_rfree', 0.5),
    ('correlation_coeff_Fo_to_Fc', 'cc_fo_fc_work', 0.75),
    ('correlation_coeff_Fo_to_Fc_free', 'cc_fo_fc_free', 0.625),
    ('pdbx_average_fsc_work', 'fsc_work', 0.5),
    ('pdbx_average_fsc_free', 'fsc_free', -0.25),
    ('correlation_coeff_I_to_Fcsqd_work', 'cc_intensity_work', 0.0),
    ('correlation_coeff_I_to_Fcsqd_free', 'cc_intensity_free', -1.0),
]
ANISO_FIELDS = [
    ('aniso_B[1][1]', 'u11', 11.125),
    ('aniso_B[2][2]', 'u22', 22.25),
    ('aniso_B[3][3]', 'u33', 33.5),
    ('aniso_B[1][2]', 'u12', -1.25),
    ('aniso_B[1][3]', 'u13', -2.5),
    ('aniso_B[2][3]', 'u23', -3.75),
]


class TestRefinementMetadata(unittest.TestCase):
    def test_pairs_roundtrip(self):
        text = "data_refinement\n_refine.pdbx_refine_id xray\n"
        for tag, _, value in SCALAR_FIELDS + ANISO_FIELDS:
            text += '_refine.' + tag + ' ' + gemmi.cif.quote(str(value)) + '\n'
        text += '_refine.ls_d_res_high 1.25\n_refine.ls_R_factor_R_work 0.1875\n'
        st = gemmi.make_structure_from_block(gemmi.cif.read_string(text + ATOM_SITE)[0])
        self.assertEqual(st.resolution, 1.25)
        self.assertEqual(st.meta.refinement[0].r_work, 0.1875)
        block = st.make_mmcif_document()[0]
        for tag, attr, value in SCALAR_FIELDS:
            with self.subTest(tag=tag):
                self.assertEqual(getattr(st.meta.refinement[0], attr), value)
                actual = block.find_value('_refine.' + tag)
                self.assertIsNotNone(actual)
                if isinstance(value, str):
                    self.assertEqual(gemmi.cif.as_string(actual), value)
                else:
                    self.assertEqual(float(actual), value)
        for tag, attr, value in ANISO_FIELDS:
            with self.subTest(tag=tag):
                self.assertEqual(getattr(st.meta.refinement[0].aniso_b, attr),
                                 value)
                self.assertEqual(float(block.find_value('_refine.' + tag)),
                                 value)
        st2 = gemmi.make_structure_from_block(block)
        self.assertEqual(st2.meta.refinement[0].aniso_b.elements_pdb(),
                         [value for _, _, value in ANISO_FIELDS])

    def test_loop_roundtrip(self):
        text = """data_refinement
loop_
_refine.pdbx_refine_id
_refine.pdbx_ls_cross_valid_method
_refine.pdbx_R_Free_selection_details
_refine.B_iso_mean
_refine.aniso_B[1][1]
_refine.aniso_B[2][2]
_refine.aniso_B[1][3]
_refine.pdbx_overall_SU_R_Blow_DPI
_refine.pdbx_average_fsc_free
'X-RAY DIFFRACTION' ? 'random set' 0 1.25 ? -2.5 0 .
'NEUTRON DIFFRACTION' '.' ? . ? 3.5 . ? -0.25
"""
        st = gemmi.make_structure_from_block(gemmi.cif.read_string(text + ATOM_SITE)[0])
        xray, neutron = st.meta.refinement
        self.assertEqual([xray.id, neutron.id],
                         ['X-RAY DIFFRACTION', 'NEUTRON DIFFRACTION'])
        self.assertEqual(xray.cross_validation_method, '')
        self.assertEqual(neutron.cross_validation_method, '.')
        self.assertEqual(xray.rfree_selection_method, 'random set')
        self.assertEqual(neutron.rfree_selection_method, '')
        self.assertEqual(xray.mean_b, 0)
        self.assertTrue(math.isnan(neutron.mean_b))
        self.assertEqual(xray.aniso_b.u11, 1.25)
        self.assertTrue(math.isnan(xray.aniso_b.u22))
        self.assertTrue(math.isnan(xray.aniso_b.u33))
        self.assertEqual(xray.aniso_b.u13, -2.5)
        self.assertTrue(math.isnan(neutron.aniso_b.u11))
        self.assertEqual(neutron.aniso_b.u22, 3.5)
        self.assertEqual(xray.dpi_blow_r, 0)
        self.assertTrue(math.isnan(neutron.dpi_blow_r))
        self.assertTrue(math.isnan(xray.fsc_free))
        self.assertEqual(neutron.fsc_free, -0.25)
        block = st.make_mmcif_document()[0]
        self.assertEqual([gemmi.cif.as_string(x) for x in
                          block.find_values('_refine.pdbx_refine_id')],
                         ['X-RAY DIFFRACTION', 'NEUTRON DIFFRACTION'])
        category = block.get_mmcif_category('_refine.')
        self.assertEqual(category['pdbx_ls_cross_valid_method'], [None, '.'])
        self.assertEqual(category['pdbx_R_Free_selection_details'],
                         ['random set', None])
        self.assertEqual(category['B_iso_mean'], ['0', None])
        self.assertEqual(category['aniso_B[1][1]'], ['1.25', None])
        self.assertEqual(category['aniso_B[2][2]'], [None, '3.5'])
        self.assertEqual(category['aniso_B[1][3]'], ['-2.5', None])
        self.assertEqual(category['pdbx_average_fsc_free'], [None, '-0.25'])

    def test_optional_columns(self):
        for fields, matrix in [(SCALAR_FIELDS, False), (ANISO_FIELDS, True)]:
            for tag, attr, value in fields:
                with self.subTest(tag=tag):
                    text = ('data_refinement\n_refine.pdbx_refine_id xray\n'
                            '_refine.' + tag + ' '
                            + gemmi.cif.quote(str(value)) + '\n')
                    st = gemmi.make_structure_from_block(
                        gemmi.cif.read_string(text + ATOM_SITE)[0])
                    ref = st.meta.refinement[0]
                    self.assertEqual(getattr(ref.aniso_b if matrix else ref,
                                             attr), value)
                    if matrix:
                        for _, other, _ in ANISO_FIELDS:
                            if other != attr:
                                self.assertTrue(math.isnan(
                                    getattr(ref.aniso_b, other)))

    def test_null_and_absent_columns(self):
        for null in ['?', '.']:
            for loop in [False, True]:
                with self.subTest(null=null, loop=loop):
                    tags = ['pdbx_refine_id']
                    tags += [tag for tag, _, _ in SCALAR_FIELDS + ANISO_FIELDS]
                    if loop:
                        text = ('data_refinement\nloop_\n_refine.'
                                + '\n_refine.'.join(tags) + '\nxray '
                                + ' '.join([null] * (len(tags) - 1)) + '\n')
                    else:
                        text = 'data_refinement\n_refine.pdbx_refine_id xray\n'
                        text += ''.join('_refine.' + tag + ' ' + null + '\n'
                                        for tag in tags[1:])
                    st = gemmi.make_structure_from_block(
                        gemmi.cif.read_string(text + ATOM_SITE)[0])
                    ref = st.meta.refinement[0]
                    for _, attr, value in SCALAR_FIELDS:
                        if isinstance(value, str):
                            self.assertEqual(getattr(ref, attr), '')
                        else:
                            self.assertTrue(math.isnan(getattr(ref, attr)))
                    self.assertTrue(all(math.isnan(x)
                                        for x in ref.aniso_b.elements_pdb()))
                    out = st.make_mmcif_document()[0]
                    self.assertEqual(list(out.find_values('_refine.B_iso_mean')),
                                     [])
        st = gemmi.make_structure_from_block(gemmi.cif.read_string(
            'data_refinement\n_refine.pdbx_refine_id xray\n' + ATOM_SITE)[0])
        self.assertEqual(len(st.meta.refinement), 1)
        self.assertEqual(st.meta.refinement[0].id, 'xray')
        self.assertTrue(math.isnan(st.meta.refinement[0].mean_b))


if __name__ == '__main__':
    unittest.main()
