#!/usr/bin/env python

import os
import unittest
import gemmi
from common import full_path, get_path_for_tempfile, assert_numpy_equal, numpy

def compare_maps(self, a, b, atol):
    #print(abs(numpy.array(a) - b).max())
    self.assertTrue(numpy.allclose(a, b, atol=atol, rtol=0))

def compare_asu_data(self, asu_data, data, f, phi):
    asu_dict = {tuple(a.hkl): a.value for a in asu_data}
    data_hkl = data.make_miller_array()
    self.assertEqual(data_hkl.dtype, 'int32')
    if type(data) is gemmi.ReflnBlock:
        data_f = data.make_float_array(f)
        data_phi = data.make_float_array(phi)
    else:
        data_f = data.column_with_label(f).array
        data_phi = data.column_with_label(phi).array
    self.assertTrue(len(data_f) > 100)
    asu_val = numpy.array([asu_dict[tuple(hkl)] for hkl in data_hkl])
    self.assertTrue(numpy.allclose(data_f, abs(asu_val), atol=5e-5, rtol=0))
    x180 = abs(180 - abs(data_phi - numpy.angle(asu_val, deg=True)))
    self.assertTrue(numpy.allclose(x180, 180, atol=6e-2, rtol=0.))

def fft_test(self, data, f, phi, size, order=gemmi.AxisOrder.XYZ):
    if numpy is None:
        return
    self.assertTrue(data.data_fits_into(size))
    grid_full = data.get_f_phi_on_grid(f, phi, size, half_l=False, order=order)
    self.assertEqual(grid_full.axis_order, order)
    array_full = numpy.array(grid_full, copy=False)
    map1 = gemmi.transform_f_phi_grid_to_map(grid_full)
    self.assertEqual(map1.axis_order, order)
    map2 = numpy.fft.ifftn(array_full.conj())
    map2 = numpy.real(map2) * (map2.size / grid_full.unit_cell.volume)
    compare_maps(self, map1, map2, atol=6e-7)
    map3 = data.transform_f_phi_to_map(f, phi, size, order=order)
    compare_maps(self, map1, map3, atol=6e-7)

    grid2 = gemmi.transform_map_to_f_phi(map1, half_l=False)
    self.assertFalse(grid2.half_l)
    self.assertEqual(grid2.axis_order, order)
    compare_maps(self, grid2, array_full, atol=2e-4)
    if grid2.axis_order != gemmi.AxisOrder.ZYX:
        compare_asu_data(self, grid2.prepare_asu_data(), data, f, phi)
    grid_half = data.get_f_phi_on_grid(f, phi, size, half_l=True, order=order)
    if order == gemmi.AxisOrder.ZYX:  # half_l+ZYX not supported yet
        return
    grid3 = gemmi.transform_map_to_f_phi(map1, half_l=True)
    self.assertTrue(grid3.half_l)
    self.assertEqual(grid3.axis_order, order)
    compare_maps(self, grid3, grid_half, atol=2e-4)
    compare_asu_data(self, grid3.prepare_asu_data(), data, f, phi)

    asu_data = grid_full.prepare_asu_data()
    back_grid = asu_data.get_f_phi_on_grid(size, half_l=False, order=order)
    compare_maps(self, back_grid, grid_full, atol=4e-5)

    asu_data = grid3.prepare_asu_data()
    back_grid = asu_data.get_f_phi_on_grid(size, half_l=True, order=order)
    compare_maps(self, back_grid, grid3, atol=5e-5)

    asu_data = data.get_f_phi(f, phi)
    compare_asu_data(self, asu_data, data, f, phi)
    hkl_grid = asu_data.get_f_phi_on_grid(size, order=order)
    compare_maps(self, hkl_grid, grid_full, atol=1e-4)
    if order != gemmi.AxisOrder.ZYX:
        m1 = gemmi.transform_f_phi_grid_to_map(hkl_grid)
        m2 = asu_data.transform_f_phi_to_map(exact_size=size)
        compare_maps(self, m1, m2, atol=1e-6)


class TestMtz(unittest.TestCase):
    def test_read_write(self):
        path = full_path('5e5z.mtz')
        mtz = gemmi.read_mtz_file(path)
        self.assertEqual(mtz.spacegroup.hm, 'P 1 21 1')
        out_name = get_path_for_tempfile()
        mtz.write_to_file(out_name)
        mtz2 = gemmi.read_mtz_file(out_name)
        os.remove(out_name)
        self.assertEqual(mtz2.spacegroup.hm, 'P 1 21 1')
        if numpy is not None:
            assert_numpy_equal(self, numpy.array(mtz, copy=False), mtz.array)
            assert_numpy_equal(self, mtz.array, mtz2.array)

    def test_remove_and_add_column(self):
        path = full_path('5e5z.mtz')
        col_name = 'FREE'
        mtz = gemmi.read_mtz_file(path)
        col = mtz.column_with_label(col_name)
        col_idx = col.idx
        ncol = len(mtz.columns)
        if numpy is None:
            return
        assert_numpy_equal(self, col.array, numpy.array(col, copy=False))
        arr = col.array.copy()
        mtz_data = numpy.array(mtz, copy=True)
        self.assertEqual(mtz_data.shape, (mtz.nreflections, ncol))
        mtz.remove_column(col_idx)
        self.assertEqual(len(mtz.columns), ncol-1)
        self.assertEqual(numpy.array(mtz, copy=False).shape,
                         (mtz.nreflections, ncol-1))
        col = mtz.add_column(col_name, 'I', dataset_id=0, pos=col_idx)
        numpy.array(col, copy=False)[:] = arr
        assert_numpy_equal(self, mtz_data, numpy.array(mtz, copy=False))

    def asu_data_test(self, grid):
        asu = grid.prepare_asu_data()
        d = asu.make_d_array()
        asu2 = gemmi.ComplexAsuData(asu.unit_cell,
                                    asu.spacegroup,
                                    asu.miller_array[d > 2.5],
                                    asu.value_array[d > 2.5])
        ngrid = asu2.transform_f_phi_to_map(sample_rate=1.5)
        hkl_grid = asu2.get_f_phi_on_grid([ngrid.nu, ngrid.nv, ngrid.nw])
        alt_ngrid = gemmi.transform_f_phi_grid_to_map(hkl_grid)
        compare_maps(self, ngrid, alt_ngrid, atol=1e-6)

    def test_f_phi_grid(self):
        path = full_path('5wkd_phases.mtz.gz')
        mtz = gemmi.read_mtz_file(path)
        size = mtz.get_size_for_hkl()
        for half_l in (False, True):
            grid1 = mtz.get_f_phi_on_grid('FWT', 'PHWT', size, half_l=half_l)
            grid2 = mtz.get_f_phi_on_grid('FWT', 'PHWT', size, half_l=half_l,
                                          order=gemmi.AxisOrder.ZYX)
            if numpy is None:
                continue
            array1 = numpy.array(grid1, copy=False)
            array2 = numpy.array(grid2, copy=False)
            self.assertTrue((array2 == array1.transpose(2,1,0)).all())
            self.asu_data_test(grid1)

        fft_test(self, mtz, 'FWT', 'PHWT', size)
        fft_test(self, mtz, 'FWT', 'PHWT', size, order=gemmi.AxisOrder.ZYX)

    def test_value_grid(self):
        #path = full_path('5wkd_phases.mtz.gz')
        path = full_path('5e5z.mtz')
        mtz = gemmi.read_mtz_file(path)
        size = mtz.get_size_for_hkl()
        if numpy is None:
            return
        asu = gemmi.ReciprocalAsu(mtz.spacegroup)
        mtz_data = numpy.array(mtz, copy=False)
        fp_idx = mtz.column_labels().index('FP')
        fp_map = {}
        for row in mtz_data:
            fp_map[tuple(row[0:3])] = row[fp_idx]
        for order in (gemmi.AxisOrder.XYZ, gemmi.AxisOrder.ZYX):
            for half_l in (True, False):
                grid = mtz.get_value_on_grid('FP', size,
                                             half_l=half_l, order=order)
                counter = 0
                for point in grid:
                    hkl = grid.to_hkl(point)
                    value = fp_map.get(tuple(hkl))
                    if asu.is_in(hkl):
                        if value is not None:
                            self.assertTrue(point.value == value
                                            or (numpy.isnan(point.value)
                                                and numpy.isnan(value)))
                            counter += 1
                        else:
                            self.assertEqual(point.value, 0.)
                    else:
                        self.assertIsNone(value)
                self.assertEqual(counter, mtz_data.shape[0])


class TestSfMmcif(unittest.TestCase):
    def test_reading(self):
        doc = gemmi.cif.read(full_path('r5wkdsf.ent'))
        rblock = gemmi.as_refln_blocks(doc)[0]
        self.assertEqual(rblock.spacegroup.hm, 'C 1 2 1')

        size = rblock.get_size_for_hkl()
        for order in (gemmi.AxisOrder.XYZ, gemmi.AxisOrder.ZYX):
            fft_test(self, rblock, 'pdbx_FWT', 'pdbx_PHWT', size)

    def test_scaling(self):
        doc = gemmi.cif.read(full_path('r5wkdsf.ent'))
        rblock = gemmi.as_refln_blocks(doc)[0]
        fobs_data = rblock.get_value_sigma('F_meas_au', 'F_meas_sigma_au')
        if numpy:
            self.assertEqual(fobs_data.value_array.shape, (367,))

        # without mask
        fc_data = rblock.get_f_phi('F_calc_au', 'phase_calc')
        scaling = gemmi.Scaling(fc_data.unit_cell, fc_data.spacegroup)
        scaling.prepare_points(fc_data, fobs_data)
        scaling.fit_isotropic_b_approximately()
        scaling.fit_parameters()
        #print(scaling.k_overall, scaling.b_overall)
        scaling.scale_data(fc_data)

        # with mask
        # TODO
        #fc_data = rblock.get_f_phi('F_calc_au', 'phase_calc')
        #st = gemmi.read_structure(full_path('5wkd.pdb'))
        #fmask = ...

class TestBinner(unittest.TestCase):
    def test_binner(self):
        path = full_path('5e5z.mtz')
        def check_limits_17(limits):
            self.assertEqual(len(limits), 17)
            self.assertAlmostEqual(limits[10], 0.27026277234462415)
        mtz = gemmi.read_mtz_file(path)
        binner = gemmi.Binner()
        method = gemmi.Binner.Method.Dstar3
        binner.setup(17, method, mtz)
        check_limits_17(binner.limits)
        self.assertEqual(binner.size, 17)
        binner = gemmi.Binner()
        binner.setup(17, method, mtz, cell=mtz.cell)
        check_limits_17(binner.limits)
        self.assertEqual(binner.get_bin([3,3,3]), 9)
        if numpy is None:
            return
        binner.setup(17, method, mtz.make_miller_array(), cell=mtz.cell)
        check_limits_17(binner.limits)
        binner.setup_from_1_d2(17, method, mtz.make_1_d2_array(), mtz.cell)
        check_limits_17(binner.limits)
        hkls = [[0,0,1], [3,3,3], [10,10,10]]
        bins = [0,9,16]
        self.assertEqual(list(binner.get_bins(hkls)), bins)
        inv_d2 = [mtz.cell.calculate_1_d2(h) for h in hkls]
        self.assertEqual(list(binner.get_bins_from_1_d2(inv_d2)), bins)

class TestConversion(unittest.TestCase):
    def test_4aap(self):
        def check_metadata(o, d):
            self.assertEqual(o.spacegroup.hm, 'P 32 2 1')
            self.assertEqual(o.cell.a, 68.575)
            self.assertEqual(d.cell.gamma, 120)
            self.assertEqual(d.wavelength, 0.97)
        doc = gemmi.cif.read(full_path('4aap-sf-subset.cif'))
        (rblock,) = gemmi.as_refln_blocks(doc)
        check_metadata(rblock, rblock)
        self.assertFalse(rblock.is_unmerged())
        mtz = gemmi.CifToMtz().convert_block_to_mtz(rblock)
        check_metadata(mtz, mtz.datasets[1])
        cif_string = gemmi.MtzToCif().write_cif_to_string(mtz)
        doc_out = gemmi.cif.read_string(cif_string)
        (rblock_out,) = gemmi.as_refln_blocks(doc_out)
        self.assertEqual(rblock.default_loop.tags[3:],
                         rblock_out.default_loop.tags)
        check_metadata(rblock_out, rblock_out)
        def cif_floats(rb, tag):
            return [float(x) if x != '?' else 8008.8008  # (arbitrary number)
                    for x in rb.block.find_values(tag)]
        tag = '_refln.phase_calc'
        self.assertEqual(cif_floats(rblock, tag), cif_floats(rblock_out, tag))
        tag = '_refln.pdbx_HL_A_iso'
        self.assertEqual(cif_floats(rblock, tag), cif_floats(rblock_out, tag))
        d1 = mtz.row_as_dict((7, 5, -1))
        self.assertEqual(d1['K'], 5)  # just to check row_as_dict()
        self.assertAlmostEqual(d1['FC'], 22.7, delta=1e-5)
        self.assertAlmostEqual(d1['PHIC'], 91.1, delta=1e-5)

        # test phase shift
        mtz.expand_to_p1()
        self.assertEqual(mtz.row_as_dict((7, 5, -1)), d1)
        mtz.ensure_asu()
        # (75-1) is not in P1 asu, we get Friedel mate
        self.assertEqual(mtz.row_as_dict((7, 5, -1)), {})
        d2 = mtz.row_as_dict((-7, -5, 1))
        self.assertAlmostEqual(d2['FC'], 22.7, delta=1e-5)
        self.assertAlmostEqual(d2['PHIC'], 360-91.1, delta=1e-5)
        self.assertAlmostEqual(d2['HLA'], d1['HLA'], delta=1e-5)
        self.assertAlmostEqual(d2['HLB'], -d1['HLB'], delta=1e-5)
        self.assertAlmostEqual(d2['HLC'], d1['HLC'], delta=1e-5)
        self.assertAlmostEqual(d2['HLD'], -d1['HLD'], delta=1e-5)
        # Expanding mtz file with SFTOOLS from CCP4 (command EXPAND)
        # and cctbx (miller_array.expand_to_p1().map_to_asu()) gives
        # different phases (and H-L coefficients) for some reflections,
        # including (-12, 5, 1).
        # Here we agree with cctbx (which gives PHIC -28.9).
        d3 = mtz.row_as_dict((-12, 5, 1))
        self.assertAlmostEqual(d3['FC'], 22.7, delta=1e-5)
        self.assertAlmostEqual(d3['PHIC'], 360-28.9, delta=1e-5)
        self.assertAlmostEqual(d3['HLA'], -5.28884, delta=1e-5)
        self.assertAlmostEqual(d3['HLB'], 10.0395, delta=1e-4)
        self.assertAlmostEqual(d3['HLC'], 1.53099, delta=1e-5)
        self.assertAlmostEqual(d3['HLD'], 4.64824, delta=1e-5)

if __name__ == '__main__':
    unittest.main()
