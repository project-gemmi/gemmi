#!/usr/bin/env python

import math
import sys
import unittest
import zlib
import gemmi
from common import full_path, get_path_for_tempfile, assert_numpy_equal, numpy

class TestFloatGrid(unittest.TestCase):
    def test_reading(self):
        m = gemmi.read_ccp4_map(full_path('5i55_tiny.ccp4'))
        self.assertEqual(m.grid.nu, 8)
        self.assertEqual(m.grid.nv, 6)
        self.assertEqual(m.grid.nw, 10)
        self.assertEqual(m.header_i32(28), 0)
        m.set_header_i32(28, 20140)  # set NVERSION
        self.assertEqual(m.header_i32(28), 20140)
        dmax = m.header_float(21)
        self.assertEqual(dmax, max(p.value for p in m.grid))
        self.assertNotEqual(m.grid.axis_order, gemmi.AxisOrder.XYZ)
        m.setup(float('nan'))
        self.assertEqual(m.grid.axis_order, gemmi.AxisOrder.XYZ)
        self.assertEqual(m.grid.nu, 60)
        self.assertEqual(m.grid.nv, 24)
        self.assertEqual(m.grid.nw, 60)
        self.assertEqual(m.grid.point_count, 60 * 24 * 60)
        self.assertEqual(m.header_float(14), 90.0)  # 14 - alpha angle
        self.assertEqual(m.grid.unit_cell.alpha, 90.0)
        self.assertEqual(m.grid.spacegroup.ccp4, 4)  # P21

        extent = m.get_extent()
        self.assertEqual(extent.minimum.tolist(), [-1e-9]*3)
        nu = m.grid.nu
        self.assertTrue((nu - 1.) / nu < extent.maximum.x < 1)

        pos = gemmi.Position(19.4, 3., 21.)
        frac = m.grid.unit_cell.fractionalize(pos)
        pos_value = 2.1543798446655273
        self.assertAlmostEqual(m.grid.interpolate_value(pos), pos_value)
        self.assertAlmostEqual(m.grid.interpolate_value(frac), pos_value)

        # this spacegroup has symop -x, y+1/2, -z
        m.grid.set_value(60-3, 24//2+4, 60-5, 100)  # image of (3, 4, 5)
        self.assertEqual(m.grid.get_value(60-3, 24//2+4, 60-5), 100)
        self.assertTrue(math.isnan(m.grid.get_value(3, 4, 5)))
        m.grid.symmetrize_max()
        self.assertEqual(m.grid.get_value(3, 4, 5), 100)
        m.grid.set_value(3, 4, 5, float('nan'))
        self.assertTrue(math.isnan(m.grid.get_value(3, 4, 5)))
        m.grid.symmetrize_min()
        self.assertEqual(m.grid.get_value(3, 4, 5), 100)
        m.grid.set_value(60-3, 24//2+4, 60-5, 80)
        m.grid.symmetrize_avg()
        self.assertEqual(m.grid.get_value(3, 4, 5), 90)
        m.grid.set_value(60-3, 24//2+4, 60-5, float('nan'))
        m.grid.symmetrize_max()
        self.assertEqual(m.grid.get_value(60-3, 24//2+4, 60-5), 90)
        if numpy:
            arr = m.grid.array
            self.assertEqual(arr.shape, (60, 24, 60))
            self.assertEqual(arr[3][4][5], 90)
            grid2 = gemmi.FloatGrid(arr)
            self.assertTrue(numpy.allclose(m.grid, grid2, atol=0.0, rtol=0,
                                           equal_nan=True))

    def test_new(self):
        N = 24
        m = gemmi.FloatGrid(N, N, N)
        self.assertEqual(m.nu, N)
        self.assertEqual(m.nv, N)
        self.assertEqual(m.nw, N)
        m.set_value(1,2,3, 1.0)
        self.assertEqual(m.sum(), 1.0)
        m.spacegroup = gemmi.find_spacegroup_by_name('C2')
        self.assertEqual(m.spacegroup.number, 5)
        m.symmetrize_max()
        self.assertEqual(m.sum(), 4.0)
        m.get_point(0, 0, 0).value += 1
        self.assertEqual(m.sum(), 5.0)
        m.fill(2.0)
        m.spacegroup = gemmi.find_spacegroup_by_name('P 62 2 2')
        self.assertEqual(len(m.spacegroup.operations()), 12)
        m.set_value(1, 2, 3, 0.0)
        m.symmetrize_min()
        self.assertEqual(m.sum(), 2 * N * N * N - 2 * 12)

    def test_grid_size(self):
        # original cell from 4a0g, and a cell with a <-> b
        cell = gemmi.UnitCell(79.442, 80.066, 136.939, 99.96, 107.12, 97.25)
        cell2 = gemmi.UnitCell(80.066, 79.442, 136.939, 99.96, 107.12, 97.25)
        self.assertAlmostEqual(cell.calculate_d([4, -32, 9]), 2.501956204)
        dmin = 2.5
        self.assertEqual(cell.get_hkl_limits(dmin), [31, 32, 54])
        self.assertEqual(cell2.get_hkl_limits(dmin), [32, 31, 54])
        # In prepare_asu_data(), max_k == (grid.nv - 1) / 2,
        # so for k=31 we need nv>=63, for k=32 we need nv>=65.
        # But when choosing grid size, we take the same n for almost equal
        # lengths (a and b here). So both nu and nv are 72 instead of 64.
        grid = gemmi.FloatGrid()
        grid.spacegroup = gemmi.SpaceGroup('P 1')
        grid.unit_cell = cell
        grid.set_size_from_spacing(dmin / 2, gemmi.GridSizeRounding.Up)
        self.assertTrue([grid.nu, grid.nv, grid.nw] == [72, 72, 120])
        grid.unit_cell = cell2
        grid.set_size_from_spacing(dmin / 2, gemmi.GridSizeRounding.Up)
        self.assertTrue([grid.nu, grid.nv, grid.nw] == [72, 72, 120])

    def test_interpolation(self):
        # Setup grids for testing
        cell = gemmi.UnitCell(10.0,10.0,10.0, 90.0, 90.0, 90.0)
        moving_grid = gemmi.FloatGrid(10, 10, 10)
        moving_grid.spacegroup = gemmi.SpaceGroup('P 1')
        moving_grid.set_unit_cell(cell)
        moving_grid.set_value(0,0,0,1.0)
        moving_grid.set_value(5,5,5,1.0)

        interpolated_grid = gemmi.FloatGrid(10, 10, 10)
        interpolated_grid.spacegroup = gemmi.SpaceGroup('P 1')
        interpolated_grid.set_unit_cell(cell)

        # Test array interpolation
        positions = numpy.array([[0.0,0.0,0.0], [5.0,5.0,5.0]], 
                                dtype=numpy.float64)
        values = moving_grid.interpolate_position_array(positions)
        self.assertAlmostEqual(values[0], 1.0)
        self.assertAlmostEqual(values[1], 1.0)

        # Test map morphing
        # A simple single solid translation of the cell, limited to two points
        point_array_list = [
            numpy.array([[1,1,1], [6,6,6]],dtype=numpy.int32),]
        position_array_list = [
            numpy.array([[1.0,1.0,1.0], [6.0,6.0,6.0]],dtype=numpy.float32),]
        ts = [gemmi.Transform(),]
        com_m = [-1.0,-1.0,-1.0,]
        com_r = [0.0,0.0,0.0,]
        
        ts[0].vec.fromlist(
            (gemmi.Vec3(*com_m) - ts[0].mat.multiply(gemmi.Vec3(*com_r))
             ).tolist()
        )

        interpolated_grid.interpolate_grid_flexible(
            moving_grid, 
            point_array_list,
            position_array_list,
            ts,
        )
        self.assertAlmostEqual(interpolated_grid.get_value(1,1,1), 1.0)
        self.assertAlmostEqual(interpolated_grid.get_value(6,6,6), 1.0)
        
        ...


class TestCcp4Map(unittest.TestCase):
    @unittest.skipIf(numpy is None, "NumPy not installed.")
    def test_567_map(self):
        # make a small, contrived map
        data = numpy.arange(5*6*7, dtype=numpy.float32).reshape((5,6,7))
        cell = gemmi.UnitCell(150, 132, 140, 90, 90, 90)
        m = gemmi.Ccp4Map()
        m.grid = gemmi.FloatGrid(data, cell, gemmi.SpaceGroup('P 1'))
        m.update_ccp4_header()

        # write, read and compare
        tmp_path = get_path_for_tempfile(suffix='.ccp4')
        m.write_ccp4_map(tmp_path)
        with open(tmp_path, 'rb') as f:
            tmp_bytes = f.read()
        crc_mod_2_32 = zlib.crc32(tmp_bytes) % 4294967296
        if sys.byteorder == 'little':
            self.assertEqual(crc_mod_2_32, 4078044323)
        elif sys.byteorder == 'big':
            self.assertEqual(crc_mod_2_32, 372922578)
        self.assertEqual(tmp_bytes,
                         m.ccp4_header + m.grid.array.tobytes(order='A'))

        box = gemmi.FractionalBox()
        box.minimum = gemmi.Fractional(0.5/5, 1.5/6, 3.5/7)
        box.maximum = gemmi.Fractional(4.5/5, 2.5/6, 5.5/7)
        m.set_extent(box)
        cut_data = data[1:5, 2:3, 4:6]
        self.assertEqual(cut_data.shape, (4, 1, 2))
        self.assertTrue(numpy.array_equal(m.grid.array, cut_data))

        # CCP4 MAPMASK generated tests/iota_yzx.ccp4.gz from map m:
        #   mapmask mapin iota_full.ccp4 mapout iota_yzx.ccp4 << eof
        #   XYZLIM 0.3 0.9 3.4 3.45 -0.5 -0.3
        #   AXIS Y Z X
        #   MODE mapin
        #   eof
        yzx_path = full_path('iota_yzx.ccp4.gz')
        expanded_data = numpy.full(data.shape, float('nan'), dtype=data.dtype)
        expanded_data[1:5, 2:3, 4:6] = cut_data

        mcut = gemmi.read_ccp4_map(yzx_path, setup=False)
        self.assertEqual(mcut.axis_positions(), [2, 0, 1])
        mcut.setup(float('nan'), gemmi.MapSetup.ReorderOnly)
        self.assertEqual(mcut.axis_positions(), [0, 1, 2])
        self.assertFalse(mcut.full_cell())
        self.assertEqual(mcut.grid.axis_order, gemmi.AxisOrder.Unknown)
        self.assertTrue(numpy.array_equal(mcut.grid.array, cut_data))

        mcut = gemmi.read_ccp4_map(full_path(yzx_path), setup=False)
        self.assertFalse(mcut.full_cell())
        self.assertEqual(mcut.grid.axis_order, gemmi.AxisOrder.Unknown)
        mcut.setup(float('nan'), gemmi.MapSetup.NoSymmetry)
        self.assertTrue(mcut.full_cell())
        self.assertEqual(mcut.grid.axis_order, gemmi.AxisOrder.XYZ)
        assert_numpy_equal(self, mcut.grid.array, expanded_data)

        mcut = gemmi.read_ccp4_map(full_path(yzx_path), setup=True)
        self.assertTrue(mcut.full_cell())
        self.assertEqual(mcut.grid.axis_order, gemmi.AxisOrder.XYZ)
        assert_numpy_equal(self, mcut.grid.array, expanded_data)

    def test_normalize(self):
        yzx_path = full_path('iota_yzx.ccp4.gz')
        m = gemmi.read_ccp4_map(full_path(yzx_path), setup=True)
        m.grid.normalize()
        m.update_ccp4_header()
        mean = m.header_float(22)
        rms = m.header_float(55)
        self.assertAlmostEqual(mean, 0)
        self.assertAlmostEqual(rms, 1)

    @unittest.skipIf(numpy is None, "NumPy not installed.")
    def test_get_subarray(self):
        data = numpy.arange(5*6*7, dtype=numpy.float32).reshape((5,6,7))
        grid = gemmi.FloatGrid(data)
        sub = grid.get_subarray([4,-3,20], [5,10,4])
        self.assertEqual(sub[0][0][0], 195.)
        self.assertEqual(sub[1][2][3], 37.)
        self.assertEqual(sub[-1][-1][-1], 128.)
        grid.set_subarray(-sub, [4,-3,20])
        sub2 = grid.get_subarray([4,-3,20], [5,10,4])
        assert_numpy_equal(self, -sub, sub2)

    @unittest.skipIf(numpy is None, "NumPy not installed.")
    def test_setup_nosymmetry(self):
        m = gemmi.read_ccp4_map(full_path('5i55_tiny.ccp4'))
        orig_point_count = m.grid.point_count
        m.setup(0, gemmi.MapSetup.Full)
        full_nonzero = numpy.count_nonzero(m.grid.array)
        self.assertTrue(full_nonzero > orig_point_count)
        m = gemmi.read_ccp4_map(full_path('5i55_tiny.ccp4'))
        m.setup(0, gemmi.MapSetup.NoSymmetry)
        nosym_nonzero = numpy.count_nonzero(m.grid.array)
        self.assertEqual(full_nonzero, nosym_nonzero * 2)
        nonzero_ext = m.grid.get_nonzero_extent()
        span = nonzero_ext.maximum - nonzero_ext.minimum
        volume = span[0] * span[1] * span[2]
        self.assertAlmostEqual(orig_point_count / m.grid.point_count, volume)

if __name__ == '__main__':
    unittest.main()
