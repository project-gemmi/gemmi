#!/usr/bin/env python

import os
import pathlib
import tempfile
import unittest

import gemmi
from common import full_path


class PathLike:
    """A path-like object that is not pathlib.Path."""
    def __init__(self, path):
        self.path = path
    def __fspath__(self):
        return self.path


class TestPathLike(unittest.TestCase):
    """File names can be str or os.PathLike, such as pathlib.Path (#101)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = pathlib.Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def variants(self, name):
        path = full_path(name)
        return [path, pathlib.Path(path), PathLike(path)]

    def test_cif_read(self):
        n_blocks = len(gemmi.cif.read(full_path('misc.cif')))
        for path in self.variants('misc.cif'):
            self.assertEqual(len(gemmi.cif.read(path)), n_blocks)
            self.assertEqual(len(gemmi.cif.read_file(path)), n_blocks)
            doc = gemmi.cif.Document()
            doc.parse_file(path)
            self.assertEqual(len(doc), n_blocks)
        for path in self.variants('1pfe.json'):
            self.assertEqual(gemmi.cif.read_mmjson(path).sole_block().name, '1PFE')

    def test_cif_write(self):
        doc = gemmi.cif.read_string('data_a\n_item 1\n')
        block = doc[0]
        for i, wrap in enumerate([str, pathlib.Path, PathLike]):
            out = str(self.dir / f'doc{i}.cif')
            doc.write_file(wrap(out))
            self.assertEqual(gemmi.cif.read(out)[0].find_value('_item'), '1')
            out = str(self.dir / f'block{i}.cif')
            block.write_file(wrap(out))
            self.assertEqual(gemmi.cif.read(out)[0].find_value('_item'), '1')

    def test_read_structure(self):
        for path in self.variants('5e5z.pdb'):
            self.assertEqual(len(gemmi.read_structure(path)), 1)
            self.assertEqual(len(gemmi.read_pdb(path)), 1)
        for path in self.variants('5i55.cif'):
            self.assertEqual(gemmi.read_structure(path).info['_entry.id'], '5I55')

    def test_write_pdb(self):
        st = gemmi.read_structure(full_path('5e5z.pdb'))
        for i, wrap in enumerate([str, pathlib.Path, PathLike]):
            out = str(self.dir / f'out{i}.pdb')
            st.write_pdb(wrap(out), gemmi.PdbWriteOptions(minimal=True))
            self.assertEqual(len(gemmi.read_structure(out)), 1)
            out = str(self.dir / f'min{i}.pdb')
            st.write_minimal_pdb(wrap(out))
            self.assertEqual(len(gemmi.read_structure(out)), 1)
            out = str(self.dir / f'kwargs{i}.pdb')
            st.write_pdb(wrap(out), minimal=True)  # deprecated form
            self.assertEqual(len(gemmi.read_structure(out)), 1)

    def test_small_structure_and_gz_size(self):
        for path in self.variants('2013551.cif'):
            self.assertEqual(len(gemmi.read_small_structure(path).sites),
                             len(gemmi.read_small_structure(full_path('2013551.cif')).sites))
        for path in self.variants('1lzh.pdb.gz'):
            self.assertEqual(gemmi.estimate_uncompressed_size(path),
                             gemmi.estimate_uncompressed_size(full_path('1lzh.pdb.gz')))

    def test_ccp4(self):
        shape = gemmi.read_ccp4_map(full_path('5i55_tiny.ccp4')).grid.shape
        for path in self.variants('5i55_tiny.ccp4'):
            self.assertEqual(gemmi.read_ccp4_map(path).grid.shape, shape)
            self.assertEqual(gemmi.read_ccp4_header(path).header_i32(1), shape[0])
        mask_shape = gemmi.read_ccp4_mask(full_path('1pfe_asu.msk.gz')).grid.shape
        for path in self.variants('1pfe_asu.msk.gz'):
            self.assertEqual(gemmi.read_ccp4_mask(path).grid.shape, mask_shape)
        ccp4 = gemmi.read_ccp4_map(full_path('5i55_tiny.ccp4'))
        for i, wrap in enumerate([str, pathlib.Path, PathLike]):
            out = str(self.dir / f'map{i}.ccp4')
            ccp4.write_ccp4_map(wrap(out))
            self.assertEqual(gemmi.read_ccp4_map(out).grid.shape, shape)

    def test_mtz(self):
        for path in self.variants('5e5z.mtz'):
            mtz = gemmi.read_mtz_file(path)
            self.assertEqual(mtz.nreflections, gemmi.read_mtz_file(full_path('5e5z.mtz')).nreflections)
        for i, wrap in enumerate([str, pathlib.Path, PathLike]):
            out = str(self.dir / f'out{i}.mtz')
            mtz.write_to_file(wrap(out))
            self.assertEqual(gemmi.read_mtz_file(out).nreflections, mtz.nreflections)

    def test_other_readers(self):
        for path in self.variants('ener_lib.cif'):
            enerlib = gemmi.EnerLib()
            enerlib.read(path)
        for path in self.variants('INTEGRATE-tiny.HKL'):
            self.assertGreater(gemmi.read_xds_ascii(path).data_size, 0)

    def test_rejected_types(self):
        path = full_path('misc.cif')
        for bad in [os.fsencode(path), 5, None, PathLike(os.fsencode(path))]:
            with self.assertRaises(TypeError):
                gemmi.cif.read_file(bad)
        # bytes were not accepted before either; a missing file is not a TypeError
        with self.assertRaises(Exception) as ctx:
            gemmi.cif.read_file(pathlib.Path(path + '.missing'))
        self.assertNotIsInstance(ctx.exception, TypeError)

    def test_unicode_path(self):
        path = self.dir / 'päth α.cif'
        path.write_text('data_a\n_item 1\n', encoding='utf-8')
        self.assertEqual(gemmi.cif.read_file(path)[0].find_value('_item'), '1')
        self.assertEqual(gemmi.cif.read_file(str(path))[0].find_value('_item'), '1')


if __name__ == '__main__':
    unittest.main()
