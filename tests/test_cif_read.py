#!/usr/bin/env python

import gc
import os
import unittest
from gemmi import cif

class TestBlock(unittest.TestCase):
    def test_find(self):
        block = cif.read_string("""
            data_test
            _nonloop_a alpha
            _nonloop_b beta
            loop_ _la _lb A B C D
        """).sole_block()
        gc.collect()
        rows = list(block.find('_la'))
        self.assertEqual(list(rows[0]), ['A'])

        rows = list(block.find(['_lb', '_la']))  # changed order
        gc.collect()
        self.assertEqual(list(rows[0]), ['B', 'A'])
        self.assertEqual(rows[1][1], 'C')

        rows = list(block.find('_nonloop_', ['a', 'b']))
        gc.collect()
        self.assertEqual([list(r) for r in rows], [['alpha', 'beta']])

    def test_find_values(self):
        v1 = cif.read_string('data_a _v1 one')[0].find_values('_v1')
        gc.collect()
        self.assertListEqual(list(v1), ['one'])
        v2 = cif.read_string('data_a loop_ _v2 a b')[0].find_values('_v2')
        gc.collect()
        self.assertListEqual(list(v2), ['a', 'b'])

    def test_reading_gzipped_file(self):
        path = os.path.join(os.path.dirname(__file__), '1pfe.cif.gz')
        cif_doc = cif.read(path)
        categories = cif_doc.sole_block().get_mmcif_category_names()
        self.assertEqual(categories[0], '_entry.')
        self.assertEqual(len(categories), 72)

if __name__ == '__main__':
    unittest.main()
