#!/usr/bin/env python

import gc
import os
import unittest
from gemmi import cif

class TestBlock(unittest.TestCase):
    def test_find(self):
        block = cif.read_string("""
            data_test
            _one 1 _two 2 _three 3
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

    def test_set_pair(self):
        doc = cif.read_string('data_a _a 1 _b 2 _c 3')
        block = doc[0]
        block.set_pair('_d', '9')
        block.set_pair('_b', '8')
        self.assertEqual(tuple(block.find_pair('_a')), ('_a', '1'))
        self.assertEqual(block.find_value('_b'), '8')
        self.assertEqual(block.find_value('_c'), '3')
        self.assertEqual(block.find_value('_d'), '9')

    def test_reading_gzipped_file(self):
        path = os.path.join(os.path.dirname(__file__), '1pfe.cif.gz')
        cif_doc = cif.read(path)
        categories = cif_doc.sole_block().get_mmcif_category_names()
        self.assertEqual(categories[0], '_entry.')
        self.assertEqual(len(categories), 72)

if __name__ == '__main__':
    unittest.main()
