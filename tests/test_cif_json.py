#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import io
import json
import unittest

from gemmi import cif

def nonempty_lines(s):
    return [line for line in s.splitlines() if line and line[0] != '#']

class TestCifAsJson(unittest.TestCase):
    def setUp(self):
        self.basename = os.path.join(os.path.dirname(__file__), 'misc')
    def test_misc(self):
        cif_doc = cif.read_file(self.basename + '.cif')
        json_str = cif_doc.as_json()
        json_from_cif = json.loads(json_str)
        with io.open(self.basename + '.json', encoding='utf-8') as f:
            reference_json = json.load(f)
        self.assertEqual(json_from_cif, reference_json)

    def test_cif_as_string(self):
        with open(self.basename + '.cif', 'rb') as f:
            cif_orig = f.read().decode('utf-8').replace('\r\n', '\n')
        cif_doc = cif.read_string(cif_orig)
        formatting_changes = {
            '# comment\n': '',
            '_x _y _z': '_x\n_y\n_z',
            '_b 2 3 4': '_b\n2\n3\n4',
            'pdbx_details    ': 'pdbx_details ',
        }
        for k, v in formatting_changes.items():
            cif_orig = cif_orig.replace(k, v)
        cif_out = cif_doc.as_string()
        self.assertListEqual(nonempty_lines(cif_orig), nonempty_lines(cif_out))

class TestSpecialFilename(unittest.TestCase):
    def test_open(self):
        if os.name == 'nt': return
        path = os.path.join(os.path.dirname(__file__), u'gęś¼.cif')
        doc = cif.read_file(path)
        self.assertEqual(doc[0].find_value('_tricky'), 'filename')

if __name__ == '__main__':
    unittest.main()
