#!/usr/bin/env python

import os
import json
import unittest

import gemmi

def nonempty_lines(s):
  return [line for line in s.splitlines() if line and line[0] != '#']

class TestCifAsJson(unittest.TestCase):
  def setUp(self):
    self.basename = os.path.join(os.path.dirname(__file__), "misc")
  def test_misc(self):
    cif_doc = gemmi.cif.read_file(self.basename + ".cif")
    json_str = cif_doc.as_json()
    json_from_cif = json.loads(json_str)
    with open(self.basename + ".json") as f:
      reference_json = json.load(f)
    self.assertEqual(json_from_cif, reference_json)

  def test_cif_as_string(self):
    with open(self.basename + ".cif", 'rb') as f:
      cif_orig = f.read().decode('utf-8')
    cif_doc = gemmi.cif.read_string(cif_orig)
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

if __name__ == '__main__':
  unittest.main()

# vim:sw=2:ts=2:et
