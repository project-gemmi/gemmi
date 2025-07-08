#!/usr/bin/env python

import unittest

import gemmi
import gemmi.fetch
from common import full_path

class TestFetch(unittest.TestCase):
    def test_fetch(self):
        cell = gemmi.read_structure(full_path('5i55.cif')).cell
        for site in "REJ":
            for use_cif in [True, False]:
                st = gemmi.fetch.read_structure_with_code('5i55', pdb_site=site,
                                                          use_cif=use_cif)
                self.assertEqual(cell, st.cell)
