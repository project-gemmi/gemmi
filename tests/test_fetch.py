#!/usr/bin/env python

import unittest
import urllib.error
import http.client

import gemmi
import gemmi.fetch
from common import full_path

class TestFetch(unittest.TestCase):
    def fetch(self, site):
        cell = gemmi.read_structure(full_path('5i55.cif')).cell
        for use_cif in [True, False]:
            try:
                st = gemmi.fetch.read_structure_with_code('5i55', pdb_site=site,
                                                          use_cif=use_cif)
                self.assertEqual(cell, st.cell)
            except urllib.error.URLError as e:
                self.skipTest(f"Cannot connect to {site} ({e})")
            except http.client.RemoteDisconnected as e:
                self.skipTest(f"Couldn't fetch file form {site}: ({e})")
    @unittest.skipIf(__name__ != '__main__', "Skip during unittest discover")
    def test_rcsb(self):
        self.fetch('R')
    @unittest.skipIf(__name__ != '__main__', "Skip during unittest discover")
    def test_pdbe(self):
        self.fetch('E')
    @unittest.skipIf(__name__ != '__main__', "Skip during unittest discover")
    def test_pdbj(self):
        self.fetch('J')

if __name__ == '__main__':
    unittest.main()
