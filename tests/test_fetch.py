#!/usr/bin/env python

import unittest
import urllib.error
import http.client

import gemmi
import gemmi.fetch
from common import full_path

class TestFetch(unittest.TestCase):
    def test_fetch(self):
        cell = gemmi.read_structure(full_path('5i55.cif')).cell
        for site in "REJ":
            for use_cif in [True, False]:
                try:
                    st = gemmi.fetch.read_structure_with_code('5i55', pdb_site=site,
                                                              use_cif=use_cif)
                    self.assertEqual(cell, st.cell)
                except urllib.error.URLError as e:
                    self.skipTest(f"Cannot connect to {site} ({e})")
                except http.client.RemoteDisconnected as e:
                    self.skipTest(f"Couldn't fetch file form {site}: ({e})")
