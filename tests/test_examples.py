#!/usr/bin/env python

import os
import sys
import unittest

TOP_DIR = os.path.join(os.path.dirname(__file__), "..")
CIF_FILE = os.path.join(TOP_DIR, "tests", "5i55.cif")
JSON_FILE = os.path.join(TOP_DIR, "tests", "1pfe.json")
EXAMPLE_DIR = os.path.join(TOP_DIR, "examples")
sys.path.append(EXAMPLE_DIR)

class TestExamples(unittest.TestCase):
    def setUp(self):
        sys.argv = ['example', CIF_FILE]
        sys.stdout = open(os.devnull, 'w')
    def tearDown(self):
        sys.stdout.close()
        sys.stdout = sys.__stdout__
    def test_aafreq(self):
        import aafreq  # noqa: E401 - imported but unused
    def test_col_order(self):
        import col_order  # noqa: E401
    def test_hello(self):
        import hello  # noqa: E401
    def test_matthews(self):
        import matthews  # noqa: E401
    def test_monomers(self):
        import monomers  # noqa: E401
    def test_simple_search(self):
        import simple_search  # noqa: E401
    def test_weight(self):
        import weight
        weight.main()

class TestFromJson(unittest.TestCase):
    sys.argv = ['from_json.py', JSON_FILE, os.devnull]
    import from_json

if __name__ == '__main__':
    unittest.main()
