#!/usr/bin/env python

import contextlib
import io
import os
import sys
import unittest
try:
    import numpy
except ImportError:
    numpy = None

TOP_DIR = os.path.join(os.path.dirname(__file__), "..")
PDB_FILE = os.path.join(TOP_DIR, "tests", "1orc.pdb")
CIF_FILE = os.path.join(TOP_DIR, "tests", "5i55.cif")
JSON_FILE = os.path.join(TOP_DIR, "tests", "1pfe.json")
MAP_FILE = os.path.join(TOP_DIR, "tests", "5i55_tiny.ccp4")
PDB_1PFE_FILE = os.path.join(TOP_DIR, "tests", "1pfe.cif.gz")
MASK_1PFE_FILE = os.path.join(TOP_DIR, "tests", "1pfe_asu.msk.gz")
EXAMPLE_DIR = os.path.join(TOP_DIR, "examples")
sys.path.insert(0, EXAMPLE_DIR)

class TestExamples(unittest.TestCase):
    def setUp(self):
        sys.argv = ['example', CIF_FILE]
        sys.stdout = open(os.devnull, 'w')
    def tearDown(self):
        sys.stdout.close()
        sys.stdout = sys.__stdout__
    def test_aafreq(self):
        import aafreq  # noqa: F401 - imported but unused
    def test_col_order(self):
        import col_order  # noqa: F401
    def test_hello(self):
        import hello  # noqa: F401
    def test_matthews(self):
        import matthews  # noqa: F401
    def test_monomers(self):
        import monomers  # noqa: F401
    def test_simple_search(self):
        import simple_search  # noqa: F401
    def test_weight(self):
        import weight
        weight.main()
    def test_long_geom(self):
        import long_geom
        self.assertEqual(long_geom.run(PDB_FILE), 0)
        self.assertEqual(long_geom.run(CIF_FILE), 1)

class TestExamples2(unittest.TestCase):
    def test_from_json(self):
        sys.argv = ['from_json.py', JSON_FILE, os.devnull]
        import from_json  # noqa: F401
    def test_map2mtz(self):
        sys.argv = ['map2mtz.py', MAP_FILE, os.devnull]
        import map2mtz  # noqa: F401
    @unittest.skipIf(sys.version_info[0] == 2, 'this example is Py3 only')
    def test_multiproc(self):
        sys.stdout = open(os.devnull, 'w')
        import multiproc  # noqa: F401
        multiproc.main(PDB_FILE)

    # In this example we use file produced with two commands:
    # $ cctbx.python -m mmtbx.command_line.mask 1pfe.pdb
    # $ mapmask mapin mask.ccp4 mskout 1pfe_asu.msk << eof
    # XYZLIM ASU
    # MASK CUT 0.5
    # AXIS X Z Y
    # MODE mapin
    # eof
    @unittest.skipIf(numpy is None or sys.version_info[0] == 2,
                     'this example is Py3 only, requires NumPy')
    def test_maskcheck(self):
        import maskcheck
        with io.StringIO() as buf, contextlib.redirect_stdout(buf):
            maskcheck.maskcheck(MASK_1PFE_FILE, PDB_1PFE_FILE)
            output = buf.getvalue()
        expected_output = '''\
Generating CCTBX-compatible mask ...
Size: 72 x 72 x 150, total 777600 points
File-Gemmi Count Fraction
0-0       509676  65.54%
1-1       265788  34.18%
0-1         1452   0.19%
1-0          684   0.09%
'''
        self.assertEqual(output.splitlines(), expected_output.splitlines())

if __name__ == '__main__':
    unittest.main()
