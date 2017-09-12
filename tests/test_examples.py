#!/usr/bin/env python

import os
import sys
import unittest

import gemmi

TOP_DIR = os.path.join(os.path.dirname(__file__), "..")
CIF_FILE = os.path.join(TOP_DIR, "tests", "5i55.cif")
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
    import aafreq
  def test_col_order(self):
    import col_order
  def test_hello(self):
    import hello
  def test_matthews(self):
    import matthews
  def test_monomers(self):
    import monomers  # nothing is run here, we only test import
  def test_simple_search(self):
    import simple_search
  def test_weight(self):
    import weight
    weight.main()

if __name__ == '__main__':
  unittest.main()

# vim:sw=2:ts=2:et
