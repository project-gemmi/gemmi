#!/usr/bin/env python

import os
import sys
import unittest

import gemmi

TOP_DIR = os.path.join(os.path.dirname(__file__), "..")
CIF_FILE = os.path.join(TOP_DIR, "1YJP.cif")
EXAMPLE_DIR = os.path.join(TOP_DIR, "examples")
sys.path.append(EXAMPLE_DIR)


class TestExamples(unittest.TestCase):
  def test_aafreq(self):
    sys.argv = ['aafreq.py', CIF_FILE]
    import aafreq
  def test_col_order(self):
    pass
  def test_hello(self):
    sys.argv = ['hello.py', CIF_FILE]
    import hello
  def test_matthews(self):
    pass
  def test_monomers(self):
    pass
  def test_simple_search(self):
    pass
  def test_weight(self):
    pass

if __name__ == '__main__':
  unittest.main()

# vim:sw=2:ts=2:et
