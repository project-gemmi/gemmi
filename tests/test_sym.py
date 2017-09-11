#!/usr/bin/env python

import unittest
from gemmi import sym

CANONICAL_SINGLES = {
 # all crystallographic possibilities
 "x"     : [ 1, 0, 0,  0],
 "y"     : [ 0, 1, 0,  0],
 "z"     : [ 0, 0, 1,  0],
 "-x"    : [-1, 0, 0,  0],
 "-y"    : [ 0,-1, 0,  0],
 "-z"    : [ 0, 0,-1,  0],
 "x-y"   : [ 1,-1, 0,  0],
 "-x+y"  : [-1, 1, 0,  0],
 "x+1/2" : [ 1, 0, 0,  6],
 "x+1/4" : [ 1, 0, 0,  3],
 "x+3/4" : [ 1, 0, 0,  9],
 "y+1/2" : [ 0, 1, 0,  6],
 "y+1/4" : [ 0, 1, 0,  3],
 "y+3/4" : [ 0, 1, 0,  9],
 "z+1/2" : [ 0, 0, 1,  6],
 "z+1/3" : [ 0, 0, 1,  4],
 "z+1/4" : [ 0, 0, 1,  3],
 "z+1/6" : [ 0, 0, 1,  2],
 "z+2/3" : [ 0, 0, 1,  8],
 "z+3/4" : [ 0, 0, 1,  9],
 "z+5/6" : [ 0, 0, 1, 10],
 "-x+1/2": [-1, 0, 0,  6],
 "-x+1/4": [-1, 0, 0,  3],
 "-x+3/4": [-1, 0, 0,  9],
 "-y+1/2": [ 0,-1, 0,  6],
 "-y+1/4": [ 0,-1, 0,  3],
 "-y+3/4": [ 0,-1, 0,  9],
 "-z+1/2": [ 0, 0,-1,  6],
 "-z+1/3": [ 0, 0,-1,  4],
 "-z+1/4": [ 0, 0,-1,  3],
 "-z+1/6": [ 0, 0,-1,  2],
 "-z+2/3": [ 0, 0,-1,  8],
 "-z+3/4": [ 0, 0,-1,  9],
 "-z+5/6": [ 0, 0,-1, 10],
}

OTHER_SINGLES = {
 # order and letter case may vary
 "Y-x"   : [-1, 1, 0,  0],
 "-X"    : [-1, 0, 0,  0],
 "-1/2+Y": [ 0, 1, 0,  -6],
 # we want to handle non-crystallographic translations
 "x+3"   : [ 1, 0, 0,  3*12],
 "1+Y"   : [ 0, 1, 0,  1*12],
 "-2+Y"  : [ 0, 1, 0, -2*12],
 "-z-5/6": [ 0, 0,-1, -10],
}

class TestSymmetry(unittest.TestCase):
    def test_parse_triplet_part(self):
        for single, row in CANONICAL_SINGLES.items():
            calculated = sym.parse_triplet_part(single)
            self.assertEqual(calculated, row)
        for single, row in OTHER_SINGLES.items():
            calculated = sym.parse_triplet_part(single)
            self.assertEqual(calculated, row)

    def test_make_triplet_part(self):
        for single, row in CANONICAL_SINGLES.items():
            calculated = sym.make_triplet_part(*row)
            self.assertEqual(calculated, single)


    def test_parse_triplet(self):
        # randomly combine SINGLES
        pass

if __name__ == '__main__':
    unittest.main()
