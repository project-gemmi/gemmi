#!/usr/bin/env python
# Read SF-mmCIF file and plot I/sigma as a function of 1/d^2.

import sys
from matplotlib import pyplot
import gemmi

for path in sys.argv[1:]:
    doc = gemmi.cif.read(path)
    rblock = gemmi.as_refln_blocks(doc)[0]
    intensity = rblock.make_float_array('intensity_meas')
    sigma = rblock.make_float_array('intensity_sigma')
    x = rblock.make_1_d2_array()
    y = intensity / sigma
    pyplot.figure(figsize=(5,5))
    pyplot.hexbin(x, y, gridsize=70, bins='log', cmap='gnuplot2')
    pyplot.show()
