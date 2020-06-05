#!/usr/bin/env python3
# IF YOU ADD OR REMOVE LINES, ADJUST :lines: in docs/hkl.rst

import gemmi
import numpy
import pandas
from matplotlib import pyplot

MTZ_PATH = 'tests/5wkd_phases.mtz.gz'
SFCIF_PATH = 'tests/r5wkdsf.ent'

# make DataFrame from MTZ file
mtz = gemmi.read_mtz_file(MTZ_PATH)
mtz_data = numpy.array(mtz, copy=False)
mtz_df = pandas.DataFrame(data=mtz_data, columns=mtz.column_labels())
# (optional) store Miller indices as integers
mtz_df = mtz_df.astype({label: 'int32' for label in 'HKL'})

# make DataFrame from mmCIF file
cif_doc = gemmi.cif.read(SFCIF_PATH)
rblock = gemmi.as_refln_blocks(cif_doc)[0]
cif_df = pandas.DataFrame(data=rblock.make_miller_array(),
                          columns=['H','K','L'])
cif_df['F_meas_au'] = rblock.make_float_array('F_meas_au')
cif_df['d'] = rblock.make_d_array()

# merge DataFrames
df = pandas.merge(mtz_df, cif_df, on=['H', 'K', 'L'])


# plot FP from MTZ as a function of F_meas_au from mmCIF
pyplot.rc('font', size=8)
pyplot.figure(figsize=(2, 2))
pyplot.scatter(x=df['F_meas_au'], y=df['FP'],
               marker=',', s=1, linewidths=0)
pyplot.xlim(xmin=0)
pyplot.ylim(ymin=0)
pyplot.show()


# plot the ratio FP : F_meas_au as a function of 1/d
pyplot.figure(figsize=(3, 3))
pyplot.scatter(x=1/df['d'],
               y=df['FP']/df['F_meas_au'],
               c=df['K'],  # color by index k
               marker='.', s=16, linewidths=0)
pyplot.xlim(xmin=0)
pyplot.show()
