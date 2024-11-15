#!/usr/bin/env python

import sys
import gemmi
import numpy as np
import statsmodels.api as sm

# Usage example: qq.py path/to/df.map 1.5  >plot.dat
assert len(sys.argv) == 3
df_map_path = sys.argv[1]
d_min = float(sys.argv[2])

# read map
df = gemmi.read_ccp4_map(df_map_path, setup=True)
grid = df.grid

# get 1D array of grid points corresponding to a single ASU
map_values = grid.array[grid.masked_asu().mask_array == 0]

# resample data to get as many points as we'd have at the Shannon limit
current_voxel_vol = grid.unit_cell.volume / grid.point_count
shannon_voxel_vol = (d_min / 2)**3
npoints = int(round(current_voxel_vol / shannon_voxel_vol * len(map_values)))

# get quantiles
map_values.sort()
quantile_levels = np.linspace(0, 1, npoints)
data = np.quantile(map_values, quantile_levels)

# print data for the plot
pplot = sm.ProbPlot(data, fit=True)
for x in zip(pplot.theoretical_quantiles,
             pplot.sample_quantiles - pplot.theoretical_quantiles):
    print('%g %g' % x)
