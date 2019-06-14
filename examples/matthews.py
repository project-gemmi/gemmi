#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Usage:
#  ./matthews.py $PDB_DIR/structures/divided/mmCIF | tee data.tsv
#  ./matthews.py plot data.tsv
# Plotting uses numpy, statsmodels, pandas, matplotlib, seaborn.
# There is also a "check" command that was used to estimate data quality.

from __future__ import print_function
import datetime
import sys
import csv
from gemmi import cif
import util


def parse_date(date_str):
    year, month, day = date_str.split('-')
    return datetime.date(int(year), int(month), int(day))


def gather_data():
    "read mmCIF files and write down a few numbers (one file -> one line)"
    writer = csv.writer(sys.stdout, dialect='excel-tab')
    writer.writerow(['code', 'na_chains', 'vs', 'vm', 'd_min', 'date', 'group'])
    for path in util.get_file_paths_from_args():
        block = cif.read(path).sole_block()
        code = cif.as_string(block.find_value('_entry.id'))
        na = sum('nucleotide' in t[0]
                 for t in block.find_values('_entity_poly.type'))
        vs = block.find_value('_exptl_crystal.density_percent_sol')
        vm = block.find_value('_exptl_crystal.density_Matthews')
        d_min = block.find_value('_refine.ls_d_res_high')
        dep_date_tag = '_pdbx_database_status.recvd_initial_deposition_date'
        dep_date = parse_date(block.find_values(dep_date_tag).str(0))
        group = block.find_value('_pdbx_deposit_group.group_id')
        writer.writerow([code, na, vs, vm, d_min, dep_date, group])


def plot(our_csv):
    from numpy import array
    import seaborn as sns
    x, y = [], []
    with open(our_csv) as csvfile:
        for row in csv.DictReader(csvfile, dialect='excel-tab'):
            try:
                vs = float(row['vs'])
                vm = float(row['vm'])
                d_min = float(row['d_min'])
            except ValueError:
                continue
            # skip entries with inconsistent Vs and Vm
            if vm == 0 or abs(vs - 100 * (1 - 1.23 / vm)) > 1.0:
                continue
            if row['date'] >= '2015-01-01':  # and not row['group']:
                x.append(d_min)
                y.append(vs)
    print('Plotting kernel density estimation from', len(x), 'points.')
    sns.set_style("whitegrid")
    sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 1.5})
    g = sns.JointGrid(array(x), array(y), space=0, xlim=(0, 4.5), ylim=(20, 90))
    # calculation time is proportional to gridsize^2; gridsize=100 is default
    g = g.plot_joint(sns.kdeplot, n_levels=30, bw=(0.05, 0.3), gridsize=300,
                     shade=True, shade_lowest=False, cmap="gnuplot2_r")
    sns.plt.sca(g.ax_marg_x)
    sns.kdeplot(g.x, shade=True, bw=0.03, gridsize=10000, vertical=False)
    sns.plt.sca(g.ax_marg_y)
    sns.kdeplot(g.y, shade=True, bw=0.2, gridsize=5000, vertical=True)
    g.ax_joint.set(xlabel=u'$d_{min}$ [Å]', ylabel='solvent volume [%]')
    #g.annotate(lambda a, b: 0, template="1970s - 2014")
    g.annotate(lambda a, b: 0, template="2015 - 4/2017")
    sns.plt.show()


# This function compares our data with the CSV file available from B.Rupp's
# website. It also does some other checks on both our and BR's CSV.
# Spoilers:
# The reso column in pdb_02_06_2013_pro_sorted_flagged_highest_cs.csv
# is differs from _refine.ls_d_res_high in about 50 cases.
# RB's data has entries with negative volume of solvent.
# PDB data has about 5000 entries with inconsistent Vm and Vs.
def check_with_rupps_data(our_csv, rupps_csv):
    with open(rupps_csv) as csvfile:
        rb_data = {row['code']: row for row in csv.DictReader(csvfile)}
    # check for suspiciously low Vs
    for r in rb_data.values():
        if float(r['vs']) < 5:
            print('Suspicious RB data: %s has Vs: %s' % (r['code'], r['vs']))

    with open(our_csv) as csvfile:
        for row in csv.DictReader(csvfile, dialect='excel-tab'):
            rb = rb_data.get(row['code'])
            try:
                vs = float(row['vs'])
                vm = float(row['vm'])
                if vm == 0:
                    print('Bad PDB data: Vm=0 in', row['code'])
                vs_calc = 100 * (1 - 1.23 / vm)
                if abs(vs - vs_calc) > 1.0:
                    print('Inconsistent PDB Vs/Vm: %s %.1f (%s->)%.1f Δ=%+.1f'
                          % (row['code'], vs, row['vm'], vs_calc, vs_calc-vs),
                          end='')
                    print('  (RB: %s->%s)' % (rb['vm'], rb['vs']) if rb else '')
            except (ValueError, ZeroDivisionError):
                pass
            if rb:
                try:  # d_min can be absent in pdb data (NMR entries)
                    d_min = float(row['d_min'])
                except ValueError:
                    continue
                diff = d_min - float(rb['reso'])
                if abs(diff) > 0.01:
                    print('PDB and RB differ: ', row['code'], row['d_min'],
                          rb['reso'], '%+.2f' % diff)


if sys.argv[1] == 'check':
    check_with_rupps_data(our_csv=sys.argv[2], rupps_csv=sys.argv[3])
elif sys.argv[1] == 'plot':
    plot(sys.argv[2])
else:
    gather_data()
