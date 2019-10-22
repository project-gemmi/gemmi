#!/usr/bin/env python3

# Plots the data from rama_gather.py

import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def plot(data_file, label, output=None):
    x, y = [], []
    for line in open(data_file):
        phi, psi = line.split()
        if phi != 'nan' and psi != 'nan':
            x.append(float(phi))
            y.append(float(psi))
    print('Plotting %d points for %s' % (len(x), label))

    plt.figure(figsize=(5.5, 5.5))
    plt.title('%s, %d points.' % (label, len(x)), fontsize=12)
    plt.xlim([-180, 180])
    plt.ylim([-180, 180])
    ax = plt.gca()
    ax.xaxis.set_major_locator(MultipleLocator(60))
    ax.yaxis.set_major_locator(MultipleLocator(60))
    plt.xlabel(r'$\phi$', fontsize=14)
    plt.ylabel(r'$\psi$', fontsize=14, labelpad=0)
    plt.grid(color='#AAAAAA', linestyle='--')
    plt.hexbin(x, y, gridsize=2*180, bins='log', cmap='Blues')
    if output:
        plt.savefig(output, dpi=300)  # dpi=70 for small images in docs
    else:
        plt.show()

for aa in sys.argv[1:]:
    plot('ramas/%s.tsv' % aa, aa)
    #plot('ramas/%s.tsv' % aa, aa, 'ramas/%s.png' % aa)
