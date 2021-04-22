#!/usr/bin/env python3

import glob

SEPARATOR_LINES = ('//\n', '\n')
for header in sorted(glob.glob('include/gemmi/*.hpp')):
    without_include = header.partition('/')[2]
    print('\n' + without_include)
    with open(header) as f:
        line = f.readline()
        if line.startswith('// Copyright'):
            line = f.readline()
            assert line in SEPARATOR_LINES, line
            line = f.readline()
        while line not in SEPARATOR_LINES:
            assert line[2] == ' ', line
            print('    ' + line[3:-1])
            line = f.readline()
