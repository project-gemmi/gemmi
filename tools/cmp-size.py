#!/usr/bin/env python3
# Script for storing and comparing binary sizes.
# Uses $PWD/.size-save.txt to store output of the size command.

import os
import sys
import subprocess

SAVE = '.size-save.txt'
paths = sys.argv[1:]

old_output = None
if os.path.exists(SAVE):
    with open(SAVE) as f:
        old_output = f.read()

new_output = subprocess.check_output(['size'] + paths)
with open(SAVE, "wb") as f:
    f.write(new_output)

if old_output is None:
    sys.exit()

def parse_size(output):
    lines = output.splitlines()[1:]
    values = []
    for line in lines:
        text, data, bss, _ = line.split(None, 3)
        values.append([int(text), int(data), int(bss)])
    return values

for old, new, p in zip(parse_size(old_output), parse_size(new_output), paths):
    print("size of %s --> " % os.path.basename(p), end="")
    if old == new:
        print("the same")
    else:
        changes = tuple(a - b for a, b in zip(new, old))
        print("text: %+d  data: %+d  bss: %+d" % changes)
