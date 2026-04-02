#!/usr/bin/env python3
# Script for storing and comparing binary sizes.
# Uses $PWD/.size-save.txt to store output keyed by file path.

import json
import os
import subprocess
import sys

SAVE = '.size-save.txt'


def load_sizes():
    if not os.path.exists(SAVE):
        return None
    with open(SAVE) as f:
        try:
            data = json.load(f)
        except ValueError:
            return None
    return data if isinstance(data, dict) else None


def measure_size(path):
    if path.endswith('.a'):
        print('Skipping static archive:', path, file=sys.stderr)
        return None
    result = subprocess.run(['size', path],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True)
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip()
        if message:
            print(message, file=sys.stderr)
        print('Skipping unsupported file:', path, file=sys.stderr)
        return None
    lines = result.stdout.splitlines()
    if len(lines) < 2:
        print('Unexpected size output for:', path, file=sys.stderr)
        return None
    try:
        text, data, bss, _ = lines[1].split(None, 3)
    except ValueError:
        print('Unexpected size output for:', path, file=sys.stderr)
        return None
    return [int(text), int(data), int(bss)]


old_sizes = load_sizes()
new_sizes = {}
for path in sys.argv[1:]:
    size_triplet = measure_size(path)
    if size_triplet is not None:
        new_sizes[path] = size_triplet

with open(SAVE, 'w') as f:
    json.dump(new_sizes, f, sort_keys=True)

if old_sizes is None:
    sys.exit()

for path in sys.argv[1:]:
    if path not in old_sizes or path not in new_sizes:
        continue
    old = old_sizes[path]
    new = new_sizes[path]
    print('size of %s --> ' % os.path.basename(path), end='')
    if old == new:
        print('the same')
    else:
        changes = tuple(a - b for a, b in zip(new, old))
        print('text: %+d  data: %+d  bss: %+d' % changes)
