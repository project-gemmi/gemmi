#!/usr/bin/env python
# Check the syminfo.lib file. Usage: syminfo.py /path/to/syminfo.lib

from __future__ import print_function
import sys
from gemmi import sym

def parse_syminfo(path):
    data = []
    cur = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line[0] == '#':
                continue
            #print('"%s"' % line)
            if line == 'begin_spacegroup':
                assert cur is None, line
                cur = {'symops': [], 'cenops': []}
                continue
            assert cur is not None, line
            if line == 'end_spacegroup':
                for must_have in ['basisop', 'ccp4', 'number', 'hall', 'xhm']:
                    assert must_have in cur, must_have
                for must_have_list in ['symops', 'cenops']:
                    assert len(cur[must_have_list]) > 0
                data.append(cur)
                cur = None
            elif line.startswith('number '):
                cur['number'] = int(line[6:])
            elif line.startswith('basisop '):
                cur['basisop'] = line[8:]
            elif line.startswith('symbol ccp4 '):
                cur['ccp4'] = int(line[12:])
            elif line.startswith('symbol xHM '):
                cur['xhm'] = line[11:].strip(" '")
            elif line.startswith('symbol Hall '):
                cur['hall'] = line[12:].strip(" '")
            elif line.startswith('symop '):
                cur['symops'].append(sym.Op(line[6:]))
            elif line.startswith('cenop '):
                cur['cenops'].append(sym.Op(line[6:]))
    return data

def main():
    syminfo = parse_syminfo(sys.argv[1])
    print(len(syminfo))
    for entry in syminfo:
        ccp4 = entry['ccp4']
        hall = entry['hall']
        assert ccp4 == 0 or ccp4 % 1000 == entry['number']
        try:
            hall_ops = sym.symops_from_hall(hall)
        except RuntimeError:  # a few entries with denominator 8 in basisop
            continue
        assert len(hall_ops.cen_ops) == len(entry['cenops'])
        generated = set(hall_ops)
        given = set(s * c for s in entry['symops'] for c in entry['cenops'])
        #assert len(generated) == len(given)
        if given != generated:
            print(hall)
            print('common:',  '  '.join(x.triplet() for x in given & generated))
            print('given:    ', '  '.join(x.triplet() for x in given - generated))
            print('generated:', '  '.join(x.triplet() for x in generated - given))
            print()


main()
