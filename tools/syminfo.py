#!/usr/bin/env python
# Check the syminfo.lib file. Usage: syminfo.py /path/to/syminfo.lib

from __future__ import print_function
import sys
import gemmi

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
                cur['xhm'] = line[11:].strip(" '").replace(' :', ':')
            elif line.startswith('symbol Hall '):
                cur['hall'] = line[12:].strip(" '")
            elif line.startswith('symop '):
                cur['symops'].append(gemmi.Op(line[6:]))
            elif line.startswith('cenop '):
                cur['cenops'].append(gemmi.Op(line[6:]))
    return data

def read_ref():
    hall_ref = {}
    with open('../tools/hall-symbols.txt') as f:
        for line in f:
            hm = line[11:25].strip().replace(':r', ':R').replace(':h', ':H')
            assert hm not in hall_ref
            num = int(line[:11].split(':')[0])
            hall_ref[hm] = (num, line[25:].strip())
    return hall_ref

def main():
    syminfo = parse_syminfo(sys.argv[1])
    seen_nums = set()
    for entry in syminfo:
        ccp4 = entry['ccp4']
        hall = entry['hall']
        basisop = entry['basisop']
        num = entry['number']
        xhm = entry['xhm']
        if num not in seen_nums:
            seen_nums.add(num)
            assert ccp4 == num
            if basisop != 'x,y,z':
                print(num, xhm, basisop)
        assert ccp4 == 0 or ccp4 % 1000 == num
        if ccp4 == num:
            if basisop != 'x,y,z':
                pass  # print(ccp4, basisop)
        if basisop != 'x,y,z':
            if '(%s)' % basisop not in hall:
                print('Hall symbol "%s" w/o basisop: %s' % (hall, basisop))
        hall_ops = gemmi.symops_from_hall(hall)
        assert len(hall_ops.cen_ops) == len(entry['cenops'])
        assert (set(gemmi.Op().translated(tr) for tr in hall_ops.cen_ops)
                == set(entry['cenops']))
        assert len(hall_ops.sym_ops) == len(entry['symops'])
        # symops differ in about dozen cases but are the same modulo
        # centering vectors
        generated = set(hall_ops)
        given = set(s * c for s in entry['symops'] for c in entry['cenops'])
        assert generated == given, entry
    print('OK. %d entries.' % len(syminfo))

    hall_ref = read_ref()
    xhm_set = set()
    for d in syminfo:
        xhm = d['xhm']
        if xhm and xhm in xhm_set:
            print('dup xHM:', xhm)
        xhm_set.add(xhm)
        ref = hall_ref.get(xhm)
        if ref is not None:
            assert d['number'] == ref[0]
            hall1 = d['hall']
            hall2 = ref[1]
            sym1 = gemmi.symops_from_hall(hall1)
            sym2 = gemmi.symops_from_hall(hall2)
            assert set(sym1) == set(sym2), (hall1, hall2)
        else:
            print('extra:', xhm, '  (%d)' % d['number'], d['hall'])
    missing = set(hall_ref.keys()) - xhm_set
    for d in sorted(missing, key=lambda m: hall_ref[m][0]):
        print('missing:', d, '  (%d)' % hall_ref[d][0])


main()
