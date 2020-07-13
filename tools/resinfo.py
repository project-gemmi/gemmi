#!/usr/bin/env python3

# Read CCD and fill-in one-letter-codes in gemmi/resinfo.hpp
# Usage: ./tools/resinfo.py > resinfo.hpp-new

import re
from sys import stderr
import gemmi
from gemmi import cif

ccd = cif.read('components.cif.gz')

STANDARD = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLN', 'GLU',
            'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
            'PRO', 'SER', 'THR', 'TRP', 'TYR', 'UNK', 'VAL',
            'SEC', 'PYL',
            'A', 'C', 'G', 'I', 'U',
            'DA', 'DC', 'DG', 'DI', 'DT', 'DU']

def calculate_formula_weight(formula):
    total = 0.
    for elem_count in formula.split():
        if elem_count.isalpha():
            elem = elem_count
            count = 1
        else:
            n = 2 if elem_count[1].isalpha() else 1
            elem = elem_count[:n]
            count = int(elem_count[n:])
        total += count * gemmi.Element(elem).weight
    return total


pattern = re.compile(r'case ID\("([A-Z0-9]+)"\): '
                     r"return .*'(.)',( +\d+), [\d.]+f ")
#pattern = re.compile(r"case '([A-Z])': return .*'(.)',( +\d+), [\d.]+f ")
for line in open('include/gemmi/resinfo.hpp'):
    m = pattern.search(line)
    if m:
        name = m.group(1)
        #name = 'D' + name
        block = ccd[name]

        # modify one letter code
        olc = block.find_value('_chem_comp.one_letter_code')
        if olc != '?':
            code = olc[0]
            if name not in STANDARD:
                code = code.lower()
            pos = m.start(2)
            line = line[:pos] + code + line[pos+1:]
        if len(olc) > 1:
            print('One-letter-codes for %s: %s' % (name, olc), file=stderr)
        parent = block.find_value('_chem_comp.mon_nstd_parent_comp_id')
        if parent and parent != '?' and ',' not in parent:
            parent_olc = ccd[parent].find_value('_chem_comp.one_letter_code')
            if parent_olc != olc:
                print(name, olc, '<>', parent, parent_olc, file=stderr)

        # modify hydrogen count
        old_nh = int(m.group(3))
        formula = cif.as_string(block.find_value('_chem_comp.formula'))
        obj = re.search(r'\b[HD](\d+)\b', formula)
        nh = int(obj.group(1)) if obj else 0
        weight = calculate_formula_weight(formula)
        line_middle = "', %3d, %#.6gf" % (nh, weight)
        line = line[:m.end(2)] + line_middle + line[m.end(3):]

        # check residue type
        res_type = cif.as_string(block.find_value('_chem_comp.type')).upper()
        has_AAD = 'RI::AAD' in line
        expected_AAD = res_type.startswith('D-PEPTIDE')
        if expected_AAD and not has_AAD:
            print('D-PEPTIDE:', name, file=stderr)
        if not expected_AAD and has_AAD:
            print('not D-PEPTIDE:', name, file=stderr)

    print(line, end='')
