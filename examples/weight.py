#!/usr/bin/env python
# This script reads mmCIF file(s) and compares weights:
# * _entity.formula_weight with _entity_poly_seq and _chem_comp.formula_weight,
# * _chem_comp.formula_weight with _chem_comp.formula and atomic weights.

import os
import argparse
import collections
import itertools
from gemmi import cif, Element, expand_if_pdb_code


H2O_MASS = Element('O').weight + 2 * Element('H').weight
# but why PO2 and not PO3Hx?
PO2_MASS = Element('P').weight + 2 * Element('O').weight


# this function is also used in examples/monomers.py
def formula_to_dict(formula):
    '"O4 P -3" -> {O:4, P:1}'
    fdict = {}
    for elnum in formula.split():
        na = sum(e.isalpha() for e in elnum)
        if na == len(elnum):
            fdict[elnum] = 1
        elif na != 0:
            fdict[elnum[:na]] = int(elnum[na:])
    return fdict

def check_chem_comp_formula_weight(block):
    for cc in block.find('_chem_comp.', ['id', 'formula', 'formula_weight']):
        if cc.str(0) in ('UNX', 'UNL'):  # unknown residue or ligand
            continue
        fdict = formula_to_dict(cc.str(1).title())
        calc_weight = sum(n * Element(e).weight for (e, n) in fdict.items())
        diff = calc_weight - cif.as_number(cc[2])
        if not (abs(diff) < 0.1):  # also true if diff is NaN
            print('%s %s  %-16s % 9.3f - %9s = %+.3f' %
                  (block.name, cc[0], cc.str(1), calc_weight, cc[2], diff))


def check_entity_formula_weight(block):
    # read chem_comp weights from the file, e.g. THR: 119.120
    cc_weights = {cc.str(0): cif.as_number(cc[1]) for cc in
                  block.find('_chem_comp.', ['id', 'formula_weight'])}

    # read polymer types from _entity_poly.type
    poly_types = {ep.str(0): ep.str(1) for ep in
                  block.find('_entity_poly.', ['entity_id', 'type'])}

    # read and sort sequences from _entity_poly_seq
    entity_seq = collections.defaultdict(list)
    for m in block.find('_entity_poly_seq.', ['entity_id', 'num', 'mon_id']):
        entity_seq[m.str(0)].append((cif.as_int(m[1]), cif.as_string(m[2])))
    for seq in entity_seq.values():
        seq.sort(key=lambda p: p[0])

    # read _entity.formula_weight values
    f_weights = block.find('_entity.', ['id', 'formula_weight'])
    entity_weights = {ent.str(0): cif.as_number(ent[1]) for ent in f_weights}

    # calculate weight from sequences and compare with _entity.formula_weight
    for ent, seq in entity_seq.items():
        # in case of microheterogeneity take the first conformer
        main_seq = [next(g) for _, g in itertools.groupby(seq, lambda p: p[0])]
        weight = sum(cc_weights[mon_id] for (num, mon_id) in main_seq)
        weight -= (len(main_seq) - 1) * H2O_MASS
        if 'ribonucleotide' in poly_types[ent]:  # DNA or RNA
            weight -= PO2_MASS
        diff = weight - entity_weights[ent]
        if abs(diff) > max(0.1, 3e-5 * weight):
            print('%4s entity_id: %2s %10.2f - %10.2f = %+8.3f' %
                  (block.name, ent, weight, entity_weights[ent], diff))

# the same function as in examples/matthews.py
def get_file_paths_from_args():
    """\
    Process arguments as filenames or directories with .cif(.gz) files,
    and yield the file paths.
    Normally we first test our scripts on a few files:
      ./myscript 1mru.cif another.cif
    and then do pdb-wide analysis:
      ./myscript $PDB_DIR/structures/divided/mmCIF
    If $PDB_DIR is set you can check the specifies PDB entries:
      ./myscript 1ABC 2def
    If you have a list of PDB codes to analyze (one code per line, the code
    must be the first word, but may be followed by others), do:
      ./myscript --only=my-list.txt $PDB_DIR/structures/divided/mmCIF
    """
    parser = argparse.ArgumentParser(usage='%(prog)s [options] path [...]')
    parser.add_argument('path', nargs='+', help=argparse.SUPPRESS)
    parser.add_argument('--only', metavar='LIST',
                        help='Use only files that match names in this file')
    args = parser.parse_args()
    only = None
    if args.only:
        with open(args.only) as list_file:
            only = set(line.split()[0].lower() for line in list_file
                       if line.strip())
    for arg in args.path:
        if os.path.isdir(arg):
            for root, dirs, files in os.walk(arg):
                dirs.sort()
                for name in sorted(files):
                    for ext in ['.cif', '.cif.gz']:
                        if name.endswith(ext):
                            if not only or name[:-len(ext)].lower() in only:
                                yield os.path.join(root, name)
                                break
        else:
            yield expand_if_pdb_code(arg)


def main():
    for path in get_file_paths_from_args():
        block = cif.read(path).sole_block()
        check_chem_comp_formula_weight(block)
        check_entity_formula_weight(block)

if __name__ == '__main__':
    main()
