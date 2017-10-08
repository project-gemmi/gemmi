#!/usr/bin/env python
import collections
import itertools
from gemmi import cif, Element
import util


H2O_MASS = Element('O').weight + 2 * Element('H').weight
# but why PO2 and not PO3Hx?
PO2_MASS = Element('P').weight + 2 * Element('O').weight


def check_chem_comp_formula_weight(block):
    for cc in block.find('_chem_comp.', ['id', 'formula', 'formula_weight']):
        if cc.str(0) in ('UNX', 'UNL'):  # unknown residue or ligand
            continue
        fdict = util.formula_to_dict(cc.str(1).title())
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


def main():
    for path in util.get_file_paths_from_args():
        block = cif.read(path).sole_block()
        check_chem_comp_formula_weight(block)
        check_entity_formula_weight(block)

if __name__ == '__main__':
    main()
