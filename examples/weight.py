#!/usr/bin/env python
import collections
import itertools
from gemmi import cif
import util


H2O_MASS = 15.9994 + 2 * 1.00794
PO2_MASS = 30.973762 + 2 * 15.9994  # but why PO2 and not PO3Hx?

ELEMENT_MASS = {
    'H': 1.00794, 'He': 4.0026, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811,
    'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'F': 18.998403, 'Ne': 20.1797,
    'Na': 22.98977, 'Mg': 24.305, 'Al': 26.981539, 'Si': 28.0855,
    'P': 30.973761, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,
    'K': 39.0983, 'Ca': 40.078, 'Sc': 44.95591, 'Ti': 47.867, 'V': 50.9415,
    'Cr': 51.9961, 'Mn': 54.93805, 'Fe': 55.845, 'Co': 58.9332,
    'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.723, 'Ge': 72.64,
    'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798,
    'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585, 'Zr': 91.224, 'Nb': 92.9064,
    'Mo': 95.95, 'Tc': 98, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
    'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.76,
    'Te': 127.6, 'I': 126.90447, 'Xe': 131.293, 'Cs': 132.905, 'Ba': 137.327,
    'La': 138.905, 'Ce': 140.116, 'Pr': 140.908, 'Nd': 144.24, 'Pm': 145,
    'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.925, 'Dy': 162.5,
    'Ho': 164.93, 'Er': 167.259, 'Tm': 168.934, 'Yb': 173.05, 'Lu': 174.967,
    'Hf': 178.49, 'Ta': 180.948, 'W': 183.84, 'Re': 186.207, 'Os': 190.23,
    'Ir': 192.217, 'Pt': 195.084, 'Au': 196.967, 'Hg': 200.59, 'Tl': 204.383,
    'Pb': 207.2, 'Bi': 208.98, 'Po': 209, 'At': 210, 'Rn': 222,
    'Fr': 223, 'Ra': 226, 'Ac': 227, 'Th': 232.038, 'Pa': 231.036,
    'U': 238.029, 'Np': 237, 'Pu': 244, 'Am': 243, 'Cm': 247,
    'Bk': 247, 'Cf': 251, 'Es': 252, 'Fm': 257, 'Md': 258,
    'No': 259, 'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271,
    'Bh': 272, 'Hs': 270, 'Mt': 276,
    'D': 2.0141, 'X': 0,
}


def check_chem_comp_formula_weight(block):
    for cc in block.find('_chem_comp.', ['id', 'formula', 'formula_weight']):
        if cc.as_str(0) in ('UNX', 'UNL'):  # unknown residue or ligand
            continue
        fdict = util.formula_to_dict(cc.as_str(1).title())
        calc_weight = sum(n * ELEMENT_MASS[e] for (e, n) in fdict.items())
        diff = calc_weight - cc.as_num(2)
        if not (abs(diff) < 0.1):  # also true if diff is NaN
            print('%s %s  %-16s % 9.3f - %9s = %+.3f' %
                  (block.name, cc[0], cc.as_str(1), calc_weight, cc[2], diff))


def check_entity_formula_weight(block):
    # read chem_comp weights from the file, e.g. THR: 119.120
    cc_weights = {cc.as_str(0): cc.as_num(1) for cc in
                  block.find('_chem_comp.', ['id', 'formula_weight'])}

    # read polymer types from _entity_poly.type
    poly_types = {ep.as_str(0): ep.as_str(1) for ep in
                  block.find('_entity_poly.', ['entity_id', 'type'])}

    # read and sort sequences from _entity_poly_seq
    entity_seq = collections.defaultdict(list)
    for m in block.find('_entity_poly_seq.', ['entity_id', 'num', 'mon_id']):
        entity_seq[m.as_str(0)].append((m.as_int(1), m.as_str(2)))
    for seq in entity_seq.values():
        seq.sort(key=lambda p: p[0])

    # read _entity.formula_weight values
    entity_weights = {ent.as_str(0): ent.as_num(1) for ent in
                      block.find('_entity.', ['id', 'formula_weight'])}

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
        block = cif.read_any(path).sole_block()
        check_chem_comp_formula_weight(block)
        check_entity_formula_weight(block)

if __name__ == '__main__':
    main()
