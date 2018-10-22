#!/usr/bin/env python3

# Compares graphs of molecules from cif files (Refmac dictionary or similar)
# with CCD entries.

import sys
import networkx
from networkx.algorithms import isomorphism
import gemmi

CCD_PATH = 'components.cif.gz'

def graph_from_chemcomp(cc):
    G = networkx.Graph()
    for atom in cc.atoms:
        G.add_node(atom.id, Z=atom.el.atomic_number)
    for bond in cc.rt.bonds:
        G.add_edge(bond.id1.atom, bond.id2.atom)
    return G

def compare(cc1, cc2):
    s1 = {a.id for a in cc1.atoms}
    s2 = {a.id for a in cc2.atoms}
    b1 = {b.lexicographic_str() for b in cc1.rt.bonds}
    b2 = {b.lexicographic_str() for b in cc2.rt.bonds}
    if s1 == s2 and b1 == b2:
        #print(cc1.name, "the same")
        return
    G1 = graph_from_chemcomp(cc1)
    G2 = graph_from_chemcomp(cc2)
    node_match = isomorphism.categorical_node_match('Z', 0)
    GM = isomorphism.GraphMatcher(G1, G2, node_match=node_match)
    if GM.is_isomorphic():
        print(cc1.name, 'is isomorphic')
        # we could use GM.match(), but here we try to find the shortest diff
        short_diff = None
        for n, mapping in enumerate(GM.isomorphisms_iter()):
            diff = {k: v for k, v in mapping.items() if k != v}
            if short_diff is None or len(diff) < len(short_diff):
                short_diff = diff
            if n == 10000:  # don't spend too much here
                print(' (it may not be the simplest isomorphism)')
                break
        for id1, id2 in short_diff.items():
            print('\t', id1, '->', id2)
    else:
        print(cc1.name, 'differs')
        if s2 - s1:
            print('\tmissing:', ' '.join(s2 - s1))
        if s1 - s2:
            print('\textra:  ', ' '.join(s1 - s2))

def main():
    ccd = gemmi.cif.read(CCD_PATH)
    absent = 0
    for f in sys.argv[1:]:
        block = gemmi.cif.read(f)[-1]
        cc1 = gemmi.make_chemcomp_from_block(block)
        try:
            block2 = ccd[cc1.name]
        except KeyError:
            absent += 1
            #print(cc1.name, 'not in CCD')
            continue
        cc2 = gemmi.make_chemcomp_from_block(block2)
        cc1.remove_hydrogens()
        cc2.remove_hydrogens()
        compare(cc1, cc2)
    if absent != 0:
        print(absent, 'of', len(sys.argv) - 1, 'monomers not found in CCD')

main()
