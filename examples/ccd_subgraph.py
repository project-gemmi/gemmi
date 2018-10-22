#!/usr/bin/env python3

# List CCD entries that contain the specified entry as a substructure.
# Ignoring hydrogens and bond types.

import sys
import networkx
from networkx.algorithms import isomorphism
import gemmi

CCD_PATH = 'components.cif.gz'

def graph_from_block(block):
    cc = gemmi.make_chemcomp_from_block(block)
    cc.remove_hydrogens()
    G = networkx.Graph()
    for atom in cc.atoms:
        G.add_node(atom.id, Z=atom.el.atomic_number)
    for bond in cc.rt.bonds:
        G.add_edge(bond.id1.atom, bond.id2.atom)
    return G

def main():
    assert len(sys.argv) == 2, "Usage: ccd_subgraph.py three-letter-code"
    ccd = gemmi.cif.read(CCD_PATH)
    entry = sys.argv[1]
    pattern = graph_from_block(ccd[entry])
    pattern_nodes = networkx.number_of_nodes(pattern)
    pattern_edges = networkx.number_of_edges(pattern)
    node_match = isomorphism.categorical_node_match('Z', 0)
    for block in ccd:
        G = graph_from_block(block)
        GM = isomorphism.GraphMatcher(G, pattern, node_match=node_match)
        if GM.subgraph_is_isomorphic():
            print(block.name, '\t +%d nodes, +%d edges' % (
                networkx.number_of_nodes(G) - pattern_nodes,
                networkx.number_of_edges(G) - pattern_edges))

main()
