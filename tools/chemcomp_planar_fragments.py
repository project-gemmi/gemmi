#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict, deque
from pathlib import Path

from chemcomp_planar_cores import detect_planar_cores, load_gemmi


def parse_args():
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description='Decompose planar cores into overlapping plane/ring fragments and articulation graph.')
    parser.add_argument('inputs', nargs='+', help='Chemcomp CIF files.')
    parser.add_argument('--core-index', type=int, default=0, help='Which planar core to inspect.')
    parser.add_argument('--min-size', type=int, default=4, help='Minimum planar core size.')
    return parser.parse_args(), repo_root


def heavy_atoms(cc):
    return {a.id for a in cc.atoms if not a.is_hydrogen()}


def bond_graph(cc, allowed=None):
    graph = defaultdict(set)
    for bond in cc.rt.bonds:
        a = bond.id1.atom
        b = bond.id2.atom
        if allowed is not None and (a not in allowed or b not in allowed):
            continue
        graph[a].add(b)
        graph[b].add(a)
    return graph


def canonical_cycle(cycle):
    seq = list(cycle)
    n = len(seq)
    rots = []
    for i in range(n):
        rots.append(tuple(seq[i:] + seq[:i]))
    rev = list(reversed(seq))
    for i in range(n):
        rots.append(tuple(rev[i:] + rev[:i]))
    return min(rots)


def find_cycles(graph, max_len=8):
    found = set()

    def dfs(start, cur, path, seen):
        for nb in sorted(graph[cur]):
            if nb == start and len(path) >= 3:
                found.add(canonical_cycle(path))
                continue
            if nb in seen or len(path) >= max_len:
                continue
            dfs(start, nb, path + [nb], seen | {nb})

    for start in sorted(graph):
        dfs(start, start, [start], {start})
    return sorted(found, key=lambda cyc: (len(cyc), cyc))


def plane_fragments(cc, core_atoms):
    aset = set(core_atoms)
    frags = []
    for plane in cc.rt.planes:
        ids = []
        seen = set()
        for atom_id in plane.ids:
            name = atom_id.atom
            if name in aset and name not in seen:
                ids.append(name)
                seen.add(name)
        if len(ids) >= 3:
            frags.append({'kind': 'plane', 'label': plane.label, 'atoms': sorted(ids)})
    return frags


def cycle_fragments(core_graph):
    frags = []
    for idx, cyc in enumerate(find_cycles(core_graph), 1):
        frags.append({'kind': 'cycle', 'label': 'cycle-{}'.format(idx), 'atoms': list(cyc)})
    return frags


def dedup_fragments(fragments):
    seen = set()
    out = []
    for frag in fragments:
        key = (frag['kind'], tuple(sorted(frag['atoms'])))
        if key in seen:
            continue
        seen.add(key)
        out.append(frag)
    return out


def articulation_atoms(fragments):
    membership = defaultdict(list)
    for frag in fragments:
        for atom in frag['atoms']:
            membership[atom].append(frag['label'])
    arts = {atom: labs for atom, labs in membership.items() if len(labs) > 1}
    return dict(sorted(arts.items()))


def fragment_graph(fragments):
    edges = []
    for i, a in enumerate(fragments):
        aset = set(a['atoms'])
        for b in fragments[i + 1:]:
            inter = sorted(aset & set(b['atoms']))
            if inter:
                edges.append((a['label'], b['label'], inter))
    return edges


def boundary_atoms(core_atoms, core_graph, fragments):
    aset = set(core_atoms)
    membership = defaultdict(list)
    for frag in fragments:
        for atom in frag['atoms']:
            membership[atom].append(frag['label'])
    out = []
    for atom in sorted(aset):
        deg = len(core_graph[atom])
        if deg <= 1:
            out.append((atom, deg, membership[atom]))
    return out


def main():
    args, repo_root = parse_args()
    gemmi = load_gemmi(repo_root)
    for item in args.inputs:
        path = Path(item)
        cc = gemmi.make_chemcomp_from_block(gemmi.cif.read(str(path)).sole_block())
        _, cores = detect_planar_cores(cc, args.min_size)
        print('{}  {}'.format(cc.name, path))
        if not cores:
            print('  no planar core\n')
            continue
        core = cores[min(args.core_index, len(cores) - 1)]
        core_atoms = core['atoms']
        core_graph = bond_graph(cc, set(core_atoms))
        fragments = dedup_fragments(plane_fragments(cc, core_atoms) + cycle_fragments(core_graph))
        arts = articulation_atoms(fragments)
        edges = fragment_graph(fragments)
        leaves = boundary_atoms(core_atoms, core_graph, fragments)
        print('  core atoms={} fragments={} articulation={} fragment-edges={}'.format(
            len(core_atoms), len(fragments), len(arts), len(edges)))
        print('  fragments:')
        for frag in fragments:
            print('    {:<8} {:<10} {}'.format(frag['kind'], frag['label'], ' '.join(sorted(frag['atoms']))))
        if arts:
            print('  articulation atoms:')
            for atom, labs in arts.items():
                print('    {:<6} {}'.format(atom, ', '.join(labs)))
        if edges:
            print('  fragment graph:')
            for a, b, inter in edges:
                print('    {} -- {} via {}'.format(a, b, ', '.join(inter)))
        if leaves:
            print('  boundary/leaves:')
            for atom, deg, labs in leaves:
                print('    {} deg={} fragments={}'.format(atom, deg, ', '.join(labs)))
        print()


if __name__ == '__main__':
    main()
