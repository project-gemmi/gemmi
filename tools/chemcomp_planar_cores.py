#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict, deque
from pathlib import Path


def load_gemmi(repo_root: Path):
    sys.path.insert(0, str(repo_root / 'build' / 'py'))
    import gemmi  # type: ignore
    return gemmi


def parse_args():
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description=(
            'Detect merged planar heavy-atom cores '
            'from chemcomp plane restraints.'
        ))
    parser.add_argument(
        'inputs',
        nargs='+',
        help='Chemcomp CIF files or directories.')
    parser.add_argument(
        '--recursive',
        action='store_true',
        help='Recurse into directories.')
    parser.add_argument('--min-size', type=int, default=4,
                        help='Minimum heavy-atom core size.')
    parser.add_argument(
        '--show-all',
        action='store_true',
        help='Show all plane-derived components, not only largest ones.')
    return parser.parse_args(), repo_root


def iter_cifs(inputs, recursive):
    for item in inputs:
        path = Path(item)
        if path.is_dir():
            pattern = '**/*.cif' if recursive else '*.cif'
            for cif in sorted(path.glob(pattern)):
                yield cif
        else:
            yield path


def atom_index_map(cc):
    return {atom.id: idx for idx, atom in enumerate(cc.atoms)}


def heavy_atom_ids(cc):
    return {atom.id for atom in cc.atoms if not atom.is_hydrogen()}


def bond_graph(cc):
    graph = defaultdict(set)
    for bond in cc.rt.bonds:
        a = bond.id1.atom
        b = bond.id2.atom
        graph[a].add(b)
        graph[b].add(a)
    return graph


def plane_groups(cc, heavy_ids):
    groups = []
    for plane in cc.rt.planes:
        ids = []
        seen = set()
        for atom_id in plane.ids:
            name = atom_id.atom
            if name in heavy_ids and name not in seen:
                ids.append(name)
                seen.add(name)
        if len(ids) >= 3:
            groups.append((plane.label, ids))
    return groups


def connected_components(graph):
    seen = set()
    comps = []
    for node in sorted(graph):
        if node in seen:
            continue
        comp = []
        q = deque([node])
        seen.add(node)
        while q:
            cur = q.popleft()
            comp.append(cur)
            for nb in sorted(graph[cur]):
                if nb not in seen:
                    seen.add(nb)
                    q.append(nb)
        comps.append(sorted(comp))
    return comps


def detect_planar_cores(cc, min_size):
    heavy_ids = heavy_atom_ids(cc)
    bonds = bond_graph(cc)
    groups = plane_groups(cc, heavy_ids)

    plane_graph = defaultdict(set)
    plane_membership = defaultdict(list)
    for label, ids in groups:
        for atom in ids:
            plane_membership[atom].append(label)
        for i, a in enumerate(ids):
            for b in ids[i + 1:]:
                if b in bonds[a]:
                    plane_graph[a].add(b)
                    plane_graph[b].add(a)

    cores = []
    for atoms in connected_components(plane_graph):
        if len(atoms) < min_size:
            continue
        atom_set = set(atoms)
        core_planes = []
        for label, ids in groups:
            ids_set = set(ids)
            if len(ids_set & atom_set) >= 2:
                core_planes.append(label)
        articulation = [a for a in atoms if len(set(plane_membership[a])) > 1]
        anchors = []
        for atom in atoms:
            outside = sorted(nb for nb in bonds[atom] if nb not in atom_set)
            if outside:
                anchors.append((atom, outside))
        inner_bonds = []
        seen_pairs = set()
        for a in atoms:
            for b in bonds[a]:
                if b in atom_set:
                    key = tuple(sorted((a, b)))
                    if key not in seen_pairs:
                        seen_pairs.add(key)
                        inner_bonds.append(key)
        cores.append({
            'atoms': atoms,
            'planes': sorted(set(core_planes)),
            'articulation': articulation,
            'anchors': anchors,
            'inner_bonds': sorted(inner_bonds),
            'plane_membership': {
                a: sorted(set(plane_membership[a]))
                for a in atoms
            },
        })
    cores.sort(key=lambda c: (-len(c['atoms']), -len(c['planes']), c['atoms']))
    return groups, cores


def summarize_core(core):
    lines = []
    lines.append(
        (
            '  core atoms={:d} planes={:d} bonds={:d} '
            'articulation={:d} anchors={:d}'
        ).format(
            len(
                core['atoms']), len(
                core['planes']), len(
                    core['inner_bonds']), len(
                        core['articulation']), len(
                            core['anchors'])))
    lines.append('    atoms: ' + ' '.join(core['atoms']))
    lines.append('    planes: ' + ' '.join(core['planes']))
    if core['articulation']:
        lines.append('    articulation: ' + ' '.join(core['articulation']))
        for atom in core['articulation']:
            lines.append(
                '      {} <= {}'.format(
                    atom, ', '.join(
                        core['plane_membership'][atom])))
    if core['anchors']:
        lines.append('    anchors:')
        for atom, outside in core['anchors']:
            lines.append('      {} -> {}'.format(atom, ', '.join(outside)))
    return '\n'.join(lines)


def main():
    args, repo_root = parse_args()
    gemmi = load_gemmi(repo_root)
    for path in iter_cifs(args.inputs, args.recursive):
        doc = gemmi.cif.read(str(path))
        cc = gemmi.make_chemcomp_from_block(doc.sole_block())
        groups, cores = detect_planar_cores(cc, args.min_size)
        print('{}  {}  atoms={} bonds={} planes={}'.format(
            cc.name, path, len(cc.atoms), len(cc.rt.bonds), len(groups)))
        if not cores:
            print('  no planar cores >= {}'.format(args.min_size))
            print()
            continue
        use = cores if args.show_all else cores[:1]
        for core in use:
            print(summarize_core(core))
        print()


if __name__ == '__main__':
    main()
