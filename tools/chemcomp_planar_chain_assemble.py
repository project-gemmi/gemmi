#!/usr/bin/env python3
import argparse
import importlib.util
import math
import sys
from collections import deque
from pathlib import Path


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_args():
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description='Prototype planar-chain assembly keeping broad plane fragments.')
    parser.add_argument('inputs', nargs='+', help='Prepared chemcomp CIF files.')
    parser.add_argument('--min-size', type=int, default=4)
    parser.add_argument('--verbose', action='store_true')
    return parser.parse_args(), repo_root


def heavy_bond_stats(cc):
    idx = {a.id: i for i, a in enumerate(cc.atoms)}
    worst = (0.0, None)
    count10 = 0
    count5 = 0
    for bond in cc.rt.bonds:
        i = idx.get(bond.id1.atom)
        j = idx.get(bond.id2.atom)
        if i is None or j is None or cc.atoms[i].is_hydrogen() or cc.atoms[j].is_hydrogen():
            continue
        if not math.isfinite(cc.atoms[i].xyz.x) or not math.isfinite(cc.atoms[j].xyz.x):
            continue
        target = bond.value if math.isfinite(bond.value) else bond.value_nucleus
        if not math.isfinite(target):
            continue
        esd = bond.esd if math.isfinite(bond.esd) and bond.esd > 0 else 0.02
        actual = cc.atoms[i].xyz.dist(cc.atoms[j].xyz)
        z = abs(actual - target) / esd
        if z > 5:
            count5 += 1
        if z > 10:
            count10 += 1
        if z > worst[0]:
            worst = (z, (bond.id1.atom, bond.id2.atom, actual, target))
    return worst, count5, count10


def assemble_with_seed(pa, cc, fragments, full_core_atoms, seed_label, bridge_planes):
    by_label = {f['label']: f for f in fragments}
    edges = pa.fragment_edges(fragments)
    core_atoms = sorted(full_core_atoms)
    core_bonds, core_adj, core_angles = pa.build_restraint_maps(cc, core_atoms)
    local = {f['label']: pa.embed_fragment(cc, f)[0] for f in fragments}
    global_pos = dict(local[seed_label])
    placed_frags = {seed_label}
    q = deque([seed_label])
    while q:
        cur = q.popleft()
        for nb_label, overlap in edges[cur]:
            if nb_label in placed_frags:
                continue
            local_pos = local[nb_label]
            fn = None
            if len(overlap) >= 2 and all(atom in global_pos and atom in local_pos for atom in overlap[:2]):
                fn = pa.best_transform(local_pos, global_pos, overlap)
            elif len(overlap) == 1:
                fn = pa.one_overlap_transform(overlap[0], by_label[nb_label], local_pos,
                                              global_pos, core_adj, core_angles, core_bonds)
            if fn is None:
                continue
            for atom, pos in local_pos.items():
                if atom not in global_pos:
                    global_pos[atom] = fn(pos)
            placed_frags.add(nb_label)
            q.append(nb_label)
    if bridge_planes:
        global_pos = pa.place_bridge_planes(bridge_planes, global_pos, core_adj, core_angles, core_bonds)
        global_pos, placed_frags = pa.attach_pending_fragments(fragments, placed_frags, local, global_pos,
                                                               core_adj, core_angles, core_bonds)
    global_pos = pa.place_terminal_leaves(core_atoms, global_pos, core_adj, core_angles, core_bonds)
    return global_pos, placed_frags


def apply_positions(gemmi, cc, positions):
    idx = {a.id: i for i, a in enumerate(cc.atoms)}
    for atom in cc.atoms:
        atom.xyz = gemmi.Position(float('nan'), float('nan'), float('nan'))
    for atom_id, vec in positions.items():
        cc.atoms[idx[atom_id]].xyz = gemmi.Position(vec.x, vec.y, 0.0)


def main():
    args, repo_root = parse_args()
    sys.path.insert(0, str(repo_root / 'tools'))
    pa = load_module(repo_root / 'tools' / 'chemcomp_planar_assemble.py', 'pa_chain')
    pf = load_module(repo_root / 'tools' / 'chemcomp_planar_fragments.py', 'pf_chain')
    gemmi = pa.load_gemmi(repo_root)
    for item in args.inputs:
        path = Path(item)
        cc = gemmi.make_chemcomp_from_block(gemmi.cif.read(str(path)).sole_block())
        _, cores = pa.detect_planar_cores(cc, args.min_size)
        if not cores:
            print(f'{cc.name}  no planar core')
            continue
        core = cores[0]
        fragments = pf.dedup_fragments(
            pf.plane_fragments(cc, core['atoms']) + pf.cycle_fragments(pf.bond_graph(cc, set(core['atoms']))))
        fragments = [f for f in fragments if 4 <= len(f['atoms']) <= 10]
        fragments.sort(key=lambda f: (f['kind'] != 'cycle', -len(f['atoms']), f['label']))
        bridge_planes = []

        best = None
        for frag in fragments:
            positions, placed_frags = assemble_with_seed(pa, cc, fragments, core['atoms'], frag['label'], bridge_planes)
            trial = gemmi.make_chemcomp_from_block(gemmi.cif.read(str(path)).sole_block())
            apply_positions(gemmi, trial, positions)
            before, b5, b10 = heavy_bond_stats(trial)
            try:
                gemmi.refine_chemcomp_xyz(trial)
            except Exception:
                pass
            after, a5, a10 = heavy_bond_stats(trial)
            key = (after[0], a10, a5, -len(positions))
            if best is None or key < best[0]:
                best = (key, frag['label'], len(positions), len(placed_frags), before, after, b5, b10, a5, a10)
        _, seed, placed_atoms, placed_frags, before, after, b5, b10, a5, a10 = best
        print(f'{cc.name}  seed={seed}  placed={placed_atoms}/{len(core["atoms"])} frags={placed_frags}/{len(fragments)}')
        print(f'  before  worst_bond_z={before[0]:.2f}  bond>5={b5} bond>10={b10}')
        if before[1]:
            a, b, actual, target = before[1]
            print(f'    worst {a}-{b} {actual:.4f} vs {target:.4f}')
        print(f'  after   worst_bond_z={after[0]:.2f}  bond>5={a5} bond>10={a10}')
        if after[1]:
            a, b, actual, target = after[1]
            print(f'    worst {a}-{b} {actual:.4f} vs {target:.4f}')
        if args.verbose:
            print('  fragments:', ', '.join(f['label'] for f in fragments))


if __name__ == '__main__':
    main()
