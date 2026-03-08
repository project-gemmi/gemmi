#!/usr/bin/env python3

import argparse
import math
import sys
from collections import defaultdict, deque
from pathlib import Path

from chemcomp_planar_cores import detect_planar_cores, load_gemmi
from chemcomp_planar_fragments import bond_graph, cycle_fragments, dedup_fragments, plane_fragments
from chemcomp_planar_embed import (Vec2, angle_between, build_restraint_maps, embed_core,
                                   get_angle,
                                   optimize_layout, refine_angles, relax_bonds)


EPS = 1e-8


def parse_args():
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description='Prototype planar fragment assembly for chemcomp cores.')
    parser.add_argument('inputs', nargs='+', help='Chemcomp CIF files.')
    parser.add_argument('--core-index', type=int, default=0, help='Which planar core to inspect.')
    parser.add_argument('--min-size', type=int, default=4, help='Minimum planar core size.')
    parser.add_argument('--bond-z-limit', type=float, default=10.0,
                        help='Fail threshold for worst bond Z (default: 10).')
    parser.add_argument('--summary-only', action='store_true',
                        help='Print only the compact status summary.')
    return parser.parse_args(), repo_root


def choose_fragments(cc, core_atoms):
    graph = bond_graph(cc, set(core_atoms))
    frags = dedup_fragments(plane_fragments(cc, core_atoms) + cycle_fragments(graph))
    bad_planes = []
    # Prefer cycles and small articulation planes. Drop broad plane supersets that are just ring+few extras.
    out = []
    cycle_sets = [set(f['atoms']) for f in frags if f['kind'] == 'cycle']
    for frag in frags:
        aset = set(frag['atoms'])
        if frag['kind'] == 'plane':
            drop = False
            for cset in cycle_sets:
                extras = aset - cset
                if cset <= aset and 0 < len(extras) <= 3:
                    drop = True
                    break
            if drop:
                continue
            local_pos, _bonds, adj, _angles = embed_fragment(cc, frag)
            ids = frag['atoms']
            bad = False
            for i, a in enumerate(ids):
                for b in ids[i + 1:]:
                    if b in adj[a]:
                        continue
                    if (local_pos[a] - local_pos[b]).length() < 0.2:
                        bad = True
                        break
                if bad:
                    break
            if bad:
                bad_planes.append(frag)
                continue
        out.append(frag)
    out.sort(key=lambda f: (f['kind'] != 'cycle', -len(f['atoms']), f['label']))
    return out, bad_planes


def fragment_edges(fragments):
    edges = defaultdict(list)
    for i, a in enumerate(fragments):
        aset = set(a['atoms'])
        for j, b in enumerate(fragments[i + 1:], i + 1):
            inter = sorted(aset & set(b['atoms']))
            if inter:
                edges[a['label']].append((b['label'], inter))
                edges[b['label']].append((a['label'], inter))
    return edges


def embed_fragment(cc, frag):
    bonds, adjacency, angles = build_restraint_maps(cc, frag['atoms'])
    placed, _cycles = embed_core(frag['atoms'], bonds, adjacency, angles)
    placed = relax_bonds(frag['atoms'], placed, bonds)
    placed = refine_angles(frag['atoms'], placed, adjacency, bonds, angles)
    placed = optimize_layout(frag['atoms'], placed, bonds, angles)
    return placed, bonds, adjacency, angles


def rigid_transform_from_two_points(src1, src2, dst1, dst2):
    vs = src2 - src1
    vd = dst2 - dst1
    ls = vs.length()
    ld = vd.length()
    if ls < EPS or ld < EPS:
        return None
    us = vs * (1.0 / ls)
    ud = vd * (1.0 / ld)
    cos_t = max(-1.0, min(1.0, us.dot(ud)))
    sin_t = us.x * ud.y - us.y * ud.x

    def apply(p):
        q = p - src1
        r = Vec2(q.x * cos_t - q.y * sin_t, q.x * sin_t + q.y * cos_t)
        return dst1 + r

    return apply


def best_transform(local, global_pos, overlap):
    a, b = overlap[:2]
    cand = []
    for pair in ((a, b), (b, a)):
        fn = rigid_transform_from_two_points(local[pair[0]], local[pair[1]],
                                             global_pos[a], global_pos[b])
        if fn is not None:
            cand.append(fn)
    if not cand:
        return None
    best_fn = None
    best_score = None
    for fn in cand:
        score = 0.0
        for atom in overlap:
            score += (fn(local[atom]) - global_pos[atom]).length()
        if best_score is None or score < best_score:
            best_score = score
            best_fn = fn
    return best_fn


def rotate(v, theta):
    c = math.cos(theta)
    s = math.sin(theta)
    return Vec2(v.x * c - v.y * s, v.x * s + v.y * c)


def local_attach_neighbor(shared, frag, local_pos, core_adj):
    for atom in sorted(core_adj[shared]):
        if atom == shared or atom not in local_pos or atom not in frag['atoms']:
            continue
        return atom
    return None


def one_overlap_transform(shared, frag, local_pos, global_pos, core_adj, core_angles, core_bonds):
    x = local_attach_neighbor(shared, frag, local_pos, core_adj)
    if x is None or shared not in global_pos:
        return None
    refs = [nb for nb in core_adj[shared] if nb in global_pos and nb not in set(frag['atoms'])]
    if not refs:
        refs = [nb for nb in core_adj[shared] if nb in global_pos and nb != x]
    if not refs:
        return None
    ref = refs[0]
    d = core_bonds.get(tuple(sorted((shared, x))))
    if d is None:
        return None
    theta = get_angle(core_angles, ref, shared, x)
    base = (global_pos[ref] - global_pos[shared]).normalized()
    local_base = (local_pos[x] - local_pos[shared]).normalized()
    candidates = []
    for sign in (1.0, -1.0):
        global_dir = rotate(base, sign * theta)
        gx = global_pos[shared] + global_dir * d
        fn = rigid_transform_from_two_points(local_pos[shared], local_pos[x],
                                             global_pos[shared], gx)
        if fn is not None:
            candidates.append(fn)
    if not candidates:
        return None
    best_fn = None
    best_score = None
    for fn in candidates:
        score = 0.0
        for atom in frag['atoms']:
            if atom in global_pos:
                score += (fn(local_pos[atom]) - global_pos[atom]).length()
        for atom in frag['atoms']:
            p = fn(local_pos[atom])
            for nb in core_adj[atom]:
                if nb in global_pos:
                    target = core_bonds.get(tuple(sorted((atom, nb))))
                    if target is not None:
                        score += abs((p - global_pos[nb]).length() - target)
        if best_score is None or score < best_score:
            best_score = score
            best_fn = fn
    return best_fn


def place_terminal_leaves(core_atoms, positions, core_adj, core_angles, core_bonds):
    added = True
    core_set = set(core_atoms)
    while added:
        added = False
        for atom in core_atoms:
            if atom in positions:
                continue
            nbs = [nb for nb in core_adj[atom] if nb in core_set]
            placed_nbs = [nb for nb in nbs if nb in positions]
            if len(nbs) != 1 or len(placed_nbs) != 1:
                continue
            center = placed_nbs[0]
            refs = [nb for nb in core_adj[center] if nb in positions and nb != atom]
            if not refs:
                continue
            ref = refs[0]
            d = core_bonds.get(tuple(sorted((atom, center))))
            if d is None:
                continue
            theta = get_angle(core_angles, ref, center, atom)
            base = (positions[ref] - positions[center]).normalized()
            candidates = []
            for sign in (1.0, -1.0):
                direction = rotate(base, sign * theta)
                candidates.append(positions[center] + direction * d)
            best = None
            best_score = None
            for cand in candidates:
                score = 0.0
                for nb in core_adj[center]:
                    if nb in positions and nb != atom:
                        target = core_angles.get((nb, center, atom))
                        if target is not None:
                            actual = angle_between(positions[nb] - positions[center], cand - positions[center])
                            score += abs(actual - target)
                if best_score is None or score < best_score:
                    best_score = score
                    best = cand
            if best is not None:
                positions[atom] = best
                added = True
    return positions


def place_bridge_planes(planes, positions, core_adj, core_angles, core_bonds):
    added = True
    while added:
        added = False
        for plane in planes:
            for atom in plane['atoms']:
                if atom in positions:
                    continue
                placed_nbs = [nb for nb in core_adj[atom] if nb in positions]
                if len(placed_nbs) >= 2:
                    a, b = placed_nbs[:2]
                    da = core_bonds.get(tuple(sorted((atom, a))))
                    db = core_bonds.get(tuple(sorted((atom, b))))
                    if da is None or db is None:
                        continue
                    dvec = positions[b] - positions[a]
                    d = dvec.length()
                    if d < EPS:
                        continue
                    x = (da * da - db * db + d * d) / (2.0 * d)
                    y_sq = max(0.0, da * da - x * x)
                    y = math.sqrt(y_sq)
                    ex = dvec * (1.0 / d)
                    ey = ex.perp()
                    cands = [positions[a] + ex * x + ey * y,
                             positions[a] + ex * x - ey * y]
                    best = None
                    best_score = None
                    for cand in cands:
                        score = 0.0
                        for ref in core_adj[atom]:
                            if ref in positions:
                                target = core_angles.get((ref, atom, a))
                                if target is not None:
                                    actual = angle_between(positions[ref] - cand, positions[a] - cand)
                                    score += abs(actual - target)
                        if best_score is None or score < best_score:
                            best_score = score
                            best = cand
                    if best is not None:
                        positions[atom] = best
                        added = True
                elif len(placed_nbs) == 1:
                    center = placed_nbs[0]
                    refs = [nb for nb in core_adj[center] if nb in positions and nb != atom]
                    if not refs:
                        continue
                    d = core_bonds.get(tuple(sorted((atom, center))))
                    if d is None:
                        continue
                    best = None
                    best_score = None
                    for ref in refs:
                        theta = get_angle(core_angles, ref, center, atom)
                        base = (positions[ref] - positions[center]).normalized()
                        for sign in (1.0, -1.0):
                            cand = positions[center] + rotate(base, sign * theta) * d
                            score = 0.0
                            for nb in core_adj[atom]:
                                if nb in positions and nb != center:
                                    target = core_bonds.get(tuple(sorted((atom, nb))))
                                    if target is not None:
                                        score += abs((cand - positions[nb]).length() - target)
                            for nb in core_adj[center]:
                                if nb in positions and nb != atom:
                                    target = core_angles.get((nb, center, atom))
                                    if target is not None:
                                        actual = angle_between(positions[nb] - positions[center],
                                                               cand - positions[center])
                                        score += abs(actual - target)
                            if best_score is None or score < best_score:
                                best_score = score
                                best = cand
                    if best is not None:
                        positions[atom] = best
                        added = True
    return positions


def attach_pending_fragments(fragments, placed_frags, local, global_pos, core_adj, core_angles, core_bonds):
    by_label = {f['label']: f for f in fragments}
    changed = True
    while changed:
        changed = False
        for frag in fragments:
            nb_label = frag['label']
            if nb_label in placed_frags:
                continue
            overlap = sorted(atom for atom in frag['atoms'] if atom in global_pos and atom in local[nb_label])
            if not overlap:
                continue
            local_pos = local[nb_label]
            fn = None
            if len(overlap) >= 2:
                fn = best_transform(local_pos, global_pos, overlap)
            else:
                shared = overlap[0]
                fn = one_overlap_transform(shared, by_label[nb_label], local_pos,
                                           global_pos, core_adj, core_angles, core_bonds)
            if fn is None:
                continue
            for atom, pos in local_pos.items():
                if atom in global_pos:
                    continue
                global_pos[atom] = fn(pos)
            placed_frags.add(nb_label)
            changed = True
    return global_pos, placed_frags


def assemble_fragments(cc, fragments, full_core_atoms=None, bridge_planes=None):
    by_label = {f['label']: f for f in fragments}
    edges = fragment_edges(fragments)
    core_atoms = sorted(full_core_atoms if full_core_atoms is not None
                        else {atom for frag in fragments for atom in frag['atoms']})
    core_bonds, core_adj, core_angles = build_restraint_maps(cc, core_atoms)
    local = {}
    for frag in fragments:
        local[frag['label']] = embed_fragment(cc, frag)[0]
    seed = fragments[0]['label']
    global_pos = dict(local[seed])
    placed_frags = {seed}
    q = deque([seed])
    while q:
        cur = q.popleft()
        for nb_label, overlap in edges[cur]:
            if nb_label in placed_frags:
                continue
            if len(overlap) < 2:
                continue
            local_pos = local[nb_label]
            if any(atom not in global_pos or atom not in local_pos for atom in overlap[:2]):
                continue
            fn = best_transform(local_pos, global_pos, overlap)
            if fn is None:
                continue
            for atom, pos in local_pos.items():
                if atom in global_pos:
                    continue
                global_pos[atom] = fn(pos)
            placed_frags.add(nb_label)
            q.append(nb_label)
        for nb_label, overlap in edges[cur]:
            if nb_label in placed_frags or len(overlap) != 1:
                continue
            local_pos = local[nb_label]
            shared = overlap[0]
            fn = one_overlap_transform(shared, by_label[nb_label], local_pos,
                                       global_pos, core_adj, core_angles, core_bonds)
            if fn is None:
                continue
            for atom, pos in local_pos.items():
                if atom in global_pos:
                    continue
                global_pos[atom] = fn(pos)
            placed_frags.add(nb_label)
            q.append(nb_label)
    if bridge_planes:
        global_pos = place_bridge_planes(bridge_planes, global_pos, core_adj, core_angles, core_bonds)
        global_pos, placed_frags = attach_pending_fragments(fragments, placed_frags, local, global_pos,
                                                            core_adj, core_angles, core_bonds)
    global_pos = place_terminal_leaves(core_atoms, global_pos, core_adj, core_angles, core_bonds)
    return global_pos, placed_frags, edges, local


def summarize(cc, core_atoms, positions):
    aset = set(core_atoms)
    _bonds, adjacency, angles = build_restraint_maps(cc, core_atoms)
    bond_res = []
    max_bond_z = 0.0
    for bond in cc.rt.bonds:
        a = bond.id1.atom
        b = bond.id2.atom
        if a in aset and b in aset and a in positions and b in positions:
            target = bond.value if math.isfinite(bond.value) else bond.value_nucleus
            actual = (positions[a] - positions[b]).length()
            diff = abs(actual - target)
            z = diff / bond.esd if bond.esd > 0 else 0.0
            max_bond_z = max(max_bond_z, z)
            bond_res.append((diff, z, a, b, actual, target))
    angle_res = []
    seen = set()
    for (a, c, b), target in angles.items():
        if (b, c, a) in seen:
            continue
        seen.add((a, c, b))
        if a in positions and b in positions and c in positions:
            actual = math.degrees(angle_between(positions[a] - positions[c], positions[b] - positions[c]))
            angle_res.append((abs(actual - math.degrees(target)), a, c, b, actual, math.degrees(target)))
    bond_res.sort(reverse=True)
    angle_res.sort(reverse=True)
    missing = sorted(a for a in core_atoms if a not in positions)
    return bond_res, angle_res, max_bond_z, missing


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
        fragments, bridge_planes = choose_fragments(cc, core['atoms'])
        pos, placed_frags, edges, _local = assemble_fragments(cc, fragments, core['atoms'], bridge_planes)
        bond_res, angle_res, max_bond_z, missing = summarize(cc, core['atoms'], pos)
        status = 'PASS' if not missing and max_bond_z <= args.bond_z_limit else 'FAIL'
        print('  status={}  worst_bond_z={:.2f}  missing={}/{}'.format(
            status, max_bond_z, len(missing), len(core['atoms'])))
        print('  fragments={} placed-fragments={} placed-atoms={}/{}'.format(
            len(fragments), len(placed_frags), len(pos), len(core['atoms'])))
        print('  placed fragment labels: {}'.format(', '.join(sorted(placed_frags))))
        if bridge_planes:
            print('  bridge planes: {}'.format(', '.join(p['label'] for p in bridge_planes)))
        if missing:
            print('  missing atoms: {}'.format(', '.join(missing)))
        if not args.summary_only:
            for diff, z, a, b, actual, target in bond_res[:5]:
                print('    bond  {:>8} {:>8}  {:.4f} vs {:.4f}  z={:.2f}'.format(
                    a, b, actual, target, z))
            for diff, a, c, b, actual, target in angle_res[:5]:
                print('    angle {:>8} {:>8} {:>8}  {:.2f} vs {:.2f}'.format(a, c, b, actual, target))
        print()


if __name__ == '__main__':
    main()
