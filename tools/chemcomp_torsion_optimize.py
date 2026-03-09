#!/usr/bin/env python3

import argparse
import math
import sys
from collections import deque
from pathlib import Path


def load_gemmi(repo_root: Path):
    sys.path.insert(0, str(repo_root / 'build' / 'py'))
    import gemmi  # type: ignore
    return gemmi


def parse_args():
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description='Prototype torsion optimizer for ChemComp XYZ generation.')
    parser.add_argument('inputs', nargs='+', help='Prepared monomer CIF files.')
    parser.add_argument('--passes', type=int, default=4, help='Greedy passes over rotatable bonds.')
    parser.add_argument('--max-subtree', type=int, default=24, help='Skip moved subtrees larger than this.')
    parser.add_argument('--verbose', action='store_true', help='Print accepted bond rotations.')
    return parser.parse_args(), repo_root


def finite(atom):
    return math.isfinite(atom.xyz.x) and math.isfinite(atom.xyz.y) and math.isfinite(atom.xyz.z)


def is_heavy(atom):
    return not atom.is_hydrogen()


def clone_pos(gemmi, pos):
    return gemmi.Position(pos.x, pos.y, pos.z)


def snapshot_positions(gemmi, cc, indices):
    return {idx: clone_pos(gemmi, cc.atoms[idx].xyz) for idx in indices}


def restore_positions(cc, saved):
    for idx, pos in saved.items():
        cc.atoms[idx].xyz = pos


def build_graphs(cc):
    atom_index = {a.id: i for i, a in enumerate(cc.atoms)}
    adjacency = [[] for _ in cc.atoms]
    bond_map = {}
    for bond in cc.rt.bonds:
        i = atom_index.get(bond.id1.atom)
        j = atom_index.get(bond.id2.atom)
        if i is None or j is None:
            continue
        adjacency[i].append(j)
        adjacency[j].append(i)
        bond_map[tuple(sorted((i, j)))] = bond
    return atom_index, adjacency, bond_map


def bond_target(bond):
    return bond.value if math.isfinite(bond.value) else bond.value_nucleus


def bond_z(cc, bond_map, i, j):
    bond = bond_map.get(tuple(sorted((i, j))))
    if bond is None or not finite(cc.atoms[i]) or not finite(cc.atoms[j]):
        return 0.0
    target = bond_target(bond)
    if not math.isfinite(target):
        return 0.0
    esd = bond.esd if math.isfinite(bond.esd) and bond.esd > 0 else 0.02
    return abs(cc.atoms[i].xyz.dist(cc.atoms[j].xyz) - target) / esd


def angle_diff_deg(actual, target):
    return abs(math.remainder(actual - target, 360.0))


def plane_from_positions(positions):
    if len(positions) < 3:
        return None
    cx = sum(p.x for p in positions) / len(positions)
    cy = sum(p.y for p in positions) / len(positions)
    cz = sum(p.z for p in positions) / len(positions)
    nx = ny = nz = 0.0
    for i, a in enumerate(positions):
        b = positions[(i + 1) % len(positions)]
        nx += (a.y - b.y) * (a.z + b.z)
        ny += (a.z - b.z) * (a.x + b.x)
        nz += (a.x - b.x) * (a.y + b.y)
    norm = math.sqrt(nx * nx + ny * ny + nz * nz)
    if norm < 1e-10:
        return None
    return (cx, cy, cz), (nx / norm, ny / norm, nz / norm)


def chiral_penalty(gemmi, cc, atom_index):
    penalty = 0.0
    for chir in cc.rt.chirs:
        if chir.sign == gemmi.ChiralityType.Both:
            continue
        ids = [atom_index.get(chir.id_ctr.atom), atom_index.get(chir.id1.atom),
               atom_index.get(chir.id2.atom), atom_index.get(chir.id3.atom)]
        if None in ids:
            continue
        atoms = [cc.atoms[i] for i in ids]
        if not all(finite(a) for a in atoms):
            continue
        v1 = atoms[1].xyz - atoms[0].xyz
        v2 = atoms[2].xyz - atoms[0].xyz
        v3 = atoms[3].xyz - atoms[0].xyz
        vol = v1.dot(v2.cross(v3))
        if chir.is_wrong(vol):
            penalty += 500.0
    return penalty


def objective(gemmi, cc, atom_index, adjacency, bond_map):
    score = 0.0
    worst_bond = (0.0, '')
    for (i, j), bond in bond_map.items():
        if not (is_heavy(cc.atoms[i]) and is_heavy(cc.atoms[j])):
            continue
        z = bond_z(cc, bond_map, i, j)
        score += min(z, 100.0) ** 2
        if z > worst_bond[0]:
            worst_bond = (z, f'{cc.atoms[i].id}-{cc.atoms[j].id}')
    for angle in cc.rt.angles:
        ids = [atom_index.get(angle.id1.atom), atom_index.get(angle.id2.atom), atom_index.get(angle.id3.atom)]
        if None in ids or not math.isfinite(angle.value):
            continue
        i, j, k = ids
        if not (finite(cc.atoms[i]) and finite(cc.atoms[j]) and finite(cc.atoms[k])):
            continue
        esd = angle.esd if math.isfinite(angle.esd) and angle.esd > 0 else 2.0
        actual = math.degrees(gemmi.calculate_angle(cc.atoms[i].xyz, cc.atoms[j].xyz, cc.atoms[k].xyz))
        score += (angle_diff_deg(actual, angle.value) / esd) ** 2
    for plane in cc.rt.planes:
        pts = []
        for atom_id in plane.ids:
            idx = atom_index.get(atom_id.atom)
            if idx is None or not finite(cc.atoms[idx]):
                pts = []
                break
            pts.append(cc.atoms[idx].xyz)
        fit = plane_from_positions(pts)
        if fit is None:
            continue
        centroid, normal = fit
        esd = plane.esd if math.isfinite(plane.esd) and plane.esd > 0 else 0.02
        for p in pts:
            dist = abs((p.x - centroid[0]) * normal[0] + (p.y - centroid[1]) * normal[1] + (p.z - centroid[2]) * normal[2])
            score += (dist / esd) ** 2
    score += chiral_penalty(gemmi, cc, atom_index)
    return score, worst_bond


def detect_small_cycles(adjacency):
    cycles = []
    seen = set()
    for a in range(len(adjacency)):
        for b in adjacency[a]:
            if a >= b:
                continue
            parent = {a: None}
            q = deque([a])
            while q and b not in parent:
                cur = q.popleft()
                for nb in adjacency[cur]:
                    if (cur == a and nb == b) or (cur == b and nb == a):
                        continue
                    if nb in parent:
                        continue
                    parent[nb] = cur
                    q.append(nb)
            if b in parent:
                path = []
                cur = b
                while cur is not None:
                    path.append(cur)
                    cur = parent[cur]
                if 3 <= len(path) <= 6:
                    cyc = tuple(sorted(path))
                    if cyc not in seen:
                        seen.add(cyc)
                        cycles.append(set(path))
    return cycles


def is_rotatable_bond(cc, adjacency, bond_map, cycles, i, j):
    bond = bond_map.get(tuple(sorted((i, j))))
    if bond is None:
        return False
    if not (is_heavy(cc.atoms[i]) and is_heavy(cc.atoms[j])):
        return False
    if any(i in cyc and j in cyc for cyc in cycles):
        return False
    if len([nb for nb in adjacency[i] if is_heavy(cc.atoms[nb])]) < 2:
        return False
    if len([nb for nb in adjacency[j] if is_heavy(cc.atoms[nb])]) < 2:
        return False
    return True


def heavy_component_excluding_edge(cc, adjacency, start, a, b):
    q = deque([start])
    seen = {start}
    comp = []
    while q:
        cur = q.popleft()
        comp.append(cur)
        for nb in adjacency[cur]:
            if not is_heavy(cc.atoms[nb]):
                continue
            if (cur == a and nb == b) or (cur == b and nb == a):
                continue
            if nb in seen:
                continue
            seen.add(nb)
            q.append(nb)
    return comp


def smaller_side(cc, adjacency, i, j, max_size):
    a = heavy_component_excluding_edge(cc, adjacency, i, i, j)
    b = heavy_component_excluding_edge(cc, adjacency, j, i, j)
    if len(a) + len(b) != sum(1 for atom in cc.atoms if is_heavy(atom)):
        return None
    comp, pivot, anchor = (a, i, j) if len(a) <= len(b) else (b, j, i)
    if len(comp) <= 1 or len(comp) > max_size:
        return None
    return comp, pivot, anchor


def rotate_vec(gemmi, v, axis, angle):
    v = gemmi.Vec3(v.x, v.y, v.z)
    axis = gemmi.Vec3(axis.x, axis.y, axis.z)
    c = math.cos(angle)
    s = math.sin(angle)
    return c * v + s * axis.cross(v) + (1.0 - c) * axis.dot(v) * axis


def rotate_subtree(gemmi, cc, component, pivot, anchor, angle):
    axis = cc.atoms[pivot].xyz - cc.atoms[anchor].xyz
    if axis.length() * axis.length() < 1e-8:
        return False
    axis = axis.normalized()
    origin = cc.atoms[pivot].xyz
    for idx in component:
        if idx == pivot:
            continue
        rel = cc.atoms[idx].xyz - origin
        cc.atoms[idx].xyz = origin + gemmi.Position(rotate_vec(gemmi, rel, axis, angle))
    return True


def torsion_targets(gemmi, cc, atom_index, i, j):
    vals = [0, math.pi/3, -math.pi/3, 2*math.pi/3, -2*math.pi/3, math.pi]
    for tor in cc.rt.torsions:
        ids = [atom_index.get(tor.id1.atom), atom_index.get(tor.id2.atom), atom_index.get(tor.id3.atom), atom_index.get(tor.id4.atom)]
        if None in ids or not math.isfinite(tor.value):
            continue
        if (ids[1], ids[2]) == (i, j) or (ids[2], ids[1]) == (i, j):
            vals.append(math.radians(tor.value))
    out = []
    for v in vals:
        if all(abs(math.remainder(v - w, 2*math.pi)) > 1e-4 for w in out):
            out.append(v)
    return out


def optimize_file(gemmi, path, passes, max_subtree, verbose):
    cc = gemmi.make_chemcomp_from_block(gemmi.cif.read(str(path)).sole_block())
    for atom in cc.atoms:
        atom.xyz = gemmi.Position(float('nan'), float('nan'), float('nan'))
    gemmi.generate_chemcomp_xyz_from_restraints(cc)
    atom_index, adjacency, bond_map = build_graphs(cc)
    cycles = detect_small_cycles(adjacency)
    start_score, start_worst = objective(gemmi, cc, atom_index, adjacency, bond_map)
    start_positions = snapshot_positions(gemmi, cc, range(len(cc.atoms)))
    accepted = 0

    for _ in range(passes):
        improved = False
        current_score, _ = objective(gemmi, cc, atom_index, adjacency, bond_map)
        for (i, j), _bond in sorted(bond_map.items()):
            if not is_rotatable_bond(cc, adjacency, bond_map, cycles, i, j):
                continue
            side = smaller_side(cc, adjacency, i, j, max_subtree)
            if side is None:
                continue
            component, pivot, anchor = side
            saved = snapshot_positions(gemmi, cc, component)
            best_local = current_score
            best_positions = None
            for target in torsion_targets(gemmi, cc, atom_index, anchor, pivot):
                for angle in [target, target + math.pi, target - math.pi]:
                    restore_positions(cc, saved)
                    if not rotate_subtree(gemmi, cc, component, pivot, anchor, angle):
                        continue
                    score, _ = objective(gemmi, cc, atom_index, adjacency, bond_map)
                    if score + 1e-6 < best_local:
                        best_local = score
                        best_positions = snapshot_positions(gemmi, cc, component)
            restore_positions(cc, saved)
            if best_positions is not None:
                restore_positions(cc, best_positions)
                actual_score, _ = objective(gemmi, cc, atom_index, adjacency, bond_map)
                if actual_score > current_score - 1e-6:
                    restore_positions(cc, saved)
                    continue
                current_score = actual_score
                accepted += 1
                improved = True
                if verbose:
                    print(f'  accepted {cc.atoms[anchor].id}-{cc.atoms[pivot].id} size={len(component)} score={current_score:.1f}')
        if not improved:
            break

    final_score, final_worst = objective(gemmi, cc, atom_index, adjacency, bond_map)
    if final_score > start_score + 1e-6:
        restore_positions(cc, start_positions)
        final_score, final_worst = objective(gemmi, cc, atom_index, adjacency, bond_map)
    return {
        'code': cc.name,
        'path': str(path),
        'start_score': start_score,
        'final_score': final_score,
        'start_worst': start_worst,
        'final_worst': final_worst,
        'accepted': accepted,
    }


def main():
    args, repo_root = parse_args()
    gemmi = load_gemmi(repo_root)
    for item in args.inputs:
        result = optimize_file(gemmi, Path(item), args.passes, args.max_subtree, args.verbose)
        print(f"{result['code']}  score {result['start_score']:.1f} -> {result['final_score']:.1f}  "
              f"worst_bond {result['start_worst'][0]:.1f}:{result['start_worst'][1]} -> "
              f"{result['final_worst'][0]:.1f}:{result['final_worst'][1]}  accepted={result['accepted']}  {result['path']}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
