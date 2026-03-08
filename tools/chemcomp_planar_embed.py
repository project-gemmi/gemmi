#!/usr/bin/env python3

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path

from chemcomp_planar_cores import detect_planar_cores, load_gemmi


EPS = 1e-8


def parse_args():
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description='Prototype 2D embedding for merged planar chemcomp cores.')
    parser.add_argument('inputs', nargs='+', help='Chemcomp CIF files.')
    parser.add_argument('--core-index', type=int, default=0, help='Which planar core to embed.')
    parser.add_argument('--min-size', type=int, default=4, help='Minimum planar core size.')
    return parser.parse_args(), repo_root


class Vec2(object):
    __slots__ = ('x', 'y')

    def __init__(self, x, y):
        self.x = float(x)
        self.y = float(y)

    def __add__(self, other):
        return Vec2(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Vec2(self.x - other.x, self.y - other.y)

    def __mul__(self, scale):
        return Vec2(self.x * scale, self.y * scale)

    __rmul__ = __mul__

    def dot(self, other):
        return self.x * other.x + self.y * other.y

    def length(self):
        return math.hypot(self.x, self.y)

    def normalized(self):
        d = self.length()
        if d < EPS:
            return Vec2(1.0, 0.0)
        return Vec2(self.x / d, self.y / d)

    def perp(self):
        return Vec2(-self.y, self.x)


def bond_value(bond):
    if math.isfinite(bond.value):
        return bond.value
    return bond.value_nucleus


def build_restraint_maps(cc, core_atoms):
    core_set = set(core_atoms)
    bonds = {}
    adjacency = defaultdict(list)
    angles = {}
    for bond in cc.rt.bonds:
        a = bond.id1.atom
        b = bond.id2.atom
        if a in core_set and b in core_set:
            key = tuple(sorted((a, b)))
            bonds[key] = bond_value(bond)
            adjacency[a].append(b)
            adjacency[b].append(a)
    for angle in cc.rt.angles:
        a = angle.id1.atom
        c = angle.id2.atom
        b = angle.id3.atom
        if a in core_set and b in core_set and c in core_set and math.isfinite(angle.value):
            angles[(a, c, b)] = math.radians(angle.value)
            angles[(b, c, a)] = math.radians(angle.value)
    return bonds, adjacency, angles


def get_angle(angles, a, c, b, default_deg=120.0):
    return angles.get((a, c, b), math.radians(default_deg))


def choose_seed_bond(adjacency, bonds):
    best = None
    best_score = None
    for a, b in bonds:
        score = (len(adjacency[a]) + len(adjacency[b]), -abs(len(adjacency[a]) - len(adjacency[b])), a, b)
        if best_score is None or score > best_score:
            best_score = score
            best = (a, b)
    return best


def angle_between(v1, v2):
    d = max(-1.0, min(1.0, v1.normalized().dot(v2.normalized())))
    return math.acos(d)


def circle_intersections(p1, r1, p2, r2):
    dvec = p2 - p1
    d = dvec.length()
    if d < EPS:
        return []
    if d > r1 + r2 + 1e-6:
        return []
    if d < abs(r1 - r2) - 1e-6:
        return []
    ex = dvec * (1.0 / d)
    x = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d)
    y_sq = max(0.0, r1 * r1 - x * x)
    y = math.sqrt(y_sq)
    base = p1 + ex * x
    ey = ex.perp()
    if y < 1e-6:
        return [base]
    return [base + ey * y, base - ey * y]


def score_candidate(atom, cand, placed, adjacency, bonds, angles):
    score = 0.0
    for nb in adjacency[atom]:
        if nb not in placed:
            continue
        dist = (cand - placed[nb]).length()
        target = bonds[tuple(sorted((atom, nb)))]
        score += 20.0 * abs(dist - target)
    placed_nbs = [nb for nb in adjacency[atom] if nb in placed]
    for i, a in enumerate(placed_nbs):
        for b in placed_nbs[i + 1:]:
            actual = angle_between(placed[a] - cand, placed[b] - cand)
            target = get_angle(angles, a, atom, b)
            score += abs(actual - target)
    return score


def place_from_one_neighbor(atom, nb, placed, adjacency, bonds, angles):
    ref = None
    for ref_atom in adjacency[nb]:
        if ref_atom != atom and ref_atom in placed:
            ref = ref_atom
            break
    p_nb = placed[nb]
    d = bonds[tuple(sorted((atom, nb)))]
    if ref is None:
        return [p_nb + Vec2(d, 0.0)]
    v = (placed[ref] - p_nb).normalized()
    theta = get_angle(angles, ref, nb, atom)
    c = math.cos(theta)
    s = math.sin(theta)
    rot1 = Vec2(v.x * c - v.y * s, v.x * s + v.y * c)
    rot2 = Vec2(v.x * c + v.y * s, -v.x * s + v.y * c)
    return [p_nb + rot1 * d, p_nb + rot2 * d]


def embed_core(core_atoms, bonds, adjacency, angles):
    placed = {}
    seed = choose_seed_bond(adjacency, bonds)
    if seed is None:
        return placed
    a, b = seed
    dab = bonds[tuple(sorted((a, b)))]
    placed[a] = Vec2(0.0, 0.0)
    placed[b] = Vec2(dab, 0.0)

    # Try to place one seed neighbor to fix orientation.
    seed_candidates = []
    for center, anchor in ((a, b), (b, a)):
        for nb in adjacency[center]:
            if nb == anchor or nb in placed:
                continue
            theta = get_angle(angles, anchor, center, nb)
            d = bonds[tuple(sorted((center, nb)))]
            base = (placed[anchor] - placed[center]).normalized()
            c = math.cos(theta)
            s = math.sin(theta)
            rot = Vec2(base.x * c - base.y * s, base.x * s + base.y * c)
            seed_candidates.append((nb, placed[center] + rot * d))
            break
        if seed_candidates:
            break
    if seed_candidates:
        nb, pos = seed_candidates[0]
        placed[nb] = pos

    improved = True
    while improved and len(placed) < len(core_atoms):
        improved = False
        for atom in core_atoms:
            if atom in placed:
                continue
            placed_nbs = [nb for nb in adjacency[atom] if nb in placed]
            if not placed_nbs:
                continue
            candidates = []
            if len(placed_nbs) >= 2:
                for i, n1 in enumerate(placed_nbs):
                    for n2 in placed_nbs[i + 1:]:
                        d1 = bonds[tuple(sorted((atom, n1)))]
                        d2 = bonds[tuple(sorted((atom, n2)))]
                        for cand in circle_intersections(placed[n1], d1, placed[n2], d2):
                            candidates.append(cand)
            if not candidates:
                for nb in placed_nbs:
                    candidates.extend(place_from_one_neighbor(atom, nb, placed, adjacency, bonds, angles))
            if not candidates:
                continue
            best = min(candidates, key=lambda cand: score_candidate(atom, cand, placed, adjacency, bonds, angles))
            placed[atom] = best
            improved = True
    return placed


def relax_bonds(core_atoms, placed, bonds, steps=200, step_size=0.25):
    names = [a for a in core_atoms if a in placed]
    if len(names) < 3:
        return placed
    fixed = set(names[:2])
    for _ in range(steps):
        force = dict((a, Vec2(0.0, 0.0)) for a in names)
        for (a, b), target in bonds.items():
            if a not in placed or b not in placed:
                continue
            delta = placed[b] - placed[a]
            dist = delta.length()
            if dist < EPS:
                continue
            err = dist - target
            corr = delta * (step_size * err / dist)
            if a not in fixed:
                force[a] = force[a] + corr
            if b not in fixed:
                force[b] = force[b] - corr
        moved = 0.0
        for a in names:
            if a in fixed:
                continue
            placed[a] = placed[a] + force[a]
            moved += force[a].length()
        if moved < 1e-6:
            break
    return placed


def summarize_embedding(core_atoms, placed, bonds, angles):
    placed_atoms = set(placed)
    bond_res = []
    for (a, b), target in bonds.items():
        if a in placed_atoms and b in placed_atoms:
            actual = (placed[a] - placed[b]).length()
            bond_res.append((abs(actual - target), a, b, actual, target))
    angle_res = []
    seen = set()
    for (a, c, b), target in angles.items():
        key = tuple(sorted((a, b)) + [c]) if False else (a, c, b)
        if (a, c, b) in seen or (b, c, a) in seen:
            continue
        seen.add((a, c, b))
        if a in placed_atoms and b in placed_atoms and c in placed_atoms:
            actual = math.degrees(angle_between(placed[a] - placed[c], placed[b] - placed[c]))
            angle_res.append((abs(actual - math.degrees(target)), a, c, b, actual, math.degrees(target)))
    bond_res.sort(reverse=True)
    angle_res.sort(reverse=True)
    return bond_res, angle_res


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
        idx = min(args.core_index, len(cores) - 1)
        core = cores[idx]
        bonds, adjacency, angles = build_restraint_maps(cc, core['atoms'])
        placed = embed_core(core['atoms'], bonds, adjacency, angles)
        placed = relax_bonds(core['atoms'], placed, bonds)
        print('  core atoms={} placed={} planes={}'.format(len(core['atoms']), len(placed), len(core['planes'])))
        bond_res, angle_res = summarize_embedding(core['atoms'], placed, bonds, angles)
        for diff, a, b, actual, target in bond_res[:5]:
            print('    bond  {:>8} {:>8}  {:.4f} vs {:.4f}'.format(a, b, actual, target))
        for diff, a, c, b, actual, target in angle_res[:5]:
            print('    angle {:>8} {:>8} {:>8}  {:.2f} vs {:.2f}'.format(a, c, b, actual, target))
        print()


if __name__ == '__main__':
    main()
