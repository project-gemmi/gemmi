#!/usr/bin/env python3

import argparse
import math
import statistics
import sys
from pathlib import Path


def load_gemmi(repo_root: Path):
    sys.path.insert(0, str(repo_root / 'build' / 'py'))
    import gemmi  # type: ignore
    return gemmi


def parse_args():
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description='Regenerate chemcomp ideal coordinates and summarize restraint deviations.')
    parser.add_argument('directory', help='Directory with monomer CIF files.')
    parser.add_argument('--recursive', action='store_true', help='Recurse into subdirectories.')
    parser.add_argument('--limit', type=int, default=0, help='Maximum number of files to process.')
    parser.add_argument('--top', type=int, default=20, help='Number of worst files to print.')
    parser.add_argument('--threshold', type=float, default=1.0,
                        help='Per-file pass threshold in ESD units (default: 1.0).')
    parser.add_argument('--bond-z-mid', type=float, default=5.0,
                        help='Mid catastrophic bond threshold in ESD units (default: 5.0).')
    parser.add_argument('--bond-z-high', type=float, default=10.0,
                        help='High catastrophic bond threshold in ESD units (default: 10.0).')
    parser.add_argument('--verbose', action='store_true', help='Print one line per file.')
    return parser.parse_args(), repo_root


def resolve_paths(root: Path, recursive: bool):
    pattern = '**/*.cif' if recursive else '*.cif'
    return sorted(root.glob(pattern))


class GeometryStats:
    def __init__(self):
        self.count = 0
        self.max_z = 0.0
        self.worst = []
        self.kind_counts = {}
        self.worst_bond_z = 0.0
        self.bonds_over_mid = 0
        self.bonds_over_high = 0

    def add(self, kind: str, label: str, z: float, detail: str):
        if not math.isfinite(z):
            return
        self.count += 1
        self.kind_counts[kind] = self.kind_counts.get(kind, 0) + 1
        if z > self.max_z:
            self.max_z = z
        if kind == 'bond' and z > self.worst_bond_z:
            self.worst_bond_z = z
        self.worst.append((z, kind, label, detail))

    def top(self, n: int):
        return sorted(self.worst, reverse=True)[:n]


def atom_map(cc):
    return {atom.id: atom for atom in cc.atoms}


def finite_xyz(atom) -> bool:
    return math.isfinite(atom.xyz.x) and math.isfinite(atom.xyz.y) and math.isfinite(atom.xyz.z)


def bond_value(bond):
    if math.isfinite(bond.value):
        return bond.value
    return bond.value_nucleus


def plane_from_positions(positions):
    if len(positions) < 3:
        return None
    cx = sum(p.x for p in positions) / len(positions)
    cy = sum(p.y for p in positions) / len(positions)
    cz = sum(p.z for p in positions) / len(positions)
    centroid = (cx, cy, cz)
    nx = ny = nz = 0.0
    for i, a in enumerate(positions):
        b = positions[(i + 1) % len(positions)]
        nx += (a.y - b.y) * (a.z + b.z)
        ny += (a.z - b.z) * (a.x + b.x)
        nz += (a.x - b.x) * (a.y + b.y)
    norm = math.sqrt(nx * nx + ny * ny + nz * nz)
    if norm < 1e-10:
        return None
    return centroid, (nx / norm, ny / norm, nz / norm)


def calc_chiral_volume(a0, a1, a2, a3):
    v1 = a1.xyz - a0.xyz
    v2 = a2.xyz - a0.xyz
    v3 = a3.xyz - a0.xyz
    return v1.dot(v2.cross(v3))


def torsion_diff_deg(gemmi, value_deg: float, target_deg: float, period: int) -> float:
    full = 360.0 / max(1, period)
    return angle_abs_diff_deg(value_deg, target_deg, full)


def angle_abs_diff_deg(value_deg: float, target_deg: float, full: float = 360.0) -> float:
    diff = math.remainder(value_deg - target_deg, full)
    return abs(diff)


def bond_label(bond) -> str:
    return f'{bond.id1.atom}-{bond.id2.atom}'


def angle_label(angle) -> str:
    return f'{angle.id1.atom}-{angle.id2.atom}-{angle.id3.atom}'


def torsion_label(tor) -> str:
    return f'{tor.id1.atom}-{tor.id2.atom}-{tor.id3.atom}-{tor.id4.atom}'


def chirality_label(chir) -> str:
    return f'{chir.id_ctr.atom},{chir.id1.atom},{chir.id2.atom},{chir.id3.atom}'


def evaluate_component(gemmi, path: Path, bond_z_mid: float, bond_z_high: float):
    cc = gemmi.make_chemcomp_from_block(gemmi.cif.read(str(path)).sole_block())
    for atom in cc.atoms:
        atom.xyz = gemmi.Position(float('nan'), float('nan'), float('nan'))
    placed = gemmi.generate_chemcomp_xyz_from_restraints(cc)
    atoms = atom_map(cc)
    stats = GeometryStats()
    missing = []

    for bond in cc.rt.bonds:
        a1 = atoms.get(bond.id1.atom)
        a2 = atoms.get(bond.id2.atom)
        target = bond_value(bond)
        if a1 is None or a2 is None or not finite_xyz(a1) or not finite_xyz(a2) or not math.isfinite(target):
            missing.append(('bond', bond_label(bond)))
            continue
        actual = a1.xyz.dist(a2.xyz)
        z = abs(actual - target) / bond.esd if bond.esd > 0 else float('inf')
        stats.add('bond', bond_label(bond), z, f'{actual:.4f} vs {target:.4f}')
        if math.isfinite(z):
            if z > bond_z_mid:
                stats.bonds_over_mid += 1
            if z > bond_z_high:
                stats.bonds_over_high += 1

    for angle in cc.rt.angles:
        a1 = atoms.get(angle.id1.atom)
        a2 = atoms.get(angle.id2.atom)
        a3 = atoms.get(angle.id3.atom)
        if a1 is None or a2 is None or a3 is None or not finite_xyz(a1) or not finite_xyz(a2) or not finite_xyz(a3) or not math.isfinite(angle.value):
            missing.append(('angle', angle_label(angle)))
            continue
        actual = math.degrees(gemmi.calculate_angle(a1.xyz, a2.xyz, a3.xyz))
        z = angle_abs_diff_deg(actual, angle.value) / angle.esd if angle.esd > 0 else float('inf')
        stats.add('angle', angle_label(angle), z, f'{actual:.4f} vs {angle.value:.4f}')

    for tor in cc.rt.torsions:
        a1 = atoms.get(tor.id1.atom)
        a2 = atoms.get(tor.id2.atom)
        a3 = atoms.get(tor.id3.atom)
        a4 = atoms.get(tor.id4.atom)
        if a1 is None or a2 is None or a3 is None or a4 is None or not finite_xyz(a1) or not finite_xyz(a2) or not finite_xyz(a3) or not finite_xyz(a4) or not math.isfinite(tor.value):
            missing.append(('torsion', torsion_label(tor)))
            continue
        actual = math.degrees(gemmi.calculate_dihedral(a1.xyz, a2.xyz, a3.xyz, a4.xyz))
        diff = torsion_diff_deg(gemmi, actual, tor.value, tor.period)
        z = diff / tor.esd if tor.esd > 0 else 0.0
        stats.add('torsion', torsion_label(tor), z, f'{actual:.4f} vs {tor.value:.4f} p{tor.period}')

    for chir in cc.rt.chirs:
        if chir.sign == gemmi.ChiralityType.Both:
            continue
        a0 = atoms.get(chir.id_ctr.atom)
        a1 = atoms.get(chir.id1.atom)
        a2 = atoms.get(chir.id2.atom)
        a3 = atoms.get(chir.id3.atom)
        if a0 is None or a1 is None or a2 is None or a3 is None or not finite_xyz(a0) or not finite_xyz(a1) or not finite_xyz(a2) or not finite_xyz(a3):
            missing.append(('chirality', chirality_label(chir)))
            continue
        vol = calc_chiral_volume(a0, a1, a2, a3)
        z = 0.0 if not chir.is_wrong(vol) else float('inf')
        stats.add('chirality', chirality_label(chir), z, f'volume={vol:.4f}')

    for plane in cc.rt.planes:
        plane_atoms = []
        for atom_id in plane.ids:
            atom = atoms.get(atom_id.atom)
            if atom is None or not finite_xyz(atom):
                plane_atoms = []
                break
            plane_atoms.append(atom)
        if len(plane_atoms) < 3:
            missing.append(('plane', plane.label))
            continue
        plane_fit = plane_from_positions([atom.xyz for atom in plane_atoms])
        if plane_fit is None:
            missing.append(('plane', plane.label))
            continue
        centroid, normal = plane_fit
        max_dist = 0.0
        for atom in plane_atoms:
            dx = atom.xyz.x - centroid[0]
            dy = atom.xyz.y - centroid[1]
            dz = atom.xyz.z - centroid[2]
            dist = abs(dx * normal[0] + dy * normal[1] + dz * normal[2])
            max_dist = max(max_dist, dist)
        z = max_dist / plane.esd if plane.esd > 0 else float('inf')
        stats.add('plane', plane.label, z, f'max_dist={max_dist:.4f}')

    return {
        'path': path,
        'code': cc.name,
        'placed': placed,
        'atoms': len(cc.atoms),
        'stats': stats,
        'missing': missing,
        'all_within_threshold': placed == len(cc.atoms) and not missing and stats.max_z < 1.0,
    }


def summarize(results, threshold: float, top_n: int,
              bond_z_mid: float, bond_z_high: float):
    processed = len(results)
    fully_placed = sum(1 for r in results if r['placed'] == r['atoms'])
    complete = sum(1 for r in results if not r['missing'])
    within = sum(1 for r in results if r['placed'] == r['atoms'] and not r['missing'] and r['stats'].max_z < threshold)
    max_zs = [r['stats'].max_z for r in results if math.isfinite(r['stats'].max_z)]
    worst_bond_zs = [r['stats'].worst_bond_z for r in results if math.isfinite(r['stats'].worst_bond_z)]
    bond_mid_ok = sum(1 for r in results if r['stats'].bonds_over_mid == 0)
    bond_high_ok = sum(1 for r in results if r['stats'].bonds_over_high == 0)
    total_bonds_over_mid = sum(r['stats'].bonds_over_mid for r in results)
    total_bonds_over_high = sum(r['stats'].bonds_over_high for r in results)
    kind_counts = {}
    for r in results:
        for kind, count in r['stats'].kind_counts.items():
            kind_counts[kind] = kind_counts.get(kind, 0) + count

    print(f'Processed files: {processed}')
    print(f'Fully placed:    {fully_placed}')
    print(f'Complete eval:   {complete}')
    print(f'All diffs < {threshold:g} esd: {within}')
    print(f'No bonds > {bond_z_mid:g} esd: {bond_mid_ok}')
    print(f'No bonds > {bond_z_high:g} esd: {bond_high_ok}')
    if max_zs:
        print(f'Max-z median:    {statistics.median(max_zs):.3f}')
        print(f'Max-z mean:      {statistics.fmean(max_zs):.3f}')
        print(f'Max-z worst:     {max(max_zs):.3f}')
    if worst_bond_zs:
        print(f'Worst bond z median: {statistics.median(worst_bond_zs):.3f}')
        print(f'Worst bond z mean:   {statistics.fmean(worst_bond_zs):.3f}')
        print(f'Worst bond z max:    {max(worst_bond_zs):.3f}')
    print(f'Bonds > {bond_z_mid:g} esd: {total_bonds_over_mid}')
    print(f'Bonds > {bond_z_high:g} esd: {total_bonds_over_high}')
    if kind_counts:
        print('Restraints:      ' + ', '.join(f'{kind}={kind_counts[kind]}' for kind in sorted(kind_counts)))

    ranked = sorted(results,
                    key=lambda r: (r['stats'].bonds_over_high,
                                   r['stats'].worst_bond_z,
                                   r['stats'].bonds_over_mid,
                                   r['stats'].max_z,
                                   len(r['missing'])),
                    reverse=True)
    print(f'\nWorst {min(top_n, len(ranked))} files:')
    for r in ranked[:top_n]:
        missing_note = f", missing={len(r['missing'])}" if r['missing'] else ''
        print(f"{r['code']:>6}  max_z={r['stats'].max_z:8.3f}  worst_bond_z={r['stats'].worst_bond_z:8.3f}  "
              f"bond>{bond_z_mid:g}={r['stats'].bonds_over_mid:>3}  "
              f"bond>{bond_z_high:g}={r['stats'].bonds_over_high:>3}  "
              f"placed={r['placed']:>3}/{r['atoms']:<3}{missing_note}  {r['path']}")
        for z, kind, label, detail in r['stats'].top(3):
            print(f'        {kind:<9} z={z:8.3f}  {label}  {detail}')
        if r['missing']:
            for kind, label in r['missing'][:3]:
                print(f'        missing    {kind} {label}')


def main():
    args, repo_root = parse_args()
    gemmi = load_gemmi(repo_root)
    root = Path(args.directory)
    if not root.is_dir():
        raise SystemExit(f'not a directory: {root}')
    paths = resolve_paths(root, args.recursive)
    if args.limit > 0:
        paths = paths[:args.limit]
    if not paths:
        raise SystemExit('no CIF files found')

    results = []
    for path in paths:
        result = evaluate_component(gemmi, path, args.bond_z_mid, args.bond_z_high)
        results.append(result)
        if args.verbose:
            status = 'OK' if result['placed'] == result['atoms'] and not result['missing'] and result['stats'].max_z < args.threshold else 'WARN'
            print(f"{status} {result['code']}: max_z={result['stats'].max_z:.3f} "
                  f"worst_bond_z={result['stats'].worst_bond_z:.3f} "
                  f"bond>{args.bond_z_mid:g}={result['stats'].bonds_over_mid} "
                  f"bond>{args.bond_z_high:g}={result['stats'].bonds_over_high} "
                  f"placed={result['placed']}/{result['atoms']} missing={len(result['missing'])}")

    summarize(results, args.threshold, args.top, args.bond_z_mid, args.bond_z_high)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
