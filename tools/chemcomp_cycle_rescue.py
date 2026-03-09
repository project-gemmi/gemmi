#!/usr/bin/env python3
import argparse, math, sys
from collections import deque
from pathlib import Path


def load_gemmi(repo_root):
    sys.path.insert(0, str(repo_root / 'build' / 'py'))
    import gemmi
    return gemmi


def finite(atom):
    return math.isfinite(atom.xyz.x) and math.isfinite(atom.xyz.y) and math.isfinite(atom.xyz.z)


def build_graph(cc):
    idx = {a.id:i for i,a in enumerate(cc.atoms)}
    adj=[[] for _ in cc.atoms]
    bonds={}
    for b in cc.rt.bonds:
        i=idx.get(b.id1.atom); j=idx.get(b.id2.atom)
        if i is None or j is None or cc.atoms[i].is_hydrogen() or cc.atoms[j].is_hydrogen():
            continue
        adj[i].append(j); adj[j].append(i)
        bonds[tuple(sorted((i,j)))] = b
    return idx, adj, bonds


def bond_target(b):
    return b.value if math.isfinite(b.value) else b.value_nucleus


def bond_z(cc, bonds, i, j):
    b=bonds[tuple(sorted((i,j)))]
    t=bond_target(b)
    if not math.isfinite(t):
        return 0.0
    esd=b.esd if b.esd>0 and math.isfinite(b.esd) else 0.02
    return abs(cc.atoms[i].xyz.dist(cc.atoms[j].xyz)-t)/esd


def worst_bond(cc, adj, bonds):
    best=(0,None)
    for (i,j),b in bonds.items():
        if not (finite(cc.atoms[i]) and finite(cc.atoms[j])):
            continue
        z=bond_z(cc,bonds,i,j)
        if z>best[0]:
            best=(z,(i,j))
    return best


def alt_path(adj, start, end):
    q=deque([start]); parent={start:None}
    while q and end not in parent:
        cur=q.popleft()
        for nb in adj[cur]:
            if (cur==start and nb==end) or (cur==end and nb==start):
                continue
            if nb not in parent:
                parent[nb]=cur; q.append(nb)
    if end not in parent:
        return None
    path=[]; cur=end
    while cur is not None:
        path.append(cur); cur=parent[cur]
    return list(reversed(path))


def circle_layout(gemmi, cc, cycle, bonds):
    n=len(cycle)
    lens=[]
    for a,b in zip(cycle, cycle[1:]+cycle[:1]):
        bt = bonds.get(tuple(sorted((a,b))))
        lens.append(bond_target(bt) if bt is not None else 1.45)
    per=sum(lens)
    R=max(per/(2*math.pi), max(lens)/1.9)
    angles=[0.0]
    acc=0.0
    for L in lens[:-1]:
        acc += 2*math.asin(max(-1.0,min(1.0,L/(2*R))))
        angles.append(acc)
    centroid=gemmi.Position(0,0,0)
    for idx in cycle:
        centroid += cc.atoms[idx].xyz
    centroid /= len(cycle)
    newpos={}
    for idx,theta in zip(cycle, angles):
        newpos[idx]=gemmi.Position(centroid.x + R*math.cos(theta), centroid.y + R*math.sin(theta), centroid.z)
    return newpos


def evaluate(gemmi, path):
    repo = Path(__file__).resolve().parents[1]
    gemmi = load_gemmi(repo)
    cc = gemmi.make_chemcomp_from_block(gemmi.cif.read(str(path)).sole_block())
    for a in cc.atoms:
        a.xyz = gemmi.Position(float('nan'), float('nan'), float('nan'))
    gemmi.generate_chemcomp_xyz_from_restraints(cc)
    idx, adj, bonds = build_graph(cc)
    before_z, pair = worst_bond(cc, adj, bonds)
    if pair is None or before_z < 10.0:
        return cc.name, before_z, before_z, []
    cycle = alt_path(adj, pair[0], pair[1])
    if cycle is None or len(cycle) < 8 or len(cycle) > 24:
        return cc.name, before_z, before_z, []
    cycle = cycle[:-1]
    newpos = circle_layout(gemmi, cc, cycle, bonds)
    saved={i: gemmi.Position(cc.atoms[i].xyz.x, cc.atoms[i].xyz.y, cc.atoms[i].xyz.z) for i in range(len(cc.atoms))}
    for i,p in newpos.items():
        cc.atoms[i].xyz = p
    try:
        gemmi.refine_chemcomp_xyz(cc)
    except Exception:
        pass
    after_z, pair2 = worst_bond(cc, adj, bonds)
    if after_z > before_z:
        for i,p in saved.items():
            cc.atoms[i].xyz = p
        after_z, pair2 = worst_bond(cc, adj, bonds)
    labels=[cc.atoms[i].id for i in cycle]
    return cc.name, before_z, after_z, labels


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('inputs', nargs='+')
    args=parser.parse_args()
    for item in args.inputs:
        name,before,after,cycle = evaluate(None, Path(item))
        print(name, f'worst_bond_z {before:.2f} -> {after:.2f}', 'cycle_len', len(cycle))
        if cycle:
            print('  cycle', ' '.join(cycle))

if __name__ == '__main__':
    main()
