import gemmi
from rdkit.Chem import AllChem as rd

BOND_TYPE_TO_RDKIT = {
    gemmi.BondType.Unspec:   rd.BondType.UNSPECIFIED,
    gemmi.BondType.Single:   rd.BondType.SINGLE,
    gemmi.BondType.Double:   rd.BondType.DOUBLE,
    gemmi.BondType.Triple:   rd.BondType.TRIPLE,
    gemmi.BondType.Aromatic: rd.BondType.AROMATIC,
    gemmi.BondType.Deloc:    rd.BondType.OTHER,
    gemmi.BondType.Metal:    rd.BondType.OTHER,
}

def chemcomp_to_rdkit(cc: gemmi.ChemComp):
    em = rd.EditableMol(rd.Mol())
    for atom in cc.atoms:
        a = rd.Atom(atom.el.atomic_number)
        a.SetMonomerInfo(rd.AtomPDBResidueInfo(atom.id))
        em.AddAtom(a)
    idx = {atom.id: n for (n, atom) in enumerate(cc.atoms)}
    for bond in cc.rt.bonds:
        order = BOND_TYPE_TO_RDKIT[bond.type]
        if bond.aromatic:
            order = rd.BondType.AROMATIC
        em.AddBond(idx[bond.id1.atom], idx[bond.id2.atom], order=order)
    return em.GetMol()


if __name__ == '__main__':
    import sys
    path = sys.argv[1]
    block = gemmi.cif.read(path)[-1]
    cc = gemmi.make_chemcomp_from_block(block)
    m = chemcomp_to_rdkit(cc)
    print(rd.MolToSmiles(m))
