# Gathers data for Ramachandran plot that can be plotted with rama_plot.py.
# See the section about torsion angles in documentation.

import sys
from math import degrees
import gemmi

ramas = {aa: [] for aa in [
    'LEU', 'ALA', 'GLY', 'VAL', 'GLU', 'SER', 'LYS', 'ASP', 'THR', 'ILE',
    'ARG', 'PRO', 'ASN', 'PHE', 'GLN', 'TYR', 'HIS', 'MET', 'CYS', 'TRP']}

for path in gemmi.CoorFileWalk(sys.argv[1]):
    st = gemmi.read_structure(path)
    if 0.1 < st.resolution < 1.5:
        model = st[0]
        for chain in model:
            for res in chain.get_polymer():
                # previous_residue() and next_residue() return previous/next
                # residue only if the residues are bonded. Otherwise -- None.
                prev_res = chain.previous_residue(res)
                next_res = chain.next_residue(res)
                if prev_res and next_res and next_res.name != 'PRO':
                    v = gemmi.calculate_phi_psi(prev_res, res, next_res)
                    try:
                        ramas[res.name].append(v)
                    except KeyError:
                        pass

# Write data to files
for aa, data in ramas.items():
    with open('ramas/' + aa + '.tsv', 'w') as f:
        for phi, psi in data:
            f.write('%.4f\t%.4f\n' % (degrees(phi), degrees(psi)))
