#!/usr/bin/env python3

# This script checks if 3x3 matrices from the MTRIX record (or _struct_ncs_oper
# in mmCIF) are rotation matrices (i.e. orthogonal with det=1).

from __future__ import print_function
import sys
import gemmi

def check_mtrix_rot(path):
    st = gemmi.read_structure(path)
    for ncs in st.ncs:
        mat = ncs.tr.mat
        m_mt = mat.multiply(mat.transpose()).tolist()
        eps = max(abs(m_mt[i][j] - (i == j))
                  for i in [0, 1, 2] for j in [0, 1, 2])
        det = mat.determinant()
        if eps < 0.5e-4 and abs(det - 1) < 0.5e-4:
            status = '    ok'
        else:
            status = '%3s orthogonal +/- %.4f,   det=%+.4f' % (
                     '' if eps < 0.01 else 'NOT', eps, det)
        print('%s %3s %c    %s' % (st.name, ncs.id, 'NY'[ncs.given], status))

def main():
    if len(sys.argv) < 2:
        sys.exit('Specify files, directories or PDB codes.')
    for arg in sys.argv[1:]:
        for path in gemmi.CoorFileWalk(gemmi.expand_if_pdb_code(arg)):
            check_mtrix_rot(path)

if __name__ == '__main__':
    main()
