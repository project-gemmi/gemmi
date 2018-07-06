#!/usr/bin/env python
# This script looks for chains that exceed the size of the unit cell (by >20%)
# in one of the a, b, c directions.

import sys
import gemmi

def run(path):
    counter = 0
    st = gemmi.read_structure(path)
    if st.cell.is_crystal():
        st.add_entity_types()
        for chain in st[0]:
            polymer = chain.get_polymer()
            if polymer:
                low_bounds = [float('+inf')] * 3
                high_bounds = [float('-inf')] * 3
                for residue in polymer:
                    for atom in residue:
                        pos = st.cell.fractionalize(atom.pos)
                        for i in range(3):
                            if pos[i] < low_bounds[i]:
                                low_bounds[i] = pos[i]
                            if pos[i] > high_bounds[i]:
                                high_bounds[i] = pos[i]
                for i in range(3):
                    delta = high_bounds[i] - low_bounds[i]
                    if delta > 1.2:  # 120% of the unit cell size
                        counter += 1
                        code = st.info['_entry.id']
                        print('%s   chain:%s   delta%c = %.3f' %
                              (code, chain.name, ord('X') + i, delta))
    return counter

def main():
    for arg in sys.argv[1:]:
        for path in gemmi.CoorFileWalk(arg):
            run(path)

if __name__ == '__main__':
    main()
