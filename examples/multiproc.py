# Example of multiprocessing: travers subdirectories and process asynchronously
# coordinate files using 4 processes in parallel.

import multiprocessing
import sys
import gemmi

def f(path):
    st = gemmi.read_structure(path)
    weight = st[0].calculate_mass()
    return (st.name, weight)

def main(top_dir):
    with multiprocessing.Pool(processes=4) as pool:
        it = pool.imap_unordered(f, gemmi.CoorFileWalk(top_dir))
        for (name, weight) in it:
            print('%s %.1f kDa' % (name, weight / 1000))

if __name__ == '__main__':
    main(sys.argv[1])
