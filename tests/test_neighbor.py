#!/usr/bin/env python

import unittest
import gemmi
from common import full_path

# In 5a11 applying NCS causes atom clashing
FRAGMENT_5A11 = """\
CRYST1   48.367   89.613   83.842  90.00 101.08  90.00 P 1 21 1      4          
MTRIX1   1 -0.999980  0.006670  0.000050       13.74419                         
MTRIX2   1  0.006260  0.941760 -0.336220       35.56956                         
MTRIX3   1 -0.002290 -0.336210 -0.941780      205.34927                         
ATOM    256  SG  CYS A  37      -1.002 -31.125  88.394  1.00 14.38           S  
ATOM   2969  SG  CYS B  37      14.582 -23.455 132.554  1.00 18.14           S  
"""  # noqa: W291 - trailing whitespace

# In 1gtv two different chains (with partial occupancy) are exactly in the
# same place
FRAGMENT_1GTV = """\
CRYST1   76.353   76.353  134.815  90.00  90.00 120.00 P 65 2 2     24
ATOM    635  SG  CYS A  85      42.948   6.483  17.913  0.48 23.86           S
ATOM   2293  SG  CYS B  85      42.948   6.483  17.913  0.52 23.86           S
"""

class TestNeighborSearch(unittest.TestCase):
    def test_5a11(self, use_populate=True):
        st = gemmi.read_pdb_string(FRAGMENT_5A11)
        a1 = st[0].sole_residue('A', gemmi.SeqId(37, ' '))[0]
        ns = gemmi.NeighborSearch(st[0], st.cell, 5)
        if use_populate:
            ns.populate()
        else:
            for n_ch, chain in enumerate(st[0]):
                for n_res, res in enumerate(chain):
                    for n_atom, atom in enumerate(res):
                        ns.add_atom(atom, n_ch, n_res, n_atom)
        marks = ns.find_atoms(a1.pos, a1.altloc, radius=3)
        m1, m2 = sorted(marks, key=lambda m: ns.dist(a1.pos, m.pos))
        self.assertAlmostEqual(ns.dist(a1.pos, m1.pos), 0, delta=5e-6)
        self.assertAlmostEqual(ns.dist(a1.pos, m2.pos), 0.13, delta=5e-3)
        cra2 = m2.to_cra(st[0])
        self.assertEqual(cra2.chain.name, 'B')
        self.assertEqual(str(cra2.residue.seqid), '37')
        self.assertEqual(cra2.atom.name, 'SG')
        marks2 = ns.find_neighbors(a1, 0.1, 3)
        self.assertEqual(len(marks2), 1)
        self.assertEqual(marks2[0], m2)

    def test_5a11_using_add_atom(self):
        self.test_5a11(use_populate=False)

    def test_1gtv(self):
        st = gemmi.read_pdb_string(FRAGMENT_1GTV)
        a1 = st[0].sole_residue('A', gemmi.SeqId(85, ' '))[0]
        ns = gemmi.NeighborSearch(st[0], st.cell, 5)
        ns.populate()
        marks = ns.find_atoms(a1.pos, a1.altloc, radius=3)
        self.assertEqual(len(marks), 2)
        for mark in marks:
            d = ns.dist(a1.pos, mark.pos)
            self.assertAlmostEqual(d, 0, delta=5e-6)
        marks2 = ns.find_neighbors(a1, 0.1, 3)
        self.assertEqual(len(marks2), 0)

        pos = gemmi.Position(100, 150, 200)
        mark = ns.find_nearest_atom(pos)
        cra = mark.to_cra(st[0])
        nim = st.cell.find_nearest_pbc_image(pos, mark.pos, 0)
        p = st.cell.find_nearest_pbc_position(cra.atom.pos, pos,
                                              mark.image_idx, inverse=True)
        self.assertAlmostEqual(nim.dist(), p.dist(cra.atom.pos))

    def test_b208(self):
        st = gemmi.read_structure(full_path('4oz7.pdb'))
        hoh208 = gemmi.Selection('B/208').copy_model_selection(st[0])
        self.assertEqual(hoh208.count_atom_sites(), 1)
        point = hoh208[0][0][0].pos
        for max_radius in [5, 20]:
            ns = gemmi.NeighborSearch(st[0], st.cell, max_radius).populate()
            marks = ns.find_atoms(point, min_dist=1e-6, radius=3)
            self.assertEqual(len(marks), 2)
            images = [m.image_idx for m in marks]
            self.assertEqual(sorted(images), [2, 3])
            image2 = st.cell.find_nearest_pbc_image(point, point, 2)
            self.assertEqual(image2.symmetry_code(), '3_545')
            self.assertAlmostEqual(image2.dist(), 2.9710496, delta=1e-6)
            image3 = st.cell.find_nearest_pbc_image(point, point, 3)
            self.assertEqual(image3.symmetry_code(), '4_355')


class TestContactSearch(unittest.TestCase):
    def test_radii_setting(self):
        cs = gemmi.ContactSearch(4.0)
        hg = gemmi.Element('Hg')
        self.assertEqual(cs.get_radius(hg), 0)
        cs.setup_atomic_radii(1, 0)
        cs.set_radius(hg, 1.5)
        self.assertEqual(cs.get_radius(hg), 1.5)

    def test_ignore_flag(self):
        st = gemmi.read_structure(full_path('4oz7.pdb'))
        st.setup_entities()
        ns = gemmi.NeighborSearch(st[0], st.cell, 5).populate()
        cs = gemmi.ContactSearch(4.0)
        cs.ignore = gemmi.ContactSearch.Ignore.SameResidue
        results = cs.find_contacts(ns)
        self.assertEqual(len(results), 607)
        cs.ignore = gemmi.ContactSearch.Ignore.SameChain
        results = cs.find_contacts(ns)
        self.assertEqual(len(results), 190)
        self.assertTrue(all(r.image_idx != 0
                            or r.partner1.chain is not r.partner2.chain
                            for r in results))


if __name__ == '__main__':
    unittest.main()
