#!/usr/bin/env python

import unittest
from gemmi import Element, IT92_get_exact


class TestElem(unittest.TestCase):
    def test_name(self):
        for name in ELEMENT_MASS:
            el = Element(name)
            self.assertEqual(el.name, name)
        self.assertEqual(Element('AL').name, 'Al')
        self.assertEqual(Element('al').name, 'Al')
        # We check only the first two characters now.
        self.assertEqual(Element('alt').name, 'Al')
        self.assertEqual(Element('Q').name, 'X')
        self.assertEqual(Element('QQ').name, 'X')
        self.assertEqual(Element('--').name, 'X')
        self.assertEqual(Element('').name, 'X')

    def test_weight(self):
        for name, mass in ELEMENT_MASS.items():
            self.assertAlmostEqual(Element(name).weight, mass, delta=1e-3)

    def test_atomic_number(self):
        self.assertEqual(Element('O').atomic_number, 8)
        self.assertEqual(Element('D').atomic_number, 1)
        self.assertEqual(Element('X').atomic_number, 0)
        for n in range(119):
            self.assertEqual(Element(n).atomic_number, n)

    def test_properties(self):
        self.assertEqual(Element('K').vdw_r, 2.75)

    def test_coef(self):
        h = Element('H')
        d = Element('D')
        fe = Element('Fe')
        cf = Element('Cf')
        es = Element('Es')
        self.assertIsNotNone(cf.it92)
        self.assertIsNotNone(cf.c4322)
        self.assertIsNone(es.it92)
        self.assertIsNone(es.c4322)
        self.assertEqual(h.it92.get_coefs(), d.it92.get_coefs())
        self.assertEqual(h.c4322.get_coefs(), d.c4322.get_coefs())
        self.assertEqual(Element('X').it92.get_coefs(),
                         Element('O').it92.get_coefs())
        fe_coefs = fe.it92.get_coefs()
        fe0_coefs = IT92_get_exact(fe, 0).get_coefs()
        fe2_coefs = IT92_get_exact(fe, 2).get_coefs()
        self.assertIsNone(IT92_get_exact(fe, 4))
        self.assertEqual(fe_coefs, fe0_coefs)
        self.assertNotEqual(fe_coefs, fe2_coefs)


ELEMENT_MASS = {
    'H': 1.00794, 'He': 4.0026, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811,
    'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'F': 18.998403, 'Ne': 20.1797,
    'Na': 22.98977, 'Mg': 24.305, 'Al': 26.981539, 'Si': 28.0855,
    'P': 30.973761, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,
    'K': 39.0983, 'Ca': 40.078, 'Sc': 44.95591, 'Ti': 47.867, 'V': 50.9415,
    'Cr': 51.9961, 'Mn': 54.93805, 'Fe': 55.845, 'Co': 58.9332,
    'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.723, 'Ge': 72.64,
    'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798,
    'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585, 'Zr': 91.224, 'Nb': 92.9064,
    'Mo': 95.95, 'Tc': 98, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
    'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.76,
    'Te': 127.6, 'I': 126.90447, 'Xe': 131.293, 'Cs': 132.905, 'Ba': 137.327,
    'La': 138.905, 'Ce': 140.116, 'Pr': 140.908, 'Nd': 144.24, 'Pm': 145,
    'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.925, 'Dy': 162.5,
    'Ho': 164.93, 'Er': 167.259, 'Tm': 168.934, 'Yb': 173.05, 'Lu': 174.967,
    'Hf': 178.49, 'Ta': 180.948, 'W': 183.84, 'Re': 186.207, 'Os': 190.23,
    'Ir': 192.217, 'Pt': 195.084, 'Au': 196.967, 'Hg': 200.59, 'Tl': 204.383,
    'Pb': 207.2, 'Bi': 208.98, 'Po': 209, 'At': 210, 'Rn': 222,
    'Fr': 223, 'Ra': 226, 'Ac': 227, 'Th': 232.038, 'Pa': 231.036,
    'U': 238.029, 'Np': 237, 'Pu': 244, 'Am': 243, 'Cm': 247,
    'Bk': 247, 'Cf': 251, 'Es': 252, 'Fm': 257, 'Md': 258,
    'No': 259, 'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271,
    'Bh': 272, 'Hs': 270, 'Mt': 276,
    'D': 2.0141, 'X': 0,
}

if __name__ == '__main__':
    unittest.main()
