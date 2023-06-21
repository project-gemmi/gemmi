#!/usr/bin/env python

import unittest
import math
import cmath
import gemmi
from common import full_path, numpy

# from 5nl9
FRAGMENT_WITH_UNK = """\
CRYST1   48.683   56.100   68.794  90.00 105.82  90.00 P 1 21 1      4
ATOM   2320  N   PRO B 145      27.193  12.782  46.813  1.00 43.92           N
ATOM   2321  CA  PRO B 145      27.971  11.557  46.488  1.00 44.82           C
ATOM   2322  C   PRO B 145      27.109  10.394  45.985  1.00 63.82           C
ATOM   2323  O   PRO B 145      25.961  10.269  46.451  1.00 67.23           O
ATOM   2324  CB  PRO B 145      28.599  11.170  47.841  1.00 46.06           C
ATOM   2325  CG  PRO B 145      28.536  12.388  48.681  1.00 49.45           C
ATOM   2326  CD  PRO B 145      27.276  13.086  48.254  1.00 44.74           C
ATOM   2327  OXT PRO B 145      27.597   9.587  45.164  1.00 83.13           O
TER    2328      PRO B 145
HETATM 2329 ZN    ZN A 201       3.818  16.117  25.579  1.00 33.14          ZN2+
HETATM 2330  K     K A 202      -6.950  11.545  40.774  0.80 41.06           K1+
HETATM 2331 ZN    ZN B 200       7.588   6.321  46.448  1.00 40.71          ZN2+
HETATM 2332  UNK UNX B 201       8.947   4.729  45.579  1.00 37.10           X
HETATM 2333  UNK UNX B 202       5.981   5.180  45.153  1.00 39.91           X
HETATM 2334  O   HOH A 301       6.823  17.093  -0.032  1.00 42.64           O
"""

class TestDensityCalculator(unittest.TestCase):
    def test_unk_does_not_cause_segfault(self):
        st = gemmi.read_pdb_string(FRAGMENT_WITH_UNK)
        for calculator in [gemmi.DensityCalculatorX,
                           gemmi.DensityCalculatorE,
                           gemmi.DensityCalculatorN]:
            dencalc = calculator()
            dencalc.d_min = 2.1
            dencalc.set_grid_cell_and_spacegroup(st)
            # we only check here that it doesn't crash
            dencalc.put_model_density_on_grid(st[0])

class TestExpandingToP1(unittest.TestCase):
    def test_1gdr(self):
        st = gemmi.read_pdb(full_path('pdb1gdr.ent'), max_line_length=72)
        sg = st.find_spacegroup()
        order = len(sg.operations())
        self.assertEqual(order, 12)
        sfcalc = gemmi.StructureFactorCalculatorX(st.cell)
        mtz = gemmi.Mtz(with_base=True)
        mtz.spacegroup = sg
        mtz.set_cell_for_all(st.cell)
        mtz.add_dataset('calc')
        mtz.add_column('FC', 'F')
        mtz.add_column('PHIC', 'P')
        hkl1 = (-7, 11, 5)
        hkl2 = (6, 8, -1)
        sf1 = sfcalc.calculate_sf_from_model(st[0], hkl1)
        sf2 = sfcalc.calculate_sf_from_model(st[0], hkl2)
        def polar(c):
            mag, phase = cmath.polar(c)
            return (mag, math.degrees(phase))
        if not numpy:
            return
        data = numpy.array([hkl1 + polar(sf1),
                            hkl2 + polar(sf2)], numpy.float32)
        mtz.set_data(data)
        mtz.expand_to_p1()
        st.transform_to_assembly('unit_cell',
                                 gemmi.HowToNameCopiedChain.AddNumber)
        # AsuData as_is is not yet in ASU.
        pre_asudata = mtz.get_f_phi('FC', 'PHIC', as_is=True)
        # Check structure factors after expand_to_p1(),
        # in particular we check here if phase shift is correct.
        for v in pre_asudata:
            usf = sfcalc.calculate_sf_from_model(st[0], v.hkl)
            self.assertAlmostEqual(v.value, usf / order, delta=1e-4)
        # The same after Mtz.ensure_asu(), which also changes phases.
        mtz.ensure_asu()
        asudata = mtz.get_f_phi('FC', 'PHIC', as_is=True)
        self.assertEqual(hkl1, tuple(asudata[0].hkl))
        for v in asudata:
            usf = sfcalc.calculate_sf_from_model(st[0], v.hkl)
            self.assertAlmostEqual(v.value, usf / order, delta=1e-4)
        # Finally, test AsuData.ensure_asu().
        pre_asudata.ensure_asu()
        for v in pre_asudata:
            usf = sfcalc.calculate_sf_from_model(st[0], v.hkl)
            self.assertAlmostEqual(v.value, usf / order, delta=1e-4)

if __name__ == '__main__':
    unittest.main()
