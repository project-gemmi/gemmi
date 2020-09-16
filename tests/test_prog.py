#!/usr/bin/env python

import os
import sys
import subprocess
import unittest

TOP_DIR = os.path.join(os.path.dirname(__file__), "..")

def has_gemmi():
    try:
        subprocess.check_output(['gemmi', '--version'], cwd=TOP_DIR)
    except OSError:
        return False
    return True

@unittest.skipIf(not has_gemmi(), "Program gemmi not found.")
class TestProg(unittest.TestCase):
    def do(self, example):
        cmd, _, rest = example.partition('\n')
        assert cmd.startswith('$ gemmi ')
        output = subprocess.check_output(cmd[2:], shell=True, cwd=TOP_DIR,
                                         stderr=subprocess.STDOUT)
        expected_lines = rest.splitlines()
        output_lines = output.decode().splitlines()
        if expected_lines[0].strip() == '[...]':
            expected_lines.pop(0)
            output_lines = output_lines[-len(expected_lines):]
        self.assertEqual(expected_lines, output_lines)

    def test_fprime1(self):
        # example from utils.rst
        self.do('''\
$ gemmi fprime --wavelength=1.2 Se
Element	 E[eV]	Wavelength[A]	   f'   	  f"
Se	10332.0	 1.2    	 -1.4186	0.72389
''')

    def test_fprime2(self):
        self.do('''\
$ gemmi fprime --energy=12345 Se
Element	 E[eV]	Wavelength[A]	   f'   	  f"
Se	12345.0	 1.00433	 -3.1985	0.52258
''')

    def test_sfcalc1(self):
        self.do('''\
$ gemmi sfcalc --compare=tests/2242624.hkl tests/2242624.cif
RMSE=0.019256  0.2252%  max|dF|=0.04784  R=0.191%  sum(F^2)_ratio=1.00094
''')

    def test_sfcalc2(self):
        self.do('''\
$ gemmi sfcalc --wavelength=0 --compare=tests/2242624.hkl tests/2242624.cif
RMSE=0.10942  1.295%  max|dF|=0.1498  R=1.279%  sum(F^2)_ratio=1.01019
''')

    def test_sfcalc3(self):
        self.do('''\
$ gemmi sfcalc --ciffp --compare=tests/2242624.hkl tests/2242624.cif
RMSE=0.019724  0.2307%  max|dF|=0.04863  R=0.196%  sum(F^2)_ratio=1.00101
''')

    def test_sfcalc4(self):
        self.do('''\
$ gemmi sfcalc -w0 --hkl=4,9,0 --hkl=5,6,4 --hkl=1,1,1 tests/2013551.cif
 (4 9 0)	1.12314668	180.000000
 (5 6 4)	1.49089619	180.000000
 (1 1 1)	11.22159034	0.000000
''')

    @unittest.skipIf(sys.platform == 'win32', 'with MSVC it differs slightly')
    def test_sfcalc_5wkd(self):
        self.do('''\
$ gemmi sfcalc --blur=12 --dmin=2.5 --rate=2.5 --rcut=1e-7 --test -v tests/5wkd.pdb
[...]
RMSE=3.1854e-05  4.493e-05%  max|dF|=0.0001630  R=0.000%  <dPhi>=4.580e-06
''')  # noqa: E501

    @unittest.skipIf(sys.platform == 'win32', 'with MSVC it differs slightly')
    def test_sfcalc_1pfe(self):
        self.do('''\
$ gemmi sfcalc --dmin=9 --rate=4 --blur=60 --rcut=1e-7 --test -v tests/1pfe.cif.gz
[...]
RMSE=0.0011646  0.0001833%  max|dF|=0.004507  R=0.000%  <dPhi>=7.624e-06
''')  # noqa: E501

    # example from utils.rst
    def test_align_text(self):
        self.do('''\
$ gemmi align -p --match=0 --gapo=0 --text-align Saturday Sunday
Score: -3   CIGAR: 1M2I5M
Saturday
|  |.|||
S--unday
''')


if __name__ == '__main__':
    unittest.main()
