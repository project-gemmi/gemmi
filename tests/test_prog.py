#!/usr/bin/env python

import os
import sys
import subprocess
import unittest
import gemmi

TOP_DIR = os.path.join(os.path.dirname(__file__), "..")

# Skip tests if the program is not installed,
# or the installed version doesn't match the library.
def has_gemmi():
    try:
        v = subprocess.check_output(['gemmi', '--version'], cwd=TOP_DIR)
    except OSError:  # usually FileNotFoundError
        return False
    except subprocess.CalledProcessError as e:
        print('Error when running gemmi --version:\n  cmd:', e.cmd,
              '\n  output:', e.output.decode(),
              '\n  returncode:', e.returncode)
        return False
    return v.split()[1].decode() == gemmi.__version__

@unittest.skipIf(not has_gemmi(), "Program gemmi not found.")
class TestProg(unittest.TestCase):
    def do(self, example):
        cmd, _, rest = example.partition('\n')
        assert cmd.startswith('$ gemmi ')
        output = subprocess.check_output(cmd[2:], shell=True, cwd=TOP_DIR,
                                         stderr=subprocess.STDOUT)
        expected_lines = rest.splitlines()
        output_lines = output.decode().splitlines()
        if expected_lines[0].startswith('[...]'):
            t = expected_lines.pop(0)
            if len(t) > 8:
                n1, n2 = t[5:].split()
                output_lines = output_lines[int(n1):int(n2)]
            else:
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
 (1 1 1)	11.22159034	360.000000
''')

    @unittest.skipIf(sys.platform == 'win32', 'with MSVC it differs slightly')
    def test_sfcalc_5wkd(self):
        self.do('''\
$ gemmi sfcalc --blur=12 --dmin=2.5 --rate=2.5 --rcut=1e-7 --test -v tests/5wkd.pdb
[...]
RMSE=3.3342e-05  4.703e-05%  max|dF|=0.0001706  R=0.000%  <dPhi>=4.718e-06
''')  # noqa: E501

    @unittest.skipIf(sys.platform == 'win32', 'with MSVC it differs slightly')
    def test_sfcalc_1pfe(self):
        self.do('''\
$ gemmi sfcalc --dmin=9 --rate=4 --blur=60 --rcut=1e-7 --test -v tests/1pfe.cif.gz
[...]
RMSE=0.0010414  0.0001639%  max|dF|=0.003530  R=0.000%  <dPhi>=8.067e-06
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

    def test_cif2mtz_5e5z(self):
        self.do('''\
$ gemmi mtz2cif tests/5e5z.mtz -
[...] -158 -155
1 2 0 46.859 0.6285 68.443 0.4591 o
1 2 1 ? ? ? ? x
1 2 2 38.249 0.3855 61.841 0.3117 f
''')

    def test_contacts_4oz7(self):
        self.do('''\
$ gemmi contact --maxdist=2.5 --ignore=4 tests/4oz7.pdb
metalc5      N   22Q A   1                CU   CU1 B 101     1555   6345  2.05
metalc6      S   22Q A   1                CU   CU1 B 101     1555   6345  2.26
metalc7      N   22Q B   1                CU   CU1 A 101     1555   6344  2.07
metalc8      S   22Q B   1                CU   CU1 A 101     1555   6344  2.22
''')

if __name__ == '__main__':
    unittest.main()
