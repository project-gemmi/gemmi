#!/usr/bin/env python3
"""Generate chemical-adjustment SVG examples with OpenBabel.

This script writes the SVG files used in docs/chemistry.rst:
  docs/img/adj_<rule>_before.svg
  docs/img/adj_<rule>_after.svg

Requirements:
  - obabel command available in PATH

Usage:
  python docs/gen_adj_svgs.py
  python docs/gen_adj_svgs.py --out-dir docs/img
"""

import argparse
import os
import shutil
import subprocess
import sys

PAIRS = [
    ("oxoacid_phosphate", "OP(=O)(O)*", "OP(=O)([O-])*"),
    ("oxoacid_sulfate", "O=S(=O)(O)O", "O=S(=O)([O-])O"),
    ("oxoacid_sulfite", "OS(=O)O", "O=S(=O)[O-]"),
    ("nitro_group", "[N+](=O)(O)C", "[N+](=O)([O-])C"),
    ("single_bond_oxide", "C[N+](C)(C)O", "C[N+](C)(C)[O-]"),
    ("hexafluorophosphate", "F[P](F)(F)(F)(F)F", "[H][P](F)(F)(F)(F)(F)F"),
    ("carboxy_asp", "O=C(O)C(*)", "O=C([O-])C(*)"),
    ("terminal_carboxylate", "O=C(O)CN(*)", "O=C([O-])CN(*)"),
    ("guanidinium", "CNC(=N)N", "CNC(=[NH2+])N"),
    ("amino_ter_amine", "NCC(=O)N(*)", "[NH3+]CC(=O)N(*)"),
    ("terminal_amine", "CCCCCCN", "CCCCCC[NH3+]"),
    ("protonated_amide_n", "CCC(=O)N", "CCC(=O)[NH2+]"),
    ("add_n_terminal_h3", "[*]C([NH2+])C(=O)[O-]", "[*]C([NH3+])C(=O)[O-]"),
]


def run_obabel(smiles, out_path):
    cmd = ["obabel", f"-:{smiles}", "-O", out_path]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        msg = proc.stderr.strip() or proc.stdout.strip() or "unknown OpenBabel error"
        raise RuntimeError(f"OpenBabel failed for {smiles!r}: {msg}")


def main():
    parser = argparse.ArgumentParser(description="Generate adjustment SVGs with OpenBabel")
    parser.add_argument("--out-dir", default="docs/img", help="output directory for SVG files")
    args = parser.parse_args()

    if shutil.which("obabel") is None:
        print("error: obabel not found in PATH", file=sys.stderr)
        return 2

    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    generated = 0
    for rule, before_smiles, after_smiles in PAIRS:
        before_path = os.path.join(out_dir, f"adj_{rule}_before.svg")
        after_path = os.path.join(out_dir, f"adj_{rule}_after.svg")
        run_obabel(before_smiles, before_path)
        run_obabel(after_smiles, after_path)
        generated += 2

    print(f"generated {generated} SVG files in {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
