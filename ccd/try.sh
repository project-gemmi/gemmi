#!/bin/bash -eu

#:rm -f gemmi/*.cif
../run-tests.sh G
echo START:

for x in "$@"; do
  code=${x^^}
  files+=("orig/$code.cif")
done
../build/gemmi drg --typeOut -v --output-dir=gemmi "${files[@]}" 2>drg.log

for x in "$@"; do
  code=${x^^}
   echo ========================$code========================
  [ -f orig/$code.cif ] || ./use.py $code
  ../build/gemmi mondiff ./acedrg/$code.cif ./gemmi/$code.cif ||:
done
echo END
