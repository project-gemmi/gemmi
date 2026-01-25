#!/bin/bash -eu

#:rm -f gemmi/*.cif
../run-tests.sh G

for x in "$@"; do
  code=${x^^}
  files+=("orig/$code.cif")
done
../build/gemmi drg -v --output-dir=gemmi "${files[@]}" 2>drg.log

for x in "$@"; do
  code=${x^^}
  # echo ========================$code========================
  [ -f orig/$code.cif ] || ./use.py $code
  echo ====================  $code ===========
  #echo ./build/gemmi mondiff ./ccd/acedrg/$code.cif ./ccd/gemmi/$code.cif
  ../build/gemmi mondiff ./acedrg/$code.cif ./gemmi/$code.cif 
done
