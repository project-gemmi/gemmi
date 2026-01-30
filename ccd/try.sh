#!/bin/bash -eu

#:rm -f gemmi/*.cif
../run-tests.sh G
echo START:

for x in "$@"; do
  code=${x^^}
  files+=("orig/$code.cif")
done
../build/gemmi drg --typeOut -v --output-dir=gemmi "${files[@]}" 2>drg.log

echo ====================  COMPARING acedrg_type ===========
for x in "$@"; do
  code=${x^^}
  echo ====================  $code ===========
  (  ../build/gemmi mondiff acedrg/${code}_atomTypes.cif gemmi/$code.cif | grep acedrg_type) ||:
done
echo ====================  COMPARING restraints ===========
for x in "$@"; do
  code=${x^^}
   echo ========================$code========================
  [ -f orig/$code.cif ] || ./use.py $code
  echo ====================  $code ===========
  #echo ./build/gemmi mondiff ./ccd/acedrg/$code.cif ./ccd/gemmi/$code.cif
  ../build/gemmi mondiff ./acedrg/$code.cif ./gemmi/$code.cif ||:
done
echo END
