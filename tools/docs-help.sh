#!/bin/bash
set -eu

cd "$(dirname "$0")"/..
BIN=.
[ -e build ] && BIN=build ||:

echo "\$ gemmi -h" > docs/gemmi-help.txt
$BIN/gemmi -h >> docs/gemmi-help.txt
for prog in align blobs cif2json cif2mtz contact contents convert ecalc \
    fprime grep h json2cif map map2sf mask merge mondiff mtz mtz2cif \
    prep reindex residues rmsz sf2map sfcalc sg tags validate wcn xds2mtz; do
  echo "\$ gemmi $prog -h" > docs/$prog-help.txt
  $BIN/gemmi $prog -h >> docs/$prog-help.txt
done

echo "\$ gemmi cif2mtz --print-spec" > docs/cif2mtz-spec.txt
$BIN/gemmi cif2mtz --print-spec >> docs/cif2mtz-spec.txt
echo "" >> docs/cif2mtz-spec.txt
echo "\$ gemmi cif2mtz --print-spec --unmerged" >> docs/cif2mtz-spec.txt
$BIN/gemmi cif2mtz --print-spec --unmerged >> docs/cif2mtz-spec.txt

$BIN/gemmi cif2json tests/misc.cif tests/misc.json

git diff --stat docs/*-help.txt tests/misc.json
