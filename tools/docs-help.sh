#!/bin/bash
set -eu

cd "$(dirname "$0")"/..
BIN=.
[ -e build ] && BIN=build ||:

echo "\$ gemmi -h" > docs/gemmi-help.txt
$BIN/gemmi -h >> docs/gemmi-help.txt
for prog in align blobs cif2json cif2mtz contact contents convert \
    fprime grep h json2cif map map2sf mask merge mondiff mtz mtz2cif \
    reindex residues rmsz sf2map sfcalc sg tags validate wcn; do
  echo "\$ gemmi $prog -h" > docs/$prog-help.txt
  $BIN/gemmi $prog -h >> docs/$prog-help.txt
done

$BIN/gemmi cif2json tests/misc.cif tests/misc.json

git diff --stat docs/*-help.txt tests/misc.json
