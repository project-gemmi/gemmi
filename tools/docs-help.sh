#!/bin/bash
set -eu

cd "$(dirname "$0")"/..
BIN=.
[ -e build ] && BIN=build ||:

echo "\$ gemmi -h" > docs/gemmi-help.txt
$BIN/gemmi -h >> docs/gemmi-help.txt
for fn in prog/*.cpp; do
  prog=$(basename $fn .cpp)
  if grep -q ${prog}_main prog/main.cpp; then
    echo "\$ gemmi $prog -h" > docs/$prog-help.txt
    $BIN/gemmi $prog -h >> docs/$prog-help.txt
  fi
done

echo "\$ gemmi cif2mtz --print-spec" > docs/cif2mtz-spec.txt
$BIN/gemmi cif2mtz --print-spec >> docs/cif2mtz-spec.txt
echo "" >> docs/cif2mtz-spec.txt
echo "\$ gemmi cif2mtz --print-spec --unmerged" >> docs/cif2mtz-spec.txt
$BIN/gemmi cif2mtz --print-spec --unmerged >> docs/cif2mtz-spec.txt

$BIN/gemmi cif2json tests/misc.cif tests/misc.json

git diff --stat docs/*-help.txt tests/misc.json
