#!/bin/bash
set -eu

cd "$(dirname "$0")"/..
BIN=.
[ -e build ] && BIN=build ||:

echo "\$ gemmi -h" > docs/gemmi-help.txt
$BIN/gemmi -h >> docs/gemmi-help.txt
for prog in bfit contact contents convert grep map mask sg validate; do
  echo "\$ gemmi $prog -h" > docs/$prog-help.txt
  $BIN/gemmi $prog -h >> docs/$prog-help.txt
done

$BIN/gemmi convert tests/misc.cif tests/misc.json

git diff --stat docs/*-help.txt tests/misc.json
