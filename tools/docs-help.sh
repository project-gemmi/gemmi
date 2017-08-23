#!/bin/bash
set -eu

cd "$(dirname "$0")"/..

for prog in validate convert grep map mask; do
  echo "\$ gemmi-$prog -h" > docs/$prog-help.txt
  src/gemmi-$prog -h >> docs/$prog-help.txt
done

src/gemmi-convert tests/misc.cif tests/misc.json

git diff --stat docs/*-help.txt tests/misc.json
