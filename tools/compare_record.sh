#!/bin/bash

set -eu
cd `dirname $0`
PDB_COPY="$PDB_DIR/structures/divided"
BIN=..
#BIN=../build
RECORD="^${1}"
# use grep "01/17$" to select month
tail -n +3 ../entries.idx | cut -f1,3 | grep "." | \
    while read -r code date; do
  code=${code,,}
  echo $code
  cif=${PDB_COPY:-/hdd}/mmCIF/${code:1:2}/${code}.cif.gz
  pdb=${PDB_COPY:-/hdd}/pdb/${code:1:2}/pdb${code}.ent.gz
  diff -U0 --label="$pdb" --label="from $cif" <(zgrep $RECORD "$pdb") \
    <($BIN/gemmi-convert --to=pdb "$cif" - | grep $RECORD) ||:
done
