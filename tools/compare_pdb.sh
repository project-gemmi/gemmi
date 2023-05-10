#!/bin/bash

# This script converts and compares PDB entries, preferably using local copy
# of both mmCIF and PDB files from the PDB.
# By default it converts mmCIF to PDB and compares if the lines are identical
# (apart from the not_identical and absent fields listed below).
#
# ./compare_pdb.sh 5moo    # entry 5moo, shows diffstat (if any)
# ./compare_pdb.sh 5moo d  # shows "diff -u".
# ./compare_pdb.sh 5moo g  # similar, but uses "git diff" for nicer output
# ./compare_pdb.sh 5moo w  # similar, but with option --word-diff=color
#
# Two env. variables can be used for different conversion route:
#
# ./compare_pdb.sh 5moo                       # converts cif -> pdb
# FROM_PDB=1 ./compare_pdb.sh 5moo            # converts pdb -> pdb
# VIA_CIF=1  ./compare_pdb.sh 5moo            # converts cif -> cif -> pdb
# FROM_PDB=1 VIA_CIF=1 ./compare_pdb.sh 5moo  # converts pdb -> cif -> pdb
#
# The final comparison is always between the PDB file from conversion
# and the original PDB with selected records removed.
#
# Index with dates helps limit the testing given time span:
# wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx
# tail -n +3 entries.idx | cut -f1,3 | grep "05/../17$" |\
# while read -r code date; do ./compare_pdb.sh $code; done

set -eu
cd `dirname $0`
LOCAL_COPY="$PDB_DIR/structures/divided"
REMOTE_COPY="http://ftp.ebi.ac.uk/pub/databases/rcsb/pdb-remediated/data/structures/divided"
TEMPDIR=${TMPDIR-/tmp}/gemmi
BIN=..
[ -e ../build ] && BIN=../build ||:
code=${1,,}
mkdir -p "$TEMPDIR"
pout="$TEMPDIR/$code-p.pdb"
gout="$TEMPDIR/$code-g.pdb"
cifout="$TEMPDIR/$code.cif"  # used with VIA_CIF=1

cif_tail="mmCIF/${code:1:2}/${code}.cif.gz"
if [[ -d ${LOCAL_COPY}/mmCIF ]]; then
  cif=$LOCAL_COPY/$cif_tail
else
  cif=$TEMPDIR/${code}.cif.gz
  if [[ ! -e $cif ]]; then
    curl $REMOTE_COPY/$cif_tail -o $cif
  fi
fi

pdb_tail="pdb/${code:1:2}/pdb${code}.ent.gz"
if [[ -d $LOCAL_COPY/pdb ]]; then
  pdb=$LOCAL_COPY/$pdb_tail
else
  pdb=$TEMPDIR/pdb${code}.ent.gz
  if [[ ! -e $pdb ]]; then
    curl $REMOTE_COPY/$pdb_tail -o $pdb
  fi
fi

# TITLE, KEYWDS: line breaks can happen in different places
# In HELIX helixID may differ.

not_identical="\
^TITLE|\
^KEYWDS|\
^SCALE|\
^HELIX "
absent="\
^AUTHOR|\
^CAVEAT|\
^COMPND|\
^CONECT|\
^FORMUL|\
^HET   |\
^HETNAM|\
^HETSYN|\
^JRNL  |\
^MASTER|\
^MDLTYP|\
^ORIGX|\
^REMARK   1|\
^REMARK   3|\
^REMARK   4|\
^REMARK 100|\
^REMARK 2..|\
^REMARK 300|\
^REMARK 375|\
^REMARK 465|\
^REMARK 470|\
^REMARK 500|\
^REMARK 525|\
^REMARK 610|\
^REMARK 620|\
^REMARK 630|\
^REMARK 800|\
^REMARK 900|\
^REMARK 999|\
^REVDAT|\
^SEQADV|\
^SITE  |\
^SLTBRG|\
^SOURCE|\
^SPRSDE"
zgrep -v -E "$not_identical|$absent" "$pdb" > "$pout"
inp="$cif"
[[ ${FROM_PDB:-} = 1 ]] && inp="$pout"
if [[ ${VIA_CIF:-} = 1 ]]; then
    echo "$(basename "$inp") -> $(basename "$cifout")"
    $BIN/gemmi convert "$inp" "$cifout"
    inp="$cifout"
fi
$BIN/gemmi convert --to=pdb "$inp" - | grep -v -E $not_identical > "$gout"
echo "Comparing ($(basename "$inp") ->) $gout vs $pout"

# Add d or w as the second arg to show diff (using git diff for colors).
if [[ ${2:-} = d ]]; then
    diff -u "$gout" "$pout"
elif [[ ${2:-} = g ]]; then
    git diff --no-index -- "$gout" "$pout"
elif [[ ${2:-} = w ]]; then
    git diff --no-index --word-diff=color -- "$gout" "$pout"
elif [[ ${2:-} = n ]]; then
    numdiff -V -r 1e-4 "$gout" "$pout"
else
    diff -u "$gout" "$pout" | diffstat -q
fi

[[ ${3:-} = c ]] && /bin/rm "$gout" "$pout"
