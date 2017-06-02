#!/bin/bash

# Index with dates helps limit the testing given time span:
#  ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx
# tail -n +3 entries.idx | cut -f1,3 | grep "05/16$" |\
# while read -r code date; do ./compare_pdb.sh $code; done

set -eu
cd `dirname $0`
code=${1,,}
tempd=/run/gemmi
pout="$tempd/$code-p.pdb"
gout="$tempd/$code-g.pdb"
cifout="$tempd/$code.cif"
cif=${PDB_COPY:-/hdd}/mmCIF/${code:1:2}/${code}.cif.gz
if [[ -d ${PDB_COPY:-/hdd}/pdb ]]; then
  pdb=${PDB_COPY:-/hdd}/pdb/${code:1:2}/pdb${code}.ent.gz
else
  pdb=$tempd/pdb${code}.ent.gz
  if [[ ! -e $pdb ]]; then
      remote=http://ftp.ebi.ac.uk/pub/databases/rcsb/pdb-remediated/data/structures/divided/pdb/${code:1:2}/pdb${code}.ent.gz
      curl $remote -o $pdb
  fi
fi

# TITLE, KEYWDS: line breaks can happen in different places
not_identical="\
^TITLE|\
^KEYWDS"
absent="\
^AUTHOR|\
^CISPEP|\
^CAVEAT|\
^COMPND|\
^CONECT|\
^DBREF|\
^FORMUL|\
^HELIX |\
^HET   |\
^HETNAM|\
^HETSYN|\
^JRNL  |\
^LINK  |\
^MASTER|\
^MDLTYP|\
^MODRES|\
^ORIGX|\
^REMARK|\
^REVDAT|\
^SEQADV|\
^SHEET |\
^SITE  |\
^SLTBRG|\
^SOURCE|\
^SPRSDE|\
^SSBOND"
zgrep -v -E "$not_identical|$absent" "$pdb" > "$pout"
inp="$cif"
[[ ${FROM_PDB:-} = 1 ]] && inp="$pout"
[[ ${PDB_BACK:-} = 1 ]] && ../gemmi-convert "$pout" "$cifout" && inp="$cifout"
../gemmi-convert --to=pdb "$inp" - | grep -v -E $not_identical > "$gout"
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
