#!/bin/bash

# Index with dates helps limit the testing given time span:
#  ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx
# tail -n +3 entries.idx | cut -f1,3 | grep "05/16$" |\
# while read -r code date; do ./compare_pdb.sh $code; done

set -eu
cd `dirname $0`
code=${1,,}
tempd=/run/gemmi
pout="$tempd/$code-p"
gout="$tempd/$code-g"
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

ignore="\
^AUTHOR|\
^CISPEP|\
^COMPND|\
^CONECT|\
^DBREF |\
^EXPDTA|\
^FORMUL|\
^HELIX |\
^HET   |\
^HETNAM|\
^HETSYN|\
^JRNL  |\
^LINK  |\
^MASTER|\
^MTRIX|\
^MODRES|\
^ORIGX|\
^REMARK|\
^REVDAT|\
^SEQADV|\
^SEQRES|\
^SHEET |\
^SITE  |\
^SLTBRG|\
^SOURCE|\
^SPRSDE|\
^SSBOND"
#zgrep -v -f ignore "$pdb" > "$pout"
zgrep -v -E "$ignore" "$pdb" > "$pout"
../gemmi-convert --to=pdb "$cif" "$gout"
echo diff -u "$gout" "$pout"
diff -u "$gout" "$pout" | diffstat -q

# Add d or w as the second arg to show diff (using git diff for colors).
[[ ${2:-} = d ]] && git_diff_opt=
[[ ${2:-} = w ]] && git_diff_opt=--word-diff=color
if [[ "${git_diff_opt+set}" ]]; then
  echo git diff --no-index $git_diff_opt -- "$gout" "$pout"
  git diff --no-index $git_diff_opt -- "$gout" "$pout"
fi
