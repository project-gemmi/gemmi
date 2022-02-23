#!/bin/bash

# Compares bulk solvent masks from Refmac and Gemmi.
# Requires PDB_DIR set to a local copy of the PDB archive,
# with directories pdb and structure_factors in $PDB_DIR/structures/divided.
# Example usage: ./compare_mask.sh 1keb

set -eu
REFMAC=${REFMAC-refmac5}
LOCAL_COPY="$PDB_DIR/structures/divided"
TEMPDIR=${TMPDIR-/tmp}/gemmi

cd `dirname $0`/..
BUILD_DIR="$(pwd)"
[ -e build ] && BUILD_DIR="$(pwd)/build"
PYTHON=`grep PYTHON_EXECUTABLE:FILEPATH= $BUILD_DIR/CMakeCache.txt | cut -d= -f2`
export PYTHONPATH=$BUILD_DIR
export PATH="$BUILD_DIR:$PATH"

code=${1,,}
pdb_tail="pdb/${code:1:2}/pdb${code}.ent.gz"
uncompressed_pdb="$TEMPDIR/$code.pdb"
mkdir -p "$TEMPDIR"
echo "create $uncompressed_pdb"
zcat "$LOCAL_COPY/$pdb_tail" > "$uncompressed_pdb"

mtz="$TEMPDIR/$code.mtz"
echo "create $mtz"
sf_tail="structure_factors/${code:1:2}/r${code}sf.ent.gz"
gemmi cif2mtz $LOCAL_COPY/$sf_tail $mtz

refmac_mask=$TEMPDIR/out.msk

extra_refmac_options=
#extra_refmac_options="solvent process islands noremove"
export CCP_SUPPRESS_HTML=1
refmac_log=$TEMPDIR/out.log
echo "running $REFMAC, output redirected to $refmac_log"
$REFMAC XYZIN $uncompressed_pdb HKLIN $mtz XYZOUT $TEMPDIR/out.pdb \
        HKLOUT $TEMPDIR/out.mtz MSKOUT $refmac_mask >$refmac_log <<eof
MAKE HYDRogens No
REFI TYPE RIGID
RIGIdbody NCYC 0
$extra_refmac_options
eof

echo "Refmac mask should be in $refmac_mask, running:"
echo "examples/maskcheck.py $refmac_mask $uncompressed_pdb"
export REFMAC
$PYTHON examples/maskcheck.py $refmac_mask $uncompressed_pdb
