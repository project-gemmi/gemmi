#!/bin/bash
# Script that runs Refmac to generate intermediate files.
# Arguments: input.pdb output.crd [libin.cif]
# Output is concatenation of crd and rst files.

set -eu

#. $HOME/ccp4/ccp4-7.1/bin/ccp4.setup-sh
export CLIBD_MON=$HOME/ccp4/checkout/monomers/
export CCP4=$HOME/ccp4/install
export CLIB=$CCP4/lib
export CLIBD=$CCP4/lib/data
export CINCL=$CCP4/include
REFMAC=$CCP4/bin/refmac5

export CCP4_SCR=${TMPDIR-/tmp}/refmac_crd
mkdir -p $CCP4_SCR
rm -f $CCP4_SCR/refmac5_*

libstr=
if (( $# > 2 )); then
    libstr="libin $3"
fi

export CCP_SUPPRESS_HTML=1
$REFMAC xyzin "$1" xyzout $CCP4_SCR/refmac5_.pdb $libstr << EOF
make exit yes
make form formatted
end
EOF
#make hydrogen no

cat $CCP4_SCR/refmac5_temp1.*_new.crd $CCP4_SCR/refmac5_temp1.*_new.rst > "$2"
