#!/bin/bash
# Script that runs Refmac to generate intermediate files.
# Arguments: input.pdb output.crd [libin.cif]
# Output is concatenation of crd and rst files.
# Alternatively, it can run 0 cycles of Refmac starting from crd.
# Arguments: input.crd input.mtz output_basename

set -eu

#. $HOME/ccp4/ccp4-7.1/bin/ccp4.setup-sh
export CLIBD_MON=$HOME/ccp4/checkout/monomers/
export CCP4=$HOME/ccp4/install
export CLIB=$CCP4/lib
export CLIBD=$CCP4/lib/data
export CINCL=$CCP4/include
export CCP_SUPPRESS_HTML=1
REFMAC=$CCP4/bin/refmac5

if [[ "$2" =~ \.crd$ ]]; then
#### Preparing crd (pdb -> crd) ####

export CCP4_SCR=${TMPDIR-/tmp}/refmac_crd
mkdir -p $CCP4_SCR
rm -f $CCP4_SCR/refmac5_*

libstr=
if (( $# > 2 )); then
    libstr="libin $3"
fi

$REFMAC xyzin "$1" xyzout $CCP4_SCR/refmac5_.pdb $libstr << EOF
${REFMAC_EXTRA_RECORDS:-}
make exit yes
make form formatted
end
EOF
#make hydrogen no

cat $CCP4_SCR/refmac5_temp1.*_new.crd $CCP4_SCR/refmac5_temp1.*_new.rst > "$2"
echo "CRD file: $2"

elif [[ "$2" =~ \.mtz$ ]]; then
#### Using crd and mtz files in Refmac ####

OUT=$3
#labin FP=F SIGFP=SIGF FREE=FreeR_flag
$REFMAC xyzin $1 hklin $2 xyzout $OUT.pdb hklout $OUT.mtz << EOF
make cr prepared
refinement type restrained
ncycle 0
end
EOF
echo "Refmac refinement output file: $OUT.*"

else
    echo "Unexpected arguments"
fi
