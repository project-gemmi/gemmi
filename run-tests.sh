#!/bin/bash

# I use this script for building, testing and updating docs.

# Note: running Python scripts with locally built gemmi module requires setting
# PYTHONPATH. Locally, I'm using bash function "gpy" that starts python with
# gemmi imported ("gpy" for interactive mode, "gpy file.py" for script):
#   echo "import gemmi" > $HOME/.import-gemmi.py
#   gpy() { PYTHONPATH=$HOME/gemmi/gemmi/build/py PYTHONSTARTUP=$HOME/.import-gemmi.py python3 -q "$@"; }

# set up variables
set -eu
cd "$(dirname "$0")"
BUILD_DIR="$(pwd)"
[ -e build ] && BUILD_DIR="$(pwd)/build"
if [ -z "${PYTHON-}" ]; then
    PYTHON=`grep '^_\?Python_EXECUTABLE:' $BUILD_DIR/CMakeCache.txt | cut -d= -f2`
fi

# Build all, except when we called with an option to avoid full compilation:
#  G - only build the program,
#  P - only build Python bindings,
#  n - do not build, only run tests.
if [ $# = 1 ] && [ $1 = G ]; then
    (cd $BUILD_DIR && make -j4 gemmi_prog)
    exit
fi
if [ $# = 1 ] && [ $1 = P ]; then
    (cd $BUILD_DIR && make -j4 gemmi_py)
    exit
fi
if [ $# != 0 ] && [ $1 = n ]; then
    shift
else
    (cd $BUILD_DIR && make -j4 all check)
    ./tools/cmp-size.py build/gemmi build/libgemmi_cpp.so build/py/gemmi/gemmi_ext*
    ./tools/docs-help.sh
fi

# Run tests and checks.
if [ $# = 0 ] || [ $1 != i ]; then
    export PYTHONPATH=$BUILD_DIR/py
    export PATH="$BUILD_DIR:$PATH"
fi
$PYTHON -m unittest discover -s tests
grep :: $BUILD_DIR/py/gemmi/*.pyi ||:

if [ -z "${NO_DOCTEST-}" ]; then
    # 'make doctest' works only if sphinx-build was installed for python3.
    #(cd docs && make doctest SPHINXOPTS="-q -n")
    (cd docs && $PYTHON -m sphinx -M doctest . _build -q -n)
fi

echo "Building docs (in background)..."
./tools/header-list.py >docs/headers.rst
(cd docs && make -j4 html SPHINXOPTS="-q -n")&

# Usually, we stop here. Below are more extensive checks below that are run
# before making a release. They are run when this script is called with 'a'
# or with an option corresponding to the check.
[ $# = 0 ] && exit;

flake8 docs/ examples/ tests/ tools/

# this check is only for pip-installed packages
if [ $1 = i ]; then
    echo 'Check if gemmi package is found by setuptools pkg_resources...'
    $PYTHON -c '__requires__ = ["gemmi"]; import pkg_resources'
    echo 'OK'
fi

# a few quick extra (x) checks
if [ $1 = x -o $1 = a ]; then
    echo 'Run codespell'
    codespell include src prog python fortran tests examples docs wasm tools ||:

    echo 'Checking serialize.hpp...'
    $PYTHON tools/check_serialize.py
fi

if [ $1 = m -o $1 = a ]; then
    echo 'Creating, compiling and removing test_mmdb{1,2}.cpp'
    echo 'Example 1: gemmi -> mmdb'
    cmd="c++ -O -Wall -Wextra -pedantic -Wshadow -Iinclude test_mmdbX.cpp \
        -lmmdb2 -Lbuild -lgemmi_cpp -lz -Wl,-rpath=build -o test_mmdbX"
    awk '/Example 1/,/^}/' include/gemmi/mmdb.hpp > test_mmdb1.cpp
    ${cmd//X/1}
    echo "Converting tests/1orc.pdb to /tmp/example1.pdb"
    ./test_mmdb1 tests/1orc.pdb /tmp/example1.pdb
    echo 'Example 2: mmdb -> gemmi'
    awk '/Example 2/,/^}/' include/gemmi/mmdb.hpp > test_mmdb2.cpp
    ${cmd//X/2}
    echo "Converting tests/1orc.pdb to /tmp/example2.pdb"
    ./test_mmdb2 tests/1orc.pdb /tmp/example2.pdb
    rm -f test_mmdb1.cpp test_mmdb2.cpp test_mmdb1 test_mmdb2
fi

if [ $1 = h -o $1 = a ]; then
    echo "check if each header can be compiled on its own"
    # skip mmdb.hpp which requires also mmdb2 headers
    for f in include/gemmi/*.hpp; do
        if [ $f != include/gemmi/mmdb.hpp ]; then
            echo -n .
            gcc -c -fsyntax-only -Dsmall=.avoid/small $f
        fi
    done
    echo
fi

if [ $1 = v -o $1 = a ]; then
    echo "running tests under valgrind"
    valgrind $BUILD_DIR/cpptest
    PYTHONMALLOC=malloc valgrind $PYTHON -m unittest discover -s tests
    (cd docs && PYTHONMALLOC=malloc valgrind --trace-children=yes \
        $PYTHON $(which sphinx-build) -M doctest . _build -q -E)
fi

if [ $1 = p -o $1 = a ]; then
    echo "check if pdb->cif conversion produces valid CIF files"
    for i in $PDB_DIR/structures/divided/pdb/z[h-z]; do
        echo $i
        for j in $i/pdb*; do
            #echo $j
            gemmi convert $j /tmp/converted.cif
            gemmi validate /tmp/converted.cif
        done
        #echo -n .
    done
fi

if [ $1 = c -o $1 = a ]; then
    echo "test with tools/compare_pdb.sh"
    echo "=== part 1: cif -> pdb, pdb -> pdb ==="
    for pdbid in 2aa2 3ee3 4ii4 5oo5 6yy6 5xx5 6dd6 7uu7 7ww7 8ff8; do
        ./tools/compare_pdb.sh $pdbid  # cif -> pdb
        FROM_PDB=1 ./tools/compare_pdb.sh $pdbid  # pdb -> pdb
    done
    echo
    echo "=== part 2: cif -> cif -> pdb, pdb -> cif -> pdb ==="
    for pdbid in 2aa2 3ee3 4ii4 5oo5 6yy6 5xx5 6dd6 7ww7 8ff8; do
        VIA_CIF=1 ./tools/compare_pdb.sh $pdbid  # cif -> cif -> pdb
        FROM_PDB=1 VIA_CIF=1 ./tools/compare_pdb.sh $pdbid  # pdb -> cif -> pdb
    done
fi

if [ $1 = w -o $1 = a ]; then
    echo "check if wasm and project-gemmi/wasm can be built"
    [ -z ${EMSDK+x} ] && . $HOME/local/emsdk/emsdk_env.sh
    (cd wasm && make clean && make -j4 && npm run test)
fi

if [ $1 = M -o $1 = a ]; then
    echo "check bulk solvent mask against cctbx"
    for code in 1keb 1mru 5oo5 6dd6; do
        echo Testing mask for $code ...
        pdb="${PDB_DIR}/structures/divided/pdb/${code:1:2}/pdb${code}.ent.gz"
        zcat $pdb > /tmp/file.pdb
        mmtbx.python -m mmtbx.command_line.mask /tmp/file.pdb
        $PYTHON examples/maskcheck.py mask.ccp4 /tmp/file.pdb
        /bin/rm mask.ccp4
    done
fi
