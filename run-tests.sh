#!/bin/sh

# I use this script for building, testing and updating docs.
# Running Python scripts with locally built gemmi module requires setting
# PYTHONPATH. I'm using the following alias (bash function actually) to set it:
# echo "import gemmi" > $HOME/.import-gemmi.py
# gpy() { PYTHONPATH=$HOME/gemmi/gemmi/build PYTHONSTARTUP=$HOME/.import-gemmi.py python3 -q "$@"; }

set -eu
cd "$(dirname "$0")"
BUILD_DIR="$(pwd)"
[ -e build ] && BUILD_DIR="$(pwd)/build"
PYTHON=`grep PYTHON_EXECUTABLE:FILEPATH= $BUILD_DIR/CMakeCache.txt | cut -d= -f2`

if [ $# = 1 ] && [ $1 = G ]; then
    (cd $BUILD_DIR && make -j4 program)
    exit
fi
if [ $# = 1 ] && [ $1 = P ]; then
    (cd $BUILD_DIR && make -j4 py)
    exit
fi
if [ $# = 0 ] || [ $1 != n ]; then
    (cd $BUILD_DIR && make -j4 all check)
    ./tools/cmp-size.py build/gemmi build/gemmi.*.so
    ./tools/docs-help.sh
fi
(cd docs && make -j4 html SPHINXOPTS="-q -n")
export PYTHONPATH=$BUILD_DIR
export PATH="$BUILD_DIR:$PATH"
$PYTHON -m unittest discover -s tests
./tools/header-list.py >docs/headers.rst
# Where python 2 and 3 output differ, the docs have output from v3.
# So 'make doctest' works only if sphinx-build was installed for python3.
(cd docs && make doctest SPHINXOPTS="-q -n -E")
flake8 docs/ examples/ tests/ tools/ setup.py
$PYTHON -m pydoc gemmi | grep :: ||:

[ $# = 0 ] && exit;

if [ $1 = m -o $1 = a ]; then
    echo 'Creating, compiling and removing test_mmdb.cpp'
    echo 'Example 1'
    awk '/Example 1/,/^}/' include/gemmi/mmdb.hpp > test_mmdb.cpp
    c++ -O -Wall -Wextra -pedantic -Wshadow -Iinclude test_mmdb.cpp -lmmdb2 -o test_mmdb
    echo 'Example 2'
    awk '/Example 2/,/^}/' include/gemmi/mmdb.hpp > test_mmdb.cpp
    c++ -O -Wall -Wextra -pedantic -Wshadow -Iinclude test_mmdb.cpp -lmmdb2 -o test_mmdb
    rm -f test_mmdb.cpp
fi

if [ $1 = h -o $1 = a ]; then
    echo "check if each header can be compiled on its own"
    # skip mmdb.hpp which requires also mmdb2 headers
    for f in include/gemmi/*.hpp; do
        if [ $f != include/gemmi/mmdb.hpp ]; then
            echo -n .
            gcc-9 -c -fsyntax-only $f
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

if [ $1 = w -o $1 = a ]; then
    [ -z ${EMSDK+x} ] && . $HOME/local/emsdk/emsdk_env.sh
    (cd wasm/mtz && ./compile.sh)
    (cd ../wasm/convert && make clean && make)
fi
