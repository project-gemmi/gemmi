#!/bin/sh

# I use this script for building, testing and updating docs.

# Note: running Python scripts with locally built gemmi module requires setting
# PYTHONPATH. Locally, I'm using bash function "gpy" that starts python with
# gemmi imported ("gpy" for interactive mode, "gpy file.py" for script):
#   echo "import gemmi" > $HOME/.import-gemmi.py
#   gpy() { PYTHONPATH=$HOME/gemmi/gemmi/build PYTHONSTARTUP=$HOME/.import-gemmi.py python3 -q "$@"; }

# set up variables
set -eu
cd "$(dirname "$0")"
BUILD_DIR="$(pwd)"
[ -e build ] && BUILD_DIR="$(pwd)/build"
if [ -z "${PYTHON-}" ]; then
    PYTHON=`grep ^_Python_EXECUTABLE: $BUILD_DIR/CMakeCache.txt | cut -d= -f2`
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
    ./tools/cmp-size.py build/gemmi build/gemmi.*.so
    ./tools/docs-help.sh
fi
(cd docs && make -j4 html SPHINXOPTS="-q -n")

# Run tests and checks.
if [ $# = 0 ] || [ $1 != i ]; then
    export PYTHONPATH=$BUILD_DIR
    export PATH="$BUILD_DIR:$PATH"
fi
$PYTHON -m unittest discover -s tests
./tools/header-list.py >docs/headers.rst

if [ -z "${NO_DOCTEST-}" ]; then
    # 'make doctest' works only if sphinx-build was installed for python3.
    (cd docs && make doctest SPHINXOPTS="-q -n -E")
fi

flake8 docs/ examples/ tests/ tools/ setup.py
# We want to avoid having '::' in pydoc from pybind11.
$PYTHON -m pydoc gemmi | grep :: ||:


# Usually, we stop here. Below are more extensive checks below that are run
# before making a release. They are run when this script is called with 'a'
# or with an option corresponding to the check.
[ $# = 0 ] && exit;

if [ $1 = i ]; then
    echo 'Check if gemmi package is found by setuptools pkg_resources...'
    $PYTHON -c '__requires__ = ["gemmi"]; import pkg_resources'
    echo 'OK'
fi

if [ $1 = m -o $1 = a ]; then
    echo 'Creating, compiling and removing test_mmdb.cpp'
    echo 'Example 1'
    cmd="c++ -O -Wall -Wextra -pedantic -Wshadow -Iinclude test_mmdb.cpp \
        -lmmdb2 -Lbuild -lgemmi_cpp -lz -o test_mmdb"
    awk '/Example 1/,/^}/' include/gemmi/mmdb.hpp > test_mmdb.cpp
    $cmd
    echo 'Example 2'
    awk '/Example 2/,/^}/' include/gemmi/mmdb.hpp > test_mmdb.cpp
    $cmd
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
    echo "check if wasm/mtz and project-gemmi/wasm can be built"
    [ -z ${EMSDK+x} ] && . $HOME/local/emsdk/emsdk_env.sh
    (cd wasm/mtz && ./compile.sh)
    (cd ../wasm/convert && make clean && make)
fi
