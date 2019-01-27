#!/bin/sh

set -eu
cd "$(dirname "$0")"
BUILD_DIR="$(pwd)"
[ -e build ] && BUILD_DIR="$(pwd)/build"
(cd $BUILD_DIR && make -j4 all gemmi-crdrst check)
./tools/docs-help.sh
(cd docs && make -j4 html SPHINXOPTS="-q -n")
export PYTHONPATH=$BUILD_DIR
python3 -m unittest discover -s tests
# Where python 2 and 3 output differ, the docs have output from v3.
# So 'make doctest' works only if sphinx-build was installed for python3.
(cd docs && make doctest SPHINXOPTS="-q -n")
flake8 docs/ examples/ tests/ tools/ setup.py

# check if each header can be compiled on its own
#for f in include/gemmi/*.hpp; do gcc-7 -Ithird_party -c -fsyntax-only $f; done

# run tests under valgrind
#valgrind $BUILD_DIR/cpptest
#PYTHONMALLOC=malloc valgrind python3.7 -m unittest discover -s tests
