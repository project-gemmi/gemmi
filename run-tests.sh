#!/bin/sh

set -eu
cd "$(dirname "$0")"
(cd src && make -j4)
(cd python && make -j4)
(cd docs && make -j4 html SPHINXOPTS="-q -n")
PYTHONPATH=python python -m unittest discover -s tests
