#!/bin/bash

wd="$PWD"

cd .. # project root

if [[ ! -d m.css ]]; then
  # TODO: switch to official repo once PR#143 is accepted
  # git clone git://github.com/mosra/m.css
  git clone --depth 1  --branch default-argument-parsing --single-branch --origin sizmailov git://github.com/sizmailov/m.css
fi

cd m.css/documentation
python python.py "$wd/api/python/conf.py"

# dirty hack:
cd "$wd/_build/html/api/python/"
cp gemmi.html index.html

# Build C++ reference
cd "$wd/../m.css/documentation"
python doxygen.py "$wd/api/c++/Doxyfile-mcss"

# same dirty hack
cd "$wd/_build/html/api/c++/"
cp namespaces.html index.html

