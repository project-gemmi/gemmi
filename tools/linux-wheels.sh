#!/bin/bash
# based on
# https://github.com/pypa/python-manylinux-demo/blob/master/travis/build-wheels.sh
set -e -u -x


function repair_wheel {
    wheel="$1"
    if ! auditwheel show "$wheel"; then
        echo "Skipping non-platform wheel $wheel"
    else
        auditwheel repair "$wheel" --plat "$PLAT" -w /io/wheelhouse/
    fi
}

#yum install -y zlib-devel

# Compile wheels
#for PYBIN in /opt/python/cp38*/bin; do
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" install pybind11
    "${PYBIN}/pip" wheel -v /io/ --no-deps -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    repair_wheel "$whl"
done

# Install packages and test
for PYBIN in /opt/python/*/bin/; do
    "${PYBIN}/pip" install gemmi --no-index -f /io/wheelhouse
    "${PYBIN}/python" -m unittest discover -v -s /io/tests/
done
