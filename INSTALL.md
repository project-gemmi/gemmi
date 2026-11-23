# Installing gemmi

The full installation guide - C++ library, Python module, conda, Fortran
and C bindings, WebAssembly, the `gemmi` program - is in the documentation:

  https://gemmi.readthedocs.io/en/latest/install.html

This file only adds a stop-gap for one case that the guide does not cover in
detail: a Python version for which no ready-made gemmi package exists yet.


## Python: no wheel for your Python version yet

`pip install gemmi` normally downloads a *wheel*: a ready-built package from
PyPI for one Python version and one platform (for example CPython 3.13 on
Linux x86-64). Installing a wheel needs no compiler and takes seconds.

When a new Python version comes out, wheels for it appear only with the next
gemmi release that is published on PyPI. Until then (for example gemmi 0.7.5
has wheels for Python 3.9 to 3.14, but not for 3.15), gemmi has to be built
from source on your machine. The commands below do that with pip; the result
is an ordinary installed gemmi module.

What you need:

* a C++ compiler: GCC or Clang on Linux, the Xcode command line tools on
  macOS (`xcode-select --install`), Visual Studio (Build Tools) on Windows;
* pip 25.3 or newer (for `--build-constraint` below):
  `python -m pip install --upgrade pip`.

pip fetches the other build tools (CMake, scikit-build-core, nanobind) by
itself. A build takes a few minutes; to let it use more cores, set
`CMAKE_BUILD_PARALLEL_LEVEL` (e.g. `export CMAKE_BUILD_PARALLEL_LEVEL=4`).

In the commands below, use the `python` of the environment you want gemmi in
(e.g. an activated virtual environment).

### A released version (e.g. 0.7.5)

    echo "nanobind<3" > build-constraints.txt
    python -m pip install --build-constraint build-constraints.txt \
        "gemmi @ git+https://github.com/project-gemmi/gemmi.git@v0.7.5"

The constraint is needed for 0.7.5 and earlier, which predate nanobind 3.
`python -m pip install --build-constraint build-constraints.txt gemmi==0.7.5`
(building from the source package on PyPI instead of the git tag) works too.

For a later release, change the tag (e.g. `@v0.7.6`). Releases made after
nanobind 3 support was added to gemmi do not need the constraint.

### The current development version (master)

    python -m pip install "gemmi @ git+https://github.com/project-gemmi/gemmi.git@master"

This gets the newest code, including changes that are not in any release yet.

### Build once, install many times

To install the same build into several environments (with the same Python
version and platform), build a wheel file first and install that:

    python -m pip wheel --no-deps -w wheelhouse \
        "gemmi @ git+https://github.com/project-gemmi/gemmi.git@master"
    python -m pip install wheelhouse/gemmi-*.whl

(For a release, use the tag and `--build-constraint` as above.) The wheel only
works for the Python version and platform it was built on; its file name says
which, e.g. `gemmi-0.7.6.dev0-cp315-cp315-linux_x86_64.whl`.

### Check the installation

    python -c "import gemmi; print(gemmi.__version__)"

Once a gemmi release with wheels for your Python version is on PyPI, a plain
`python -m pip install --upgrade gemmi` replaces the source build.
