# MTZ reader in WebAssembly

The MTZ file format is used for the storage of reflection data --
the data from diffraction experiments in macromolecular crystallography.

This library is primarily intended for web-based molecular viewers.

It can read an MTZ file and transform map coefficients (reciprocal space)
to a density map (real space).

It uses:

* C++ code from [GEMMI](https://project-gemmi.github.io/) (license: MPL 2.0),
  for [handling MTZ files](https://gemmi.readthedocs.io/en/latest/hkl.html),
* [PocketFFT](https://gitlab.mpcdf.mpg.de/mtr/pocketfft) (license: MIT)
  for discrete Fourier transform,
* and [Emscripten](https://emscripten.org/) to compile C++ to WASM + JS.

You can see this package in action on https://uglymol.github.io/view/

Source code →
[gemmi repo](https://github.com/project-gemmi/gemmi/tree/master/wasm/mtz)

Binaries → `npm i mtz`


## Usage

The module is generated with Emscripten option MODULARIZE,
which wraps everything into a function that invoked (once) gives the module:

    var module = GemmiMtz();

To wait until the module is ready-to-use call a promise-like `then` method:

    GemmiMtz().then(function(module) {
      mtz = module.readMtz(array_buffer);
      ...
      mtz.delete();
    });

The `readMtz` function above takes ArrayBuffer with the content of an MTZ file
and returns Mtz object with:

* `cell` -- object representing unit cell with properties corresponding to
  cell parameters (`a`, `b`, `c`, `alpha`, `beta`, `gamma`),
* `calculate_map_from_labels(f_label: string, phi_label: string)` --
  calculates a map by transforming map coefficients from the requested
  columns; returns Float32Array (length == nx * ny * nz) representing
  density values on the grid,
* `calculate_map(diff_map: Boolean)` -- calculates a map with default labels:
  - FWT/PHWT or 2FOFCWT/PH2FOFCWT for normal map (`diff_map===false`),
  - DELFWT/PHDELWT or FOFCWT/PHFOFCWT for difference map,
* `last_error` -- contains error message if the functions above return null,
* `nx`, `ny`, `nz` -- dimension of the last calculated map
  (necessary to interpret the flat array as a 3D array),
* `rmsd` -- also a property of the last calculated map,
* `delete()` -- frees the allocated memory (inside the WebAssembly
  instance's memory).

`calculate_map` returns a map as a Float32Array view of the WASM memory.
This view is invalidated when the next map is calculated,
so you may want to make a copy.
