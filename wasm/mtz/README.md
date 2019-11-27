
# MTZ reader in WebAssembly

The MTZ file format is used for the storage of reflection data --
the data from diffraction experiments in macromolecular crystallography.

This library is intended primarily for web-based molecular viewers.
It can read MTZ file and transform map coefficients (reciprocal space)
to a density map (real space).

* It is part of GEMMI (license: MPL 2.0).
* It uses PocketFFT (license: MIT) for discrete Fourier transform.
* It uses Emscripten to compile C++ to WASM and JavaScript.

The API is yet to be documented.
