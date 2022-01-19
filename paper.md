---
title: 'GEMMI: A library for structural biology'
tags:
  - structural biology
  - crystallography
  - CIF
authors:
  - name: Marcin Wojdyr
    orcid: 0000-0003-3980-4092
    affiliation: 1
affiliations:
 - name: Global Phasing Ltd., Cambridge, UK
   index: 1
date: 17 January 2022
bibliography: paper.bib

---

# Summary

GEMMI is a library, accompanied by a set of programs, developed primarily
for use in the field of macromolecular crystallography (MX).
Parts of this library are useful also in structural bioinformatics
and in chemical crystallography.

The library covers three main areas, which overlap and contain common elements,
such as handling of the crystallographic symmetry.

The first area is working with structural models of macromolecules.
This includes reading and writing files in the PDB and mmCIF formats,
analyzing and modifying models and working with restraint dictionaries.
The dictionaries are used to restrain geometry of a model
using prior knowledge about monomers in the model.

The second area is working with crystallographic data – experimentally
observed reflections. This includes reading and writing files in the MTZ
and mmCIF formats and performing commonly used operations on the reflections.

The third area is working with electron density maps – real or complex
values on a 3D grid. Electron density can be calculated from both the
structural model and experimental data. The functionality here includes
switching between the so-called direct space and the reciprocal space,
which involves Fourier transform.

GEMMI is written in C++. It has Python bindings and limited Fortran bindings.
Parts of the library have been compiled to WebAssembly for use in web
applications. For example, UglyMol [@Wojdyr:2017] uses it to read MTZ files.

# Statement of need

GEMMI is funded by CCP4 and Global Phasing Ltd., providers of MX software,
to deliver functionality needed in other projects of these organizations.
Initially, it was focused on working with the PDBx/mmCIF file format,
then it expanded to other areas.

The library has a significant overlap with two more mature libraries
in this field: CCTBX [@cctbx] and Clipper [@clipper].
Nevertheless, even when the same functionality is implemented again,
there is often a chance to make a different trade-off
between the speed of calculations and the accuracy of results,
or between the simplicity of the code and the number of provided options.

GEMMI is used in a number of projects, including
autoBUSTER [@autoBUSTER], CCP4i2 [@Potterton:2018],
CCP4 Cloud [@Krissinel:2018], Servalcat [@Yamashita:2021],
an analysis of covalent linkages in structural models [@Nicholls:2021],
reciprocalspaceship [@Greisman:2021], and many others.

# Acknowledgements

The library contains contributions from Keitaro Yamashita, Claus Flensburg
and other users. It uses third-party libraries:
PocketFFT [@pocketfft] for Fast Fourier Transform,
KSW2 [@ksw2] for sequence alignment,
QCProt [@qcp] for structure superposition,
Cromer-Liberman routine from Larch [@larch],
[PEGTL](https://github.com/taocpp/PEGTL) for creating PEG parsers,
as well as [sajson](https://github.com/chadaustin/sajson),
[stb_sprintf](https://github.com/nothings/stb),
[fast_float](https://github.com/fastfloat/fast_float),
[tinydir](https://github.com/cxong/tinydir),
[zlib](http://zlib.net/)
and [pybind11](https://github.com/pybind/pybind11).

The project would not be possible without Eugene Krissinel, Gerard Bricogne
and Garib Murshudov, who initiated it, and without many discussions with
users and with colleagues from Global Phasing and CCP4.

# References
