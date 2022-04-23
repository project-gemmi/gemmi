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
date: 19 January 2022
bibliography: paper.bib

---

# Summary

GEMMI is a cross-platform library, accompanied by a set of small programs,
developed primarily for use in the field of macromolecular crystallography (MX).
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
and mmCIF formats and performing various operations on the reflections.

The third area is working with electron density maps – real or complex
values on a 3D grid. Electron density can be calculated from both the
structural model and experimental data. The functionality here includes
reading and writing files in the MRC/CCP4 map format, analysing and modifying
the density, and using the fast Fourier transform to switch between
the so-called direct space and the reciprocal space.

GEMMI is written in C++. It has Python bindings and, for selected functions,
also C and Fortran bindings. Interestingly, the library can be compiled to
WebAssembly for use in web applications. For example, UglyMol [@Wojdyr:2017]
uses it to read MTZ files inside the web browser.

# Statement of need

GEMMI is funded by two organizations that develop MX software:
CCP4 and Global Phasing Ltd. The aim is to deliver functionality
needed in other projects of these organizations.
Initially, the focus was on working with the PDBx/mmCIF file format,
then the scope was expanded to other areas.

The library has a significant overlap with other libraries used
in this field: CCTBX [@cctbx] and Clipper [@clipper].
But even when two functions from different libraries have similar purpose,
they usually differ in some aspects, for example, by making a different
trade-off between the speed of calculations and the accuracy of results,
or between the simplicity of the code and the number of provided options.

GEMMI is used in a number of projects, including
autoBUSTER [@autoBUSTER], CCP4i2 [@Potterton:2018],
CCP4 Cloud [@Krissinel:2018], Servalcat [@Yamashita:2021],
an analysis of covalent linkages [@Nicholls:2021],
reciprocalspaceship [@Greisman:2021], and many others.

# Acknowledgements

The library contains contributions from Keitaro Yamashita, Claus Flensburg
and other users. It uses third-party libraries:
[PocketFFT](https://gitlab.mpcdf.mpg.de/mtr/pocketfft) for Fast Fourier Transform,
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

This project would not be possible without Eugene Krissinel, Gérard Bricogne
and Garib Murshudov, who initiated it, and without many discussions with
users and with colleagues from Global Phasing and CCP4.

# References
