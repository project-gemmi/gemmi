
Volumetric grid
###############

When working with macromolecular models we often use
3D data on an evenly spaced, rectangular grid.
The data may represent electron density, a mask of the protein area,
or any other scalar value.

In Gemmi such a data is stored, together with metadata, in a class
called Grid. Actually, it is a set of classes for storing
different types of data: floating point numbers, integers, boolean masks.

The Grid class is aware of crystallographic symmetries,
but if symmetry is not set (or is set to P1)
it works as a box with periodic boundary conditions.

We support one file format for storing the grid data on disk: MRC/CCP4 map.
The Grid class contains low-level functions for accessing
and modifying the header (metadata) of this format.

MRC/CCP4 maps
=============

Grid data can be read and written as CCP4 map file.
Gemmi handles the following modes of this format:

* 0 -- which correspond to C++ data type int8_t,
* 1 -- corresponds to int16_t,
* 2 -- float,
* and 6 -- uint16_t.

Mode 2 (float) is usually used for the electron density,
and mode 0 (int8_t) for masks, i.e. the 0/1 data that marks part of the volume
(e.g. the solvent region).

The CCP4 file format is quite flexible. The data is stored as sections,
rows and columns that correspond to any permutation of the X, Y and Z axes.
Sometimes CCP4 files contain only a part of the asymmetric unit,
and sometimes they have redundant information for the sake of programs
that do not use symmetry to expand the data (but redundant data can be
necessary if the asymmetric unit is not a parallelepiped).

In Gemmi, we can either

* store the data read from a file as it is written in the file,
* or we can expand and transpose the data to cover whole unit cell
  with the rows, columns, sections corresponding to the X, Y, Z axes.

The latter is necessary to use functions that depend
on the symmetry of the unit cell.

Functions
=========

TODO

C++
===

In C++ all functionality related to grids is contained in a single header
file::

    #include <gemmi/grid.hpp>

There we have a templated ``struct Grid``::


    template<typename T=float> struct Grid;

TODO

Python
======

TODO
