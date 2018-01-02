
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
rows and columns that can correspond to any permutation of the X, Y and Z axes.
CCP4 files may contain only a part of the asymmetric unit,
or more than an asymmetric unit (i.e. redundant data).
There are two typical approaches to generate a crystallographic map:

* old-school way: a map covering a molecule with some margin
  around it is produced using CCP4 utilities such as ``fft`` and ``mapmask``,
* or a map is made for the asymmetric unit (asu), and the program that reads
  the map is supposed to expand the symmetry. This approach is used by
  the CCP4 clipper library and by programs that use this library, such as Coot.

The latter approach generates map for exactly one asu if possible,
i.e. if the shape of the asu, in fractional coordinates,
is rectangular. Otherwise, some redundancy cannot be avoided.

The maps generated for asu tend to be smaller than the maps around
the molecule (as compared in the
`UglyMol wiki <https://github.com/uglymol/uglymol/wiki/ccp4-dsn6-mtz>`_).

In crystallography CCP4 maps are rarely used nowadays, because most
of the programs can calculate the map on the fly from the reflection data.

C++
---

To read a map from a file::

    #include <gemmi/grid.hpp>

    gemmi::Grid<> grid;
    grid.read_ccp4_map(filename);

If the type of grid data differs from the type of data in file, the library
will attempt to convert the data when reading.

To work with the mask data (``int8_t``, but typically only values 0 and 1
are used) use::

    gemmi::Grid<int8_t> grid;

The CCP4 map header is organised as 56 words followed by space for ten
80-character text labels.
The Grid functions that access the data from the map header use the word
number (as in the format description) as a location in the header::

    int32_t header_i32(int w) const;
    float header_float(int w) const;
    // ccp4 map header has mostly 80-byte strings
    std::string header_str(int w, size_t len=80) const;

    void set_header_i32(int w, int32_t value);
    void set_header_float(int w, float value);
    void set_header_str(int w, const std::string& str);

For example::

    int mode = grid.header_i32(4);
    float x = grid.header_float(11);

``read_ccp4_map()`` stores the data as it is written in the file.
In many situation, it is convenient to have the data expanded to the whole
unit cell, with axes in a specific order (X, Y, Z is the most conventional
one). For this we have a function::

    grid.setup();

(Some of the functions described later in this section require this call.)


To write a map to a file::

    // if min/max/mean values in the map header are to be updated
    grid.hstats = gemmi::calculate_grid_statistics(grid.data);

    // if the header needs to be updated
    int mode = 2;
    grid.update_ccp4_header(mode);

    // finally, write the map to a file
    grid.write_ccp4_map(filename);

Data and symmetry
=================

The actual data is a just a custom 3d array with dimensions
``nu``, ``nv`` and ``nw``.
The data can be accessed in two ways::

    // quick: for 0<=u<nu, 0<=v<nv, 0<=w<nw.
    T get_value_q(int u, int v, int w) const;
    // safe: u, v, and w and wrapped using modulo function (u mod nu, etc.)
    T get_value_s(int u, int v, int w) const;

TODO

(how to set space group, unit cell, size)

Functions
=========

TODO: Higher-level functions. set_points_around()

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

Fortran
=======

TODO
