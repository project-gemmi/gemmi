
Volumetric grid
###############

When working with macromolecular models we often use
3D data on an evenly spaced, rectangular grid.
The data may represent the electron density, a mask of the protein area,
or any other scalar data.

In Gemmi such a data is stored in a class called Grid.
Actually, it is a set of classes for storing
different types of data: floating point numbers, integers or boolean masks.
These classes also store:

* unit cell dimensions (so the grid nodes can be assigned atomic coordinates),
* crystallographic symmetry (that determines which points on the grid
  are equivalent under the symmetry),
* and the header of the MRC/CCP4 file format if the data comes from a file.

If the symmetry is not set (or is set to P1)
then we effectively have a box with periodic boundary conditions.

MRC/CCP4 maps
=============

We support one file format for storing the grid data on disk: MRC/CCP4 map.
This format has a few different modes that correspond to different data types.
Gemmi supports:

* mode 0 -- which correspond to the C++ type int8_t,
* mode 1 -- corresponds to int16_t,
* mode 2 -- float,
* and mode 6 -- uint16_t.

CCP4 programs use mode 2 (float) for the electron density,
and mode 0 (int8_t) for masks. Mask is 0/1 data that marks part of the volume
(e.g. the solvent region). Other modes are not used in crystallography,
but may be used for CryoEM data.

This file format is quite flexible. The data is stored as sections,
rows and columns that correspond to a permutation of the X, Y and Z axes
as defined in the file header.
The file can contain only a part of the asymmetric unit,
or more than an asymmetric unit (i.e. redundant data).
There are two typical approaches to generate a crystallographic map:

* old-school way: a map covering a molecule with some margin
  around it is produced using CCP4 utilities such as ``fft`` and ``mapmask``,
* or a map is made for the asymmetric unit (asu), and the program that reads
  the map is supposed to expand the symmetry. This approach is used by
  the CCP4 clipper library and by programs that use this library, such as Coot.

The latter approach generates map for exactly one asu if possible,
i.e. if the shape of the asu in fractional coordinates
is rectangular. Otherwise, redundancy cannot be avoided.

The maps generated for asu tend to be smaller than the maps around
the molecule (as compared in the
`UglyMol wiki <https://github.com/uglymol/uglymol/wiki/ccp4-dsn6-mtz>`_).

Nowadays, the CCP4 format is rarely used in crystallography.
Almost all programs read the reflection data and calculate maps internally.

C++
---

To read a map from a file::

    #include <gemmi/grid.hpp>

    gemmi::Grid<float> grid;
    grid.read_ccp4_map(filename);

To work with the mask data (``int8_t``, but typically only values 0 and 1
are used) use ``int8_t`` instead of ``float``::

    gemmi::Grid<int8_t> grid;

If the grid data type does not match the file data type, the library
will attempt to convert the data when reading.

The CCP4 map header is organised as 56 words followed by space for ten
80-character text labels.
The member functions that access the data from the map header use the word
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

    grid.setup(GridSetup::Full, NAN);  // unknown values are set to NAN

(Some of the functions described later in this section require this call.)


To write a map to a file::

    // the file header needs to be prepared/updated with an explicit call
    int mode = 2; // ccp4 file mode: 2 for floating-point data, 0 for masks
    bool update_stats = true; // update min/max/mean/rms values in the header
    grid.update_ccp4_header(mode, update_stats);

    grid.write_ccp4_map(filename);

Python
------

The Python API is similar.

.. doctest::

    >>> import gemmi
    >>>
    >>> m = gemmi.read_ccp4_map('../tests/5i55_tiny.ccp4')
    >>> m.nu, m.nv, m.nw  # small numbers as it is a toy example
    (8, 6, 10)
    >>> m.space_group
    <gemmi.SpaceGroup("P 1 21 1")>
    >>> m.unit_cell
    <gemmi.UnitCell(29.45, 10.5, 29.7, 90, 111.975, 90)>
    >>> m.setup()
    >>> m.nu, m.nv, m.nw
    (60, 24, 60)

The low-level header access has three getters and three setters,
as in the C++ version.

.. doctest::

    >>> m.header_float(20), m.header_float(21)  # dmin, dmax
    (-0.5310382843017578, 2.3988280296325684)
    >>> m.header_i32(28)
    0
    >>> m.set_header_i32(28, 20140)
    >>> m.header_str(57, 80).strip()
    'Created by MAPMAN V. 080625/7.8.5 at Wed Jan 3 12:57:38 2018 for A. Nonymous'

Data and symmetry
=================

The actual data is a 3d array with dimensions ``nu``, ``nv`` and ``nw``,
internally kept in a grid member ``std::vector<T> data``.
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
