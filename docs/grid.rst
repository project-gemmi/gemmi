
Volumetric grid
###############

3D data on an evenly spaced, rectangular grid.
The grid is aware of crystallographic symmetries,
but if symmetry is not set (or is set to P1)
it works as a box with periodic boundary conditions.

MRC/CCP4 maps
=============

Grid data can be read and written as MRC/CCP4 map file.
Gemmi handles the following modes of this format:

* 0 -- which correspond to C++ data type int8_t,
* 1 -- corresponds to int16_t,
* 2 -- float,
* and 6 -- uint16_t.

Mode 2 (float) is usually used for the electron density,
and mode 0 (int8_t) for masks, i.e. the 0/1 data that marks part of the volume
(e.g. the solvent region).

C++
===

::

    #include <gemmi/grid.hpp>

