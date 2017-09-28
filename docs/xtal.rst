
Volumetric grid
###############

TODO: describe data structure

MRC/CCP4 maps
=============

Grid data can be read and written as MRC/CCP4 map file.
Gemmi handles modes 0, 1, 2 and 6 of this format,
which correspond to C++ data types int8_t, int16_t, float and uint16_t.
Mode 2 (float) is usually used for the electron density,
and mode 0 (int8_t) for masks, i.e. the 0/1 data that marks part of the volume
(e.g. the solvent region).

This part of the library is **not finished**.

C++
===

::

    #include <gemmi/grid.hpp>


Reflection files
################

(to be done)

SF mmCIF
========

MTZ
===

