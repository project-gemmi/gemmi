
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


Utilities
=========

gemmi-map
---------

(work in progress)

.. literalinclude:: map-help.txt
   :language: console


Reflection files
################

(to be done)

SF mmCIF
========

MTZ
===

Symmetry
########

Gemmi/symmetry provides space-group related functionality
needed in other parts of the library -- when working with coordinate
files, electron density maps and reflections.

Although the Gemmi project is developed for macromolecular crystallography
for which only 65 space groups are relevant,
we include all the 230 crystallographic space groups
for the sake of completeness.

This part of Gemmi has no dependencies:
it is all in a single C++ header :file:`symmetry.hpp`.

Space group notations
=====================

Gemmi tabulates 530 settings of the 230 crystallographic space groups,
including the following informations:

* space group numbers (1-230),
* ccp4 numbers (numbers assigned to particular settings; modulo 1000
  they give space group number: 3 and 1003 correspond to
  ``P 1 2 1`` and ``P 1 1 2``),
* Herman-Mauguin (H-M) symbols a.k.a. the international notation
  (``I a -3 d``, ``C 1 2 1``),
* extensions to the H-M notations (none, ``1``, ``2``, ``H`` or ``R``)
  that make extended H-M symbols (``R 3:H``, ``P 4/n:1``),
* short H-M symbols (``Ia-3d``, ``C2``),
* and the Hall notation (``-I 4bd 2c 3``, ``C 2y``, ``P 32 2c (0 0 -1)``).

Any of the above identifiers can be used to find a space group.

TODO: example in C++ or Python

The space group numbers, extended H-M symbols and Hall symbols are the
same as in ...

CCP4 numbers are taken from the CCP4 :file:`syminfo.lib` file.

International Tables for Crystallography (ITfC) describe
each of the 230 space groups in Volume A chapter 2.3.
Additionally, Vol.B ch.1.4 (in 2010 ed.) has two tables with computer-adapted
symbols:

* the Hall symbols (530 entries corresponding to different settings of
  the 230 groups),
* and the explicit symbols introduced by Shmueli_ (306 entries).

TODO: how it is related to other data tables: ITfC, syminfo.lib, sgtbx

Triplets and matrices
=====================

Crystallographic symmetry operations have a few notations.
Gemmi undertands only coordinate triplets (sometimes called
the Jones' faithful notation) such as ``x,x-y,z+1/2``.

The operations are equivalent to

* either a 3x3 rotation matrix and a translation
  vector,
* or to a single 4x4 transformation matrix (which in crystallography
  is called *Seitz matrix*).

TODO: example in Python

The triplets and matrices may not represent crystallographic symmetry,
but for example, generation of the biological assembly.
For this reason, ``x,y+1,z`` is not reduced to the identity.
But when the operations are transformed (two operations are combined,
or the inverse is calculated) they are assumed to be symmetry operations
and the shift is reduced.

(TODO: but do we need any functions working on non-symmetry operations,
such as combining operators or inversion?).

Space groups and operations
===========================

For a software developer, the choice between having an explicit list of all
the operations and generating them from a smaller set of symbols is a
trade-off between simplicity of the code and the amount of the tabulated data.

The number of symmetry operations per space group is between 1 and 192,
but in all the programs we've seen the centering vectors are kept separately,
which leaves up to 48 symmetry operations (point group m-3m)
combined with up to 4 centering vectors.

* The simplest way of storing the operations is to list them in a text file.

  * The most popular such files are CCP4 :file:`syminfo.lib` (540 entries)
    and its predecessor :file:`symop.lib`. :file:`syminfo.lib` is the primary
    spacegroup information for csymlib_ (part of libccp4, LGPL3)
    and a few other projects.
  * The same approach is used by OpenBabel_ (chemical toolbox in C++, GPL2):
    space-group names and triplets are tabulated in :file:`space-groups.txt`
    (currently 541 entries).

* The operations can be as well tabulated in the code.
  This approach is used (with one or two layers of indirection to reduce
  the data size) by spglib_ (C library, BSD) and NGL_.

* The inversion center (if applicable) can be kept separately,
  to reduce the maximum number of stored operations to 24+1.
  This is how the symmetry data is encoded in Fityk_.

* All the operations can be generated from only up-to-3 operations
  (+inversion), at the expense of more complex code.
  This approach is used in Mantid_ (triplets are tabulated in
  :file:`SpaceGroupFactory.cpp`, although the code contains up to 5
  generators per spacegroup).

* Finally, one can use one of the two computer-adapted descriptions from ITfC.
  The so-called explicit notation (``ICC$I3Q000$P4C393$P2D933``) is the
  longer of the two, but easier to parse by the computer.
  It is used in the SPGGEN_ program.

* The Hall notation (``-I 4bd 2c 3``), first proposed by Sydney R. Hall
  in 1981, is short (3-16 characters), unambiguous and looks nice.
  The notation can be interpreted by a few libraries:

  * SgInfo_ and SgLite_ (old C libraries from Ralf W. Grosse-Kunstleve
    recently re-licensed to BSD),
  * sgtbx_ (successor of SgInfo written in C++/Python, part of cctbx),
  * CCP4 Clipper_ (in :file:`spacegroup.cpp`).

  which makes it relatively popular.
  On the bad side, the conciseness is achieved by complex
  `rules <http://cci.lbl.gov/sginfo/hall_symbols.html>`_ of interpreting
  the symbols (including a few auxiliary tables), on the top of complexity
  needed to generate all operations from the encoded generators.
  (TODO: note about different versions of the Hall notation)

In Gemmi we derive operations from the Hall symbols
(despite of mixed feelings about the overly implicit notation).

.. _SgInfo: https://github.com/rwgk/sginfo
.. _SgLite: https://github.com/rwgk/sglite
.. _sgtbx: https://github.com/rwgk/sglite
.. _csymlib: http://www.ccp4.ac.uk/html/C_library/csymlib_8h.html
.. _spglib: https://atztogo.github.io/spglib/
.. _Clipper: http://www.ysbl.york.ac.uk/~cowtan/clipper/doc/
.. _OpenBabel: https://github.com/openbabel/openbabel
.. _Mantid: https://github.com/mantidproject/mantid
.. _Shmueli: http://dx.doi.org/10.1107/S0108767384001161
.. _NGL: https://github.com/arose/ngl
.. _Fityk: https://github.com/wojdyr/fityk
.. _SPGGEN: http://dx.doi.org/10.1107/S1600576716007330

C++
===

::

    #include <gemmi/symmetry.hpp>

